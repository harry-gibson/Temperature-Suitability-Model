using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace TempSuitability_CSharp
{
    /// <summary>
    /// A Cohort is the smallest tracked unit of mosquito population. Each cohort tracks 
    /// its infectiousness and surviving proportion of mosquitoes, based on the temperatures 
    /// it has been subjected to.
    /// </summary>
    class Cohort
    {
        #region fields
        private double m_Quantity;
        private double m_degDays;
        private int m_Age;

        private readonly double m_InfectiousThreshold;
        private readonly int m_LifespanSlices;

        #endregion
        #region properties
        private double Contribution
        {
            get
            {
                if (m_degDays <= m_InfectiousThreshold || IsDead)
                {
                    return 0;
                }
                return m_Quantity;
            }
        }
        public bool IsDead
        {
            get
            {
                return m_Age >= m_LifespanSlices;
            }
        }
        #endregion

        #region construct new cohort
        public Cohort(double InfectiousnessThreshold, int LifespanSlices)
        {
            m_Quantity = 1;
            m_degDays = 0;
            m_InfectiousThreshold = InfectiousnessThreshold;
            m_LifespanSlices = LifespanSlices;
            m_Age = 0;
        }
        #endregion

        #region cohort's lifecycle
        private void Decay(double survivalRate)
        {
            //double survivalRate = 
            //    Math.Pow(
            //        (Math.Exp(-1 / (-4.4 + (1.31 * currentTemp) - (0.03 * (Math.Pow(currentTemp, 2)))))), 
            //        (sliceLengthDays));
            m_Quantity = m_Quantity * survivalRate;
        }

        public double Iterate(double degDays, double survivalRate)
        {
            // We take degree days and temperature as params rather than calculating one from 
            // the other to avoid doing this calculation for every living cohort - just do it once 
            // in the population class instead.
            // The temperature suitability contribution from this cohort is based on what's left 
            // at the beginning of the timeslice, not what's left after the dying that takes 
            // place during it
            double prevQuantity = Contribution;
            Decay(survivalRate);
            m_degDays += degDays;
            m_Age += 1;
            return prevQuantity;           
        }
        #endregion
    }

    /// <summary>
    /// Represents the population of mosquitoes living at a location, comprising a number of 
    /// Cohorts - one for each timeslice in the lifespan of a mosquito.
    /// </summary>
    class Population
    {
        /// <summary>
        /// Population is implemented as a Queue (FIFO) of Cohorts. It's private so we can ensure 
        /// that it remains of fixed size by controlling access to it.
        /// </summary>
        private Queue<Cohort> m_Population;

        private readonly int m_LifespanSlices;
        private readonly double m_SliceLengthDays;
        private readonly double m_TempThreshold;
        private readonly double m_InfectionThreshold;
        public bool IsInitialised { get; private set; }
        public Population(int LifespanSlices, double SliceLengthDays, double TempThreshold, double InfectionThreshold)
        {
            m_LifespanSlices = LifespanSlices;
            m_SliceLengthDays = SliceLengthDays;
            m_InfectionThreshold = InfectionThreshold;
            m_TempThreshold = TempThreshold;
            IsInitialised = false;
            // Create the population "queue" with length of one lifespan 
            m_Population = new Queue<Cohort>(LifespanSlices);
        }

        /// <summary>
        /// Initialise the population with an initial set of cohorts, iterating each one 
        /// every time a new one is added so they are all appropriately "decayed"
        /// </summary>
        /// <param name="SpinUpTemps"></param>
        public void Initialise(double[] SpinUpTemps)
        {
            if (SpinUpTemps.Length != m_LifespanSlices)
            {
                throw new ArgumentException("Population must be initialised with one temp for each slice in a lifespan");
            }
            for (int i = 0; i < m_LifespanSlices; i++)
            {
                double minTemp = SpinUpTemps[i];
                double survivalRate = GetSurvivingFraction(minTemp);
                double degreeDays = Math.Max(((minTemp - m_TempThreshold) * m_SliceLengthDays), 0);
                
                // we don't care what the actual result of the sum is at this point
                // It's much neater to say:
                //m_Population.Sum(coh => coh.Iterate(degreeDays, survivalRate));
                // But also much slower! Use ye olde loop instead.
                foreach (Cohort coh in m_Population)
                {
                    coh.Iterate(degreeDays, survivalRate);
                }
                m_Population.Enqueue(new Cohort(m_InfectionThreshold, m_LifespanSlices));
            }
            IsInitialised = true;
        }

        /// <summary>
        /// Move the model forward by one timeslice. 
        /// 
        /// Applies the minimum temperature to all living cohorts, and returns the summarised TS contribution 
        /// of those cohorts, calculated before 
        /// the cohorts are themselves iterated to decay according to the temperature they're subjected 
        /// to. Dead cohort is removed and a new cohort is hatched.
        /// </summary>
        /// <param name="minTemp"></param>
        /// <returns></returns>
        public double Iterate(double minTemp)
        {
            if (!IsInitialised)
            {
                throw new InvalidOperationException("Population has not yet been initialised");
            }
            // Calculate degree days once here rather than in every Cohort.
            double degreeDays = Math.Max(((minTemp - m_TempThreshold) * m_SliceLengthDays), 0);
            // Calculate survival fraction once here rather than in every Cohort
            double survivalRate = GetSurvivingFraction(minTemp);
            // Temperature suitability at this slice is simply the sum of Contributions from every Cohort.
            // This is given before applying the death / decay of this timeslice.
            // So all we have to do is this time-slice's temperature to all currently living cohorts and 
            // summarise the return values.
            double tsAtSlice = 0;
            // It would be sweeter to say:
            // double tsAtSLice = m_Population.Sum(coh => coh.Iterate(degreeDays, survivalRate));
            // But that's substantially slower thanks to the linq overhead, given the tightness of this loop.
            foreach (Cohort coh in m_Population)
            {
                tsAtSlice += coh.Iterate(degreeDays, survivalRate);
            }

            // Remove the oldest cohort, which should now be dead anyway
            Cohort tDeadCohort = m_Population.Dequeue();
            Debug.Assert(tDeadCohort.IsDead);
            // Add the new cohort
            m_Population.Enqueue(new Cohort(m_InfectionThreshold, m_LifespanSlices));
            return tsAtSlice;
        }

        /// <summary>
        /// The fraction surviving is the same for all cohorts as it only depends on temperature, so we 
        /// calculate it once at each timeslice for the whole population
        /// </summary>
        /// <param name="minTemp"></param>
        /// <returns></returns>
        private double GetSurvivingFraction(double minTemp)
        {
            return Math.Pow(
              (Math.Exp(-1 / (-4.4 + (1.31 * minTemp) - (0.03 * (Math.Pow(minTemp, 2)))))),
              (m_SliceLengthDays));
        }
    }
}
