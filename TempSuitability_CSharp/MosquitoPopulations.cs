using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Runtime.CompilerServices;
namespace TempSuitability_CSharp
{
    /// <summary>
    /// Classes implementing this interface represent a population of mosquitoes at a cell,
    ///  which comprises a number of cohorts (the mosquitoes born in a single timeslice) - 
    /// one cohort for each timeslice in the lifespan of a mosquito. The population can 
    /// be iterated (moved forward by a timestep) to run the temperature suitability model. 
    /// Interface allows for different underlying implementations of the population to test performance 
    /// of different approaches.
    /// </summary>
    interface IIterablePopulation
    {
        bool IsInitialised { get; }
        void Initialise(double[] DailyLifetimeTemps);
        double Iterate(double DailyMinTemp);
    }

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
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
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
    /// Implements a mosquito population at a cell using a normal OO approach based on 
    /// C# collection types (Queue) of a new type Cohort. 
    /// This is quite slow so not used in production but the code was written this way originally.
    /// </summary>
    class PopulationOO : IIterablePopulation
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
        private readonly double m_UpperTempLimit;

        public bool IsInitialised { get; private set; }

        public PopulationOO(int LifespanSlices, double SliceLengthDays, double TempThreshold, double InfectionThreshold, double DeathAboveTemperature)
        {
            m_LifespanSlices = LifespanSlices;
            m_SliceLengthDays = SliceLengthDays;
            m_InfectionThreshold = InfectionThreshold;
            m_TempThreshold = TempThreshold;
            m_UpperTempLimit = DeathAboveTemperature;
            IsInitialised = false;
            // Create the population "queue" with length of one lifespan 
            m_Population = new Queue<Cohort>(LifespanSlices);
        }

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
                double survivalRate = MossieMethods.GetSurvivingFraction(minTemp, m_SliceLengthDays, m_UpperTempLimit);
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
        /// Move the population forward by one timeslice. 
        /// 
        /// Applies the minimum temperature to all living cohorts, and returns the summarised TS contribution 
        /// of those cohorts, calculated before 
        /// the cohorts are themselves iterated to decay according to the temperature they're subjected 
        /// to. Dead cohort is removed and a new cohort is hatched.
        /// </summary>
        /// <param name="minTemp"></param>
        /// <returns></returns>
        /// 
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double Iterate(double minTemp)
        {
            if (!IsInitialised)
            {
                throw new InvalidOperationException("Population has not yet been initialised");
            }
            // Calculate degree days once here rather than in every Cohort.
            double degreeDays = Math.Max(((minTemp - m_TempThreshold) * m_SliceLengthDays), 0);
            // Calculate survival fraction once here rather than in every Cohort
            double survivalRate = MossieMethods.GetSurvivingFraction(minTemp, m_SliceLengthDays, m_UpperTempLimit);
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
            //Debug.Assert(tDeadCohort.IsDead);
            // Add the new cohort
            m_Population.Enqueue(new Cohort(m_InfectionThreshold, m_LifespanSlices));
            return tsAtSlice;

        }
    }

    /// <summary>
    /// Implements a mosquito population at a cell using basic array types, as opposed to Queue and Cohort, 
    /// for improved performance.
    /// The population is simply represented by two parallel arrays,
    /// one to track the contribution of all living cohorts and one to track the number of degree-days they
    /// have each been exposed to. We maintain a property to track which is the "oldest" which moves forwards 
    /// at each iteration
    /// </summary>
    class PopulationArray :IIterablePopulation
    {
        private double[] m_Contribs;
        private double[] m_DegDays;
        private int m_NextToDie = 0;

        private readonly int m_LifespanSlices;
        private readonly double m_SliceLengthDays;
        private readonly double m_TempThreshold;
        private readonly double m_InfectionThreshold;
        private readonly double m_UpperTempLimit;

        public bool IsInitialised { get; private set; }
        public PopulationArray(int LifespanSlices, double SliceLengthDays, double TempThreshold, double InfectionThreshold, double DeathAboveTemperature)
        {
            m_LifespanSlices = LifespanSlices;
            m_SliceLengthDays = SliceLengthDays;
            m_InfectionThreshold = InfectionThreshold;
            m_TempThreshold = TempThreshold;
            m_UpperTempLimit = DeathAboveTemperature;
            IsInitialised = false;
            // Create the population "queue" with length of one lifespan 
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

            m_Contribs = new double[m_LifespanSlices];
            m_DegDays = new double[m_LifespanSlices];
            double[] localContribs = m_Contribs;
            double[] localDegDays = m_DegDays;

            for (int i = 0; i < m_LifespanSlices; i++)
            {
                double minTemp = SpinUpTemps[i];
                // survival rate is zero if the upper temp threshold is exceeded, meaning all cohorts will die
                double survivalRate = MossieMethods.GetSurvivingFraction(minTemp, m_SliceLengthDays, m_UpperTempLimit);
                double degreeDays = Math.Max(((minTemp - m_TempThreshold) * m_SliceLengthDays), 0);

                // we don't care what the actual result of the sum is at this point
                for (int cohPos = 0; cohPos < i; cohPos++)
                {
                    localContribs[cohPos] *= survivalRate;
                    localDegDays[cohPos] += degreeDays;
                }
                localContribs[i] = 1;
                localDegDays[i] = 0;
            }
            m_Contribs = localContribs;
            m_DegDays = localDegDays;
            IsInitialised = true;
        }

        /// <summary>
        /// Move the population forward by one timeslice. 
        /// 
        /// Applies the minimum temperature to all living cohorts, and returns the summarised TS contribution 
        /// of those cohorts, calculated before 
        /// the cohorts are themselves iterated to decay according to the temperature they're subjected 
        /// to. Dead cohort is removed and a new cohort is hatched.
        /// </summary>
        /// <param name="minTemp"></param>
        /// <returns></returns>
        /// 
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double Iterate(double minTemp)
        {
            if (!IsInitialised)
            {
                throw new InvalidOperationException("Population has not yet been initialised");
            }
            // work on a local copy rather than referencing the fields in the tight loop;
            // this saves approx 10% of overall program runtime!
            double[] localContribs = m_Contribs;
            double[] localDegDays = m_DegDays;
            double localThreshold = m_InfectionThreshold;

            // Calculate degree days once here rather than in every Cohort.
            double degreeDays = Math.Max(((minTemp - m_TempThreshold) * m_SliceLengthDays), 0);
            // Calculate survival fraction once here rather than in every Cohort
            double survivalRate = MossieMethods.GetSurvivingFraction(minTemp, m_SliceLengthDays, m_UpperTempLimit);
            // Temperature suitability at this slice is simply the sum of Contributions from every Cohort.
            // This is given before applying the death / decay of this timeslice.
            // So all we have to do is apply this time-slice's temperature to all currently living cohorts and 
            // summarise the return values.
            double tsAtSlice = 0;
            // calculate and summarise the contribution of all living cohorts first...
            for (int cohPos = 0; cohPos < localContribs.Length; cohPos++)
            {
                if (localDegDays[cohPos] > localThreshold)
                {
                    tsAtSlice += localContribs[cohPos];
                }
                // ... BEFORE applying the decay function
                localContribs[cohPos] *= survivalRate;
                localDegDays[cohPos] += degreeDays;
            }
            // ... and BEFORE replacing the oldest dead cohort with a new one
            localContribs[m_NextToDie] = 1;
            localDegDays[m_NextToDie] = 0;

            m_Contribs = localContribs;
            m_DegDays = localDegDays;
            m_NextToDie = (m_NextToDie + 1) % m_LifespanSlices;

            return tsAtSlice;
        }
    }

    /// <summary>
    /// Implements a mosquito population at a cell using pointer types, as opposed to Queue and Cohort, 
    /// or native array types.
    /// Due to the tightness of the Iterate loop this is substantially the fastest implementation.
    /// (The iterate loop runs 372 * 12 * 365 * 17 times at every pixel location!) 
    /// The population is simply represented by two parallel arrays (which we will access via pointers),
    /// one to track the contribution of all living cohorts and one to track the number of degree-days they
    /// have each been exposed to. We maintain a property to track which is the "oldest" which moves forwards 
    /// at each iteration
    /// </summary>
    class PopulationPtr : IIterablePopulation
    {
        private double[] m_Contribs;
        private double[] m_DegDays;
        private int m_NextToDie = 0;

        private readonly int m_LifespanSlices;
        private readonly double m_SliceLengthDays;
        private readonly double m_TempThreshold;
        private readonly double m_InfectionThreshold;
        private readonly double m_UpperTempLimit;

        public bool IsInitialised { get; private set; }
        public PopulationPtr(int LifespanSlices, double SliceLengthDays, double TempThreshold, double InfectionThreshold, double DeathAboveTemperature)
        {
            m_LifespanSlices = LifespanSlices;
            m_SliceLengthDays = SliceLengthDays;
            m_InfectionThreshold = InfectionThreshold;
            m_TempThreshold = TempThreshold;
            m_UpperTempLimit = DeathAboveTemperature;
            IsInitialised = false;
        }

        /// <summary>
        /// Initialise the population with an initial set of cohorts, iterating each one 
        /// every time a new one is added so they are all appropriately "decayed"
        /// </summary>
        /// <param name="SpinUpTemps"></param>
        public unsafe void Initialise(double[] SpinUpTemps)
        {
            if (SpinUpTemps.Length != m_LifespanSlices)
            {
                throw new ArgumentException("Population must be initialised with one temp for each slice in a lifespan");
            }
            // initialise the population arrays (all values will be set to 0)
            m_Contribs = new double[m_LifespanSlices];
            m_DegDays = new double[m_LifespanSlices];
            double[] localContribs = m_Contribs;
            double[] localDegDays = m_DegDays;

            for (int i = 0; i < m_LifespanSlices; i++)
            {
                double minTemp = SpinUpTemps[i];
                double survivalRate = MossieMethods.GetSurvivingFraction(minTemp, m_SliceLengthDays, m_UpperTempLimit);
                double degreeDays = Math.Max(((minTemp - m_TempThreshold) * m_SliceLengthDays), 0);

                // decrease all the living cohorts by the survival fraction and 
                // we don't care what the actual result of the sum is at this point
                for (int cohPos = 0; cohPos < i; cohPos++)
                {
                    localContribs[cohPos] *= survivalRate;
                    localDegDays[cohPos] += degreeDays;
                }
                // "hatch" the next cohort
                localContribs[i] = 1;
                localDegDays[i] = 0;
            }
            m_Contribs = localContribs;
            m_DegDays = localDegDays;
            IsInitialised = true;
        }

        /// <summary>
        /// Move the population forward by one timeslice. 
        /// 
        /// Applies the minimum temperature to all living cohorts, and returns the summarised TS contribution 
        /// of those cohorts, calculated BEFORE the cohorts are themselves iterated to decay according to 
        /// the temperature they're subjected to. 
        /// The oldest (dead) cohort is then removed and a new cohort is hatched.
        /// </summary>
        /// <param name="SliceTemp"></param>
        /// <returns></returns>
        /// 
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe double Iterate(double SliceTemp)
        {
            if (!IsInitialised)
            {
                throw new InvalidOperationException("Population has not yet been initialised");
            }
            // Calculate degree days once here rather than in every Cohort.
            double degreeDays = (SliceTemp - m_TempThreshold) * m_SliceLengthDays;
            if (degreeDays < 0)
            {
                degreeDays = 0;
            }
            // Calculate survival fraction once here rather than in every Cohort. It may be interesting at some stage 
            // to investigate calculating a survival rate that varies with age to reflect different acceptable temp 
            // ranges at different stages of the lifecycle.
            double survivalRate = MossieMethods.GetSurvivingFraction(SliceTemp, m_SliceLengthDays, m_UpperTempLimit);
            // Temperature suitability at this slice is simply the sum of Contributions from every Cohort.
            // This is given before applying the death / decay of this timeslice.
            // So all we have to do is apply this time-slice's temperature to all currently living cohorts and 
            // summarise the return values.
            double tsAtSlice = 0;

            // work on a local copy rather than referencing the fields in the tight loop;
            // this single step saves approx 10% of overall program runtime in the array-based implementation.
            double[] localContribs = m_Contribs;
            double[] localDegDays = m_DegDays;
            double localThreshold = m_InfectionThreshold;
            int count = m_LifespanSlices;
            fixed (double* pContribs = localContribs, pDegDays = localDegDays)
            {
                double* pContrib = pContribs, pDegDay = pDegDays;
                for (int i = 0; i < count; i++)
                {
                    // go through each cohort of the current population, based on pointers
                    if (*pDegDay > localThreshold)
                    {
                        // this cohort has reached sporogenesis threshold, it contributes to the temperature suitability 
                        // of this timeslice based on the surviving fraction _before_ this timeslice's die-off
                        tsAtSlice += *pContrib;
                    }
                    // reduce the surviving fraction of this cohort (possibly to zero if mosoquito-baking temp was exceeded)
                    *pContrib *= survivalRate;
                    // add the degree-days of this slice to the total degree days experienced by this cohort (probably quicker 
                    // just to do it rather than only conditionally doing it if baking didn't occur)
                    *pDegDay += degreeDays;
                    // move pointer for contribution and degree day onto next value
                    pContrib++;
                    pDegDay++;
                }
            }
            // spawn a new one (replacing the previous value at this position, the cohort which is now dead)
            localContribs[m_NextToDie] = 1;
            localDegDays[m_NextToDie] = 0;
            // copy local variables back to field
            m_Contribs = localContribs;
            m_DegDays = localDegDays;
            m_NextToDie = (m_NextToDie + 1) % m_LifespanSlices;
            return tsAtSlice;
        }

    }

    static class MossieMethods
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double GetSurvivingFraction(double sliceTemp, double sliceLengthDays, double deathTemp)
        {
            if (sliceTemp > deathTemp)
            {
                return 0;
            }

            var f = Math.Pow(
             (Math.Exp(-1 / (-4.4 + (1.31 * sliceTemp) - (0.03 * (Math.Pow(sliceTemp, 2)))))),
             (sliceLengthDays));
            return f < 0 ? 0 : f >= 1 ? 0 : f;
        }
    }
}
