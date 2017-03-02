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
    /// which comprises a number of cohorts (the mosquitoes born in a single timeslice) - 
    /// one cohort for each timeslice in the lifespan of a mosquito. The population can 
    /// be iterated (moved forward by a timestep) by supplying a temperature value which will 
    /// be used to calculate the temperature-dependent survival of the population and thus 
    /// to run the temperature suitability model.
    /// Interface allows for different underlying implementations of the population to test performance 
    /// of different approaches.
    /// </summary>
    interface ITemperatureIterablePopulation
    {
        bool IsInitialised { get; }
        void Initialise(double[] DailyLifetimeTemps);
        double Iterate(double DailyMinTemp);
    }

    /// <summary>
    /// Classes implementing this interface represent a population of mosquitoes at a cell, 
    /// which comprises a number of cohorts (the mosquitoes born in a single timeslice) - 
    /// one cohort for each timeslice in the lifespan of a mosquito. The population can 
    /// be iterated (moved forward by a timestep) by supplying a temperature and a relative 
    /// humidity value which will be used to calculate the temperature and humidity dependent 
    /// survival of the population and thus to run the temperature suitability model.
    /// Interface allows for different underlying implementations of the population to test performance 
    /// of different approaches.
    /// </summary>
    interface ITempHumidityIterablePopulation :ITemperatureIterablePopulation
    {
        void Initialise(double[] DailyLifetimeTemps, double[] DailyLifetimeRH);
        double Iterate(double DailyMinTemp, double DailyRHPercent);
        bool GetRequiresHumidity();
    }

    /// <summary>
    /// enum representing the classes currently known to implement ITemperatureIterablePopulation
    /// </summary>
    enum PopulationTypes
    {
        /// <summary>
        /// A population that is represented by a linq Queue structure: most C#-ish but slow
        /// </summary>
        OOCohorts,
        /// <summary>
        /// A population that is represented by a series of arrays
        /// </summary>
        Arrays,
        /// <summary>
        /// A population that is represented by a series of arrays with pointer access
        /// </summary>
        Pointers
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
    class PopulationOO : ITemperatureIterablePopulation
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
        
        public PopulationOO(PopulationParams modelConfigParams)
        {
            m_LifespanSlices = modelConfigParams.LifespanSlices;
            m_SliceLengthDays = modelConfigParams.SliceLength.TotalDays;
            m_InfectionThreshold = modelConfigParams.DegreeDayThreshold;
            m_TempThreshold = modelConfigParams.MinTempThreshold;
            m_UpperTempLimit = modelConfigParams.MosquitoDeathTemperature;
            IsInitialised = false;
            m_Population = new Queue<Cohort>(m_LifespanSlices);
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
                double survivalRate = MossieMethods.GetSurvivingFraction_Martens2(minTemp, m_SliceLengthDays, m_UpperTempLimit);
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
            double survivalRate = MossieMethods.GetSurvivingFraction_Martens2(minTemp, m_SliceLengthDays, m_UpperTempLimit);
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
    class PopulationArray :ITemperatureIterablePopulation
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
    
        public PopulationArray(PopulationParams modelConfigParams)
        {
            m_LifespanSlices = modelConfigParams.LifespanSlices;
            m_SliceLengthDays = modelConfigParams.SliceLength.TotalDays;
            m_InfectionThreshold = modelConfigParams.DegreeDayThreshold;
            m_TempThreshold = modelConfigParams.MinTempThreshold;
            m_UpperTempLimit = modelConfigParams.MosquitoDeathTemperature;
            IsInitialised = false;
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
                double survivalRate = MossieMethods.GetSurvivingFraction_Martens2(minTemp, m_SliceLengthDays, m_UpperTempLimit);
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
            double survivalRate = MossieMethods.GetSurvivingFraction_Martens2(minTemp, m_SliceLengthDays, m_UpperTempLimit);
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
    class PopulationPtr : ITemperatureIterablePopulation
    {
        protected double[] m_Contribs;
        protected double[] m_DegDays;
        protected int m_NextToDie = 0;
        protected MossieMethods.SurvivalFunctions m_survivalFunc;

        protected readonly int m_LifespanSlices;
        protected readonly double m_SliceLengthDays;
        protected readonly double m_TempThreshold;
        protected readonly double m_InfectionThreshold;
        protected readonly double m_UpperTempLimit;
        
        public bool IsInitialised { get; protected set; }

        /// <summary>
        /// The survival function for this implementation of the population is fixed as the "Martens2" function
        /// </summary>
        /// <param name="modelConfigParams"></param>
        public PopulationPtr(PopulationParams modelConfigParams)
        {
            m_LifespanSlices = modelConfigParams.LifespanSlices;
            m_SliceLengthDays = modelConfigParams.SliceLength.TotalDays;
            m_InfectionThreshold = modelConfigParams.DegreeDayThreshold;
            m_TempThreshold = modelConfigParams.MinTempThreshold;
            m_UpperTempLimit = modelConfigParams.MosquitoDeathTemperature;
            IsInitialised = false;
            m_survivalFunc = MossieMethods.SurvivalFunctions.Martens2;
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
                double survivalRate = MossieMethods.GetSurvivingFraction_Martens2(minTemp, m_SliceLengthDays, m_UpperTempLimit);
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
            double survivalRate = MossieMethods.GetSurvivingFraction_Martens2(SliceTemp, m_SliceLengthDays, m_UpperTempLimit);
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

    /// <summary>
    /// Implements a mosquito population at a cell which can be iterated with a temperature+humidity pair of values
    /// to decay according to a choice of different survival functions.
    /// </summary>
    class ConfigurableDecayPopulation: PopulationPtr, ITempHumidityIterablePopulation
    {
        public ConfigurableDecayPopulation(PopulationParams modelConfigParams) : base(modelConfigParams)
        {
            m_survivalFunc = modelConfigParams.SurvivalFunction;
        }
        
        public bool GetRequiresHumidity()
        {
              switch (m_survivalFunc)
                {
                    case MossieMethods.SurvivalFunctions.Martens2:
                    case MossieMethods.SurvivalFunctions.Martens3:
                    case MossieMethods.SurvivalFunctions.BayohMordecai:
                    case MossieMethods.SurvivalFunctions.BayohErmert:
                        return false;
                    default:
                        return true;
                }
        }

        public unsafe void Initialise(double[] SpinUpTemps, double[] SpinUpHumidities)
        {
            if (SpinUpTemps.Length != m_LifespanSlices || SpinUpHumidities.Length != m_LifespanSlices)
            {
                throw new ArgumentException("Population must be initialised with one temp and humidity for each slice in a lifespan");
            }
            // initialise the population arrays (all values will be set to 0)
            m_Contribs = new double[m_LifespanSlices];
            m_DegDays = new double[m_LifespanSlices];
            double[] localContribs = m_Contribs;
            double[] localDegDays = m_DegDays;

            for (int i = 0; i < m_LifespanSlices; i++)
            {
                double sliceTemp = SpinUpTemps[i];
                double sliceHum = SpinUpHumidities[i];
                double survivalRate;
                switch (m_survivalFunc)
                {
                    case MossieMethods.SurvivalFunctions.BayohErmert:
                        survivalRate = MossieMethods.GetSurvivingFraction_BayohErmert(sliceTemp, m_SliceLengthDays, m_UpperTempLimit);
                        break;
                    case MossieMethods.SurvivalFunctions.BayohMordecai:
                        survivalRate = MossieMethods.GetSurvivingFraction_BayohMordecai(sliceTemp, m_SliceLengthDays, m_UpperTempLimit);
                        break;
                    case MossieMethods.SurvivalFunctions.Martens2:
                        survivalRate = MossieMethods.GetSurvivingFraction_Martens2(sliceTemp, m_SliceLengthDays, m_UpperTempLimit);
                        break;
                    case MossieMethods.SurvivalFunctions.Martens3:
                        survivalRate = MossieMethods.GetSurvivingFraction_Martens3(sliceTemp, m_SliceLengthDays, m_UpperTempLimit);
                        break;
                    case MossieMethods.SurvivalFunctions.BayohParham:
                        survivalRate = MossieMethods.GetSurvivingFraction_BayohParham(sliceTemp, sliceHum, m_SliceLengthDays, m_UpperTempLimit);
                        break;
                    case MossieMethods.SurvivalFunctions.BayohLunde:
                        survivalRate = MossieMethods.GetSurvivingFraction_BayohLunde(sliceTemp, sliceHum, m_SliceLengthDays, m_UpperTempLimit);
                        break;
                    default:
                        throw new InvalidOperationException("Unrecognised survival function / not implemented");
                }
                double degreeDays = Math.Max(((sliceTemp - m_TempThreshold) * m_SliceLengthDays), 0);

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

        public unsafe double Iterate (double SliceTemp, double SliceHumidity)
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
            double survivalRate;// = MossieMethods.GetSurvivingFraction_Martens2(SliceTemp, m_SliceLengthDays, m_UpperTempLimit);
            switch (m_survivalFunc)
            {
                case MossieMethods.SurvivalFunctions.BayohErmert:
                    survivalRate = MossieMethods.GetSurvivingFraction_BayohErmert(SliceTemp, m_SliceLengthDays, m_UpperTempLimit);
                    break;
                case MossieMethods.SurvivalFunctions.BayohMordecai:
                    survivalRate = MossieMethods.GetSurvivingFraction_BayohMordecai(SliceTemp, m_SliceLengthDays, m_UpperTempLimit);
                    break;
                case MossieMethods.SurvivalFunctions.Martens2:
                    survivalRate = MossieMethods.GetSurvivingFraction_Martens2(SliceTemp, m_SliceLengthDays, m_UpperTempLimit);
                    break;
                case MossieMethods.SurvivalFunctions.Martens3:
                    survivalRate = MossieMethods.GetSurvivingFraction_Martens3(SliceTemp, m_SliceLengthDays, m_UpperTempLimit);
                    break;
                case MossieMethods.SurvivalFunctions.BayohParham:
                    survivalRate = MossieMethods.GetSurvivingFraction_BayohParham(SliceTemp, SliceHumidity, m_SliceLengthDays, m_UpperTempLimit);
                    break;
                case MossieMethods.SurvivalFunctions.BayohLunde:
                    survivalRate = MossieMethods.GetSurvivingFraction_BayohLunde(SliceTemp, SliceHumidity, m_SliceLengthDays, m_UpperTempLimit);
                    break;
                default:
                    throw new InvalidOperationException("Unrecognised survival function / not implemented");
            }
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


    /// <summary>
    /// Static implementations of various functions for calculating the proportion of mosquitos to survive a given number 
    /// of days. The functions are all extracted from summaries published in :
    /// "How malaria models relate temperature to malaria transmission" (Lunde, Bayoh & Bernt Lindtjørn, 2013)
    /// DOI: 10.1186/1756-3305-6-20.
    /// All functions return a value between 0 and 1 representing a proportion of the existing population to survive 
    /// a timeslice of the given length based on the input temperature (and humidity) values. This takes no account 
    /// of natural lifespan i.e. if the requested timeslice is longer than a lifespan, that won't be accounted for.
    /// Four of the functions depend only on temperature (in degrees celsius) and the other two also require 
    /// relative humidity (in percent).
    /// All take a "death temperature" as a hard upper limit for survival, for comparison with previous results. 
    /// </summary>
    static class MossieMethods
    {
        public enum SurvivalFunctions
        {
            Martens2,
            Martens3,
            BayohMordecai,
            BayohErmert,
            BayohParham,
            BayohLunde
        }
        //TempSuitability_CSharp.MossieMethods.SurvivalFunctions
        /// <summary>
        /// Implements the "Martens equation" for mosquito survival as a function of air temperature, published 
        /// as a PhD thesis at 
        /// https://cris.maastrichtuniversity.nl/portal/files/1425195/guid-1dbc4a8c-a914-41b6-ab86-ea00cf972d9a-ASSET1.0
        ///     (Martens W: Health impacts of climate change and ozone depletion: an eco-epidemiological modelling approach. 
        ///     PhD thesis. 1997, Maastricht, Netherlands: Maastricht University)
        /// Of note is that this equation was suggested empirically based on _three_ data points from the literature 
        /// from 1949, 1949, and 1981. So it's not like it's the bible.
        /// A slightly simplified version of this function (without the e term) was the one used by Gething et al in the 
        /// original non-dynamic temperature suitability mapping, and this function as it is was used by Weiss et al.
        /// </summary>
        /// <param name="sliceTemp"></param>
        /// <param name="sliceLengthDays"></param>
        /// <param name="deathTemp"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double GetSurvivingFraction_Martens2(double sliceTemp, double sliceLengthDays, double deathTemp)
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

        /// <summary>
        /// This implements a re-parameterised version of the "Martens equation" for survival as a function of 
        /// air temperature, as published in the summary paper referenced in the class description. The equation form
        /// is identical to that implemented in GetSurvivingFraction_Martens2 but the constant values are slightly altered.
        /// </summary>
        /// <param name="sliceTemp"></param>
        /// <param name="sliceLengthDays"></param>
        /// <param name="deathTemp"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double GetSurvivingFraction_Martens3(double sliceTemp, double sliceLengthDays, double deathTemp)
        {
            if (sliceTemp > deathTemp)
            {
                return 0;
            }

            var f = Math.Pow(
             (Math.Exp(-1 / (-4.31564 + (2.19646 * sliceTemp) - (0.058276 * (Math.Pow(sliceTemp, 2)))))),
             (sliceLengthDays));
            return f < 0 ? 0 : f >= 1 ? 0 : f;
        }

        /// <summary>
        /// This implementes the "Bayou-Mordecai" function for mosquito survival as a function of temperature, as 
        /// published in the summary paper referenced in the class description. As there is no exponential term this 
        /// model is likely to suffer from inappropriately prolonged survival at high temperatures.
        /// </summary>
        /// <param name="sliceTemp"></param>
        /// <param name="sliceLengthDays"></param>
        /// <param name="deathTemp"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double GetSurvivingFraction_BayohMordecai(double sliceTemp, double sliceLengthDays, double deathTemp)
        {
            // = morP at line 24 of mordecai.R
            if (sliceTemp > deathTemp)
            {
                return 0;
            }
            var daily = -0.000828 * (Math.Pow(sliceTemp, 2)) + 0.0367 * sliceTemp + 0.522;
            var f = Math.Pow(daily, sliceLengthDays);
            return f < 0 ? 0 : f >= 1 ? 0 : f;
        }

        /// <summary>
        /// This implementes the "Bayou-Ermert" function for mosquito survival as a function of temperature, as 
        /// published in the summary paper referenced in the class description. As there is no exponential term this 
        /// model is likely to suffer from inappropriately prolonged survival at high temperatures.
        /// </summary>
        /// <param name="sliceTemp"></param>
        /// <param name="sliceLengthDays"></param>
        /// <param name="deathTemp"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double GetSurvivingFraction_BayohErmert(double sliceTemp, double sliceLengthDays, double deathTemp)
        {
            if (sliceTemp > deathTemp)
            {
                return 0;
            }
            var daily =
                -2.123e-7 * Math.Pow(sliceTemp, 5)
                + 1.961e-5 * Math.Pow(sliceTemp, 4)
                - 6.394e-4 * Math.Pow(sliceTemp, 3)
                + Math.Pow(8.217, -3) * Math.Pow(sliceTemp, 2)
                - 1.865e-2 * sliceTemp
                + 7.238e-1;
            var f = Math.Pow(daily, sliceLengthDays);
            return f < 0 ? 0 : f >= 1 ? 0 : f;
        }

        /// <summary>
        /// Implements equation 16 from 
        /// Modeling the role of environmental variables on the population dynamics of the malaria vector Anopheles gambiae sensu stricto
        /// by Parham et al (DOI: DOI: 10.1186/1475-2875-11-271)
        /// Referred to as "Bayoh-Parham" model in the summary paper referenced in class description.
        /// </summary>
        /// <param name="sliceTemp"></param>
        /// <param name="rhPercent"></param>
        /// <param name="sliceLengthDays"></param>
        /// <param name="deathTemp"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double GetSurvivingFraction_BayohParham(double sliceTemp, double rhPercent, double sliceLengthDays, double deathTemp)
        {
            if (sliceTemp > deathTemp)
            {
                return 0;
            }
            var b0 = 0.00113 * Math.Pow(rhPercent, 2) - 0.158 * rhPercent - 6.61;
            var b1 = -2.32e-4 * Math.Pow(rhPercent, 2) + 0.0515 * rhPercent + 1.06;
            var b2 = 4e-6 * Math.Pow(rhPercent, 2) - 1.09e-3 * rhPercent - 0.0255;
            var daily = Math.Exp(-1 / ((b2 * Math.Pow(sliceTemp, 2)) + (b1 * sliceTemp) + b0));
            var f = Math.Pow(daily, sliceLengthDays);
            return f < 0 ? 0 : f >= 1 ? 0 : f;
        }

        /// <summary>
        /// Implements the Bayoh-Lunde model for mosquito survival as a function of temperature and relative humidity, introduced in the paper 
        /// referenced in the class description. The term for mosquito size is not included in this implementation, corresponding to a 
        /// mosquito size of 3.05mm where the original term calculates to 1.
        /// </summary>
        /// <param name="sliceTemp"></param>
        /// <param name="rhPercent"></param>
        /// <param name="sliceLengthDays"></param>
        /// <param name="deathTemp"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double GetSurvivingFraction_BayohLunde(double sliceTemp, double rhPercent, double sliceLengthDays, double deathTemp)
        {
            if (sliceTemp > deathTemp)
            {
                return 0;
            }
            var fRH = 6.48007 + .69570 * (1 - Math.Exp(-0.06 * rhPercent));
            var cmn = 1 + (sliceTemp + 1) / 21;
            var brack1 = Math.Pow(cmn, (2.0/3));
            var brack2 = Math.Pow(cmn, 2);
            var pwr = 10 + brack1 * (brack2 - (cmn * 2) - fRH);
            var alpha = Math.Exp(pwr); // neglecting the g term for mosquito size

            var zeta = 6;
            float[] zetafacs = new float[] { 1, 1, 2, 6, 24, 120 };
            
            double daily = 0;
            for (int i = 0; i < zeta; i++)
            {
                daily += (Math.Pow(alpha, i) / zetafacs[i]) * Math.Exp(-alpha);
            }
            // given a fraction surviving one day, calculate the proportion surviving a sub-daily slice
            var f = Math.Pow(daily, sliceLengthDays);
            return f < 0 ? 0 : f >= 1 ? 0 : f;
        }
    }
}
