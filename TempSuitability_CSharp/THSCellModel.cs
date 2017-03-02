using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Interpolation;

namespace TempSuitability_CSharp
{
    class THSCellModel:TSCellModel
    {
        protected IInterpolation _HumiditySpline;
        protected bool _CanRunHum = false;
        protected new ITempHumidityIterablePopulation m_Pop;
        
        protected new bool _CanRun
        {
            get
            {
                if (m_Pop.GetRequiresHumidity())
                {
                    return _CanRunTemps && _CanRunHum;
                }
                return _CanRunTemps;
            }
        }

        public THSCellModel(PopulationParams modelConfigParams, GeographicCellLocation locationParams) :
            base(modelConfigParams, locationParams, PopulationTypes.Pointers)
        {
            m_Pop = new ConfigurableDecayPopulation(modelConfigParams);
        }

        /// <summary>
        /// Loads the humidity data when it is provided as temporally-dynamic data. 
        /// The temperature series must be loaded first! The humidity series does not have to span the same 
        /// period as the temperature series, but it is recommended to be at least as long as the temperature 
        /// series as otherwise the spline extrapolation of the humidity data into the range of the temperature 
        /// data may be dodgy!
        /// </summary>
        public bool SetHumidityData (float[] HumidityPercentValues, DateTime[] HumidityDatePoints, float HumidityNoDataValue)
        {
            if (HumidityDatePoints.Length != HumidityPercentValues.Length)
            {
                return false;
            }
            // we will require the temperatures to be loaded first, the model run period will be determined by the temperatures
            if (m_InputTemperatureDates.Length == 0)
            {
                return false;
            }
            int nPoints = HumidityPercentValues.Length;
            // check they are sorted
            for (int i = 1; i < nPoints; i++)
            {
                if (HumidityDatePoints[i - 1] > HumidityDatePoints[i])
                {
                    return false;
                }
            }
            var startDay = m_InputTemperatureDates.First();
            List<double> validHumidities, validHumidityDates;
            validHumidities = new List<double>(nPoints);
            validHumidityDates = new List<double>(nPoints);
            // strip out any nodata ones as these can't be handled by the spline object, and convert 
            // the dates to seconds since model runtime start point
            for (int i = 0; i < nPoints; i++)
            {
                if (HumidityPercentValues[i] != HumidityNoDataValue)
                {
                    var dSeconds = (HumidityDatePoints[i] - startDay).TotalSeconds;
                    validHumidities.Add(HumidityPercentValues[i]);
                    validHumidityDates.Add(dSeconds);
                }
            }
            if (validHumidities.Count / HumidityPercentValues.Length < modelParams.ValidDataProportion)
            {
                return false;
            }
            _HumiditySpline = CubicSpline.InterpolateNaturalSorted(validHumidityDates.ToArray(), validHumidities.ToArray());
            _CanRunHum = true;
            return true;
        }

        /// <summary>
        /// Loads the humidity data when it is provided as synoptic data, i.e. one point for each calendar month,
        /// all years assumed to be identical. This will repeat the supplied synoptic series over the full timespan 
        /// of the model run as defined by the temperature series. The temperature series must be loaded first!
        /// </summary>
        /// <param name="HumidityPercentValues"></param>
        /// <param name="HumidityNoDataValue"></param>
        /// <returns></returns>
        public bool SetSynopticHumidityData(float[] HumidityPercentValues, float HumidityNoDataValue)
        {
            if(HumidityPercentValues.Length != 12)
            {
                return false;
            }
            if (m_InputTemperatureDates.Length == 0)
            {
                return false;
            }
            // get a datetime object for the 15th of every month spanned by the input non-synoptic temperature data
            var tempDataMonths = m_InputTemperatureDates             // all the input dates (e.g. 8 daily)
                   .Select(d => new DateTime(d.Year, d.Month, 15))   // change each of them to the first of that month
                   .Distinct()                                      // get unique
                   .OrderBy(d => d)                                 // sort
                   .ToArray();                                      // return as array
            var nMonths = tempDataMonths.Length;
            var humRepeatedSeries = new float[nMonths];
            for (int i = 0; i < nMonths; i++)
            {
                humRepeatedSeries[i] = HumidityPercentValues[tempDataMonths[i].Month-1];
            }
            return SetHumidityData(humRepeatedSeries, tempDataMonths, HumidityNoDataValue);
        }

        public new float[] RunModel()
        {
            if (!_CanRun)
            {
                throw new InvalidOperationException("Data has not yet been loaded");
            }
            if (!_PotentiallySuitableTemperature)
            {
                // IF the temperature never enters suitable range, return a new array (which will have the 
                // default value of zero, i.e. no temperature suitability) andof the same length as the 
                // output would have been i.e. the number of months in the input series excluding the 
                // first lifespan-worth.
                // The output dates will be calculated on demand at this point.
                int nMonths = OutputDates.Count();
                return new float[nMonths];
            }
            var startDay = m_InputTemperatureDates.First();
            var endDay = m_InputTemperatureDates.Last();
            var mainRunStartDate = startDay + modelParams.MosquitoLifespanDays;

            TimeSpan oneDay = new TimeSpan(1, 0, 0, 0, 0);
            double sliceLengthDays = modelParams.SliceLength.TotalDays;
            int slicesPerDay = (int)(1.0 / sliceLengthDays);
            int sliceLengthHours = (int)modelParams.SliceLength.TotalHours;
            int slicesPerLifespan = modelParams.LifespanSlices; // this rounds away any fractional days of the lifespan

            double[] tSpinUpTempData = new double[slicesPerLifespan];
            double[] tSpinUpHumData = new double[0];
            bool requiresHumidity = m_Pop.GetRequiresHumidity();
            if (requiresHumidity)
            {
                tSpinUpHumData = new double[slicesPerLifespan];
            }

            double maxAirTemp, minAirTemp;
            double humPct = 0;
            // get a local ref to the spline objects (just fractionally more efficient)
            var maxSpline = _MaxTempSpline;
            var minSpline = _MinTempSpline;
            var humSpline = _HumiditySpline;
      
            // We will run through the overall data timespan from start to end, on a timestep of one calendar day. 
            // At each calendar day we will interpolate a max and min temperature value from the 8-daily converted 
            // input series, then we will run through the sub-daily (e.g. 2-hour) slices that comprise a day and 
            // interpolate temperature and calculate suitability at each such slice.

            // First "spin up" the model by running through the first mosquito-lifespan-worth of data so we 
            // have a stabilised population. Generate an array with a temperature for each 2-hour timeslice of a
            // a mosquito's life for a one-shot initialise call
            int spinupSliceNum = 0;
            for (DateTime currentDay = startDay; currentDay < mainRunStartDate; currentDay += oneDay)
            {
                TimeSpan timeSinceStart = (currentDay - startDay);
                double secondsSinceStart = timeSinceStart.TotalSeconds;
                int runDay = timeSinceStart.Days;
                // interpolate the temp for this calendar day, from the 8-or-whatever-daily input series
                maxAirTemp = maxSpline.Interpolate(secondsSinceStart);
                minAirTemp = minSpline.Interpolate(secondsSinceStart);
                if (requiresHumidity)
                {
                    humPct = humSpline.Interpolate(secondsSinceStart);
                }
                int julianDay = currentDay.DayOfYear;
                Tuple<double, double> sunriseSunset = modelLocation.GetSunriseSunsetTimes(julianDay);
                for (TimeSpan daySlice = new TimeSpan(0); daySlice < oneDay; daySlice += modelParams.SliceLength)
                {
                    double currentTemp = InterpolateHourlyTemperature(daySlice.Hours, sunriseSunset.Item1, sunriseSunset.Item2, minAirTemp, maxAirTemp);
                    tSpinUpTempData[spinupSliceNum] = currentTemp;
                    if (requiresHumidity)
                    {
                        tSpinUpHumData[spinupSliceNum] = humPct;
                    }
                    spinupSliceNum += 1;
                }
            }
            // spin up the population model
            if (requiresHumidity)
            {
                // calls initialise as defined in the temperature+humidity-dependent population modeller class
                m_Pop.Initialise(tSpinUpTempData, tSpinUpHumData);
            }
            else
            {
                // calls initialise as defined in the base (temperature only) pop modeller class
                m_Pop.Initialise(tSpinUpTempData);
            }

            // Now, run the model for the (entire) remainder of the time period, tracking the results on a monthly basis
            SortedList<DateTime, double> results = new SortedList<DateTime, double>();
            Dictionary<DateTime, int> daysInMonth = new Dictionary<DateTime, int>();
            // Iterate through the time period one calendar day at a time
            for (DateTime currentDay = mainRunStartDate; currentDay <= endDay; currentDay += oneDay)
            {
                double tsOfDay = 0;
                TimeSpan timeSinceStart = (currentDay - startDay);
                double secondsSinceStart = timeSinceStart.TotalSeconds;
                int runDay = timeSinceStart.Days;
                DateTime currentMonth = new DateTime(currentDay.Year, currentDay.Month, 1);

                // interpolate the temp for this calendar day, from the 8-or-whatever-daily input series
                maxAirTemp = maxSpline.Interpolate(secondsSinceStart);
                minAirTemp = minSpline.Interpolate(secondsSinceStart);
                if (requiresHumidity)
                {
                    humPct = humSpline.Interpolate(secondsSinceStart);
                }
                int julianDay = currentDay.DayOfYear;
                Tuple<double, double> sunriseSunset = modelLocation.GetSunriseSunsetTimes(julianDay);
                // go through each e.g. 12 slices of this day
                for (TimeSpan daySlice = new TimeSpan(0); daySlice < oneDay; daySlice += modelParams.SliceLength)
                {
                    double currentTemp = InterpolateHourlyTemperature(daySlice.Hours, sunriseSunset.Item1, sunriseSunset.Item2, minAirTemp, maxAirTemp);
                    if (requiresHumidity)
                    {
                        tsOfDay += m_Pop.Iterate(currentTemp, humPct);
                    }
                    else
                    {
                        tsOfDay += m_Pop.Iterate(currentTemp);
                    }
                }
                // for each calendar month record the total (sum) temp suitability of all the 2-hour slices in that month
                if (results.ContainsKey(currentMonth))
                {
                    results[currentMonth] += tsOfDay;
                    // also count the number of calendar days in each month that we've calculated a value for (not necessarily a whole 
                    // month, at the start and end, and february varies too)
                    daysInMonth[currentMonth] += 1;
                }
                else
                {
                    results[currentMonth] = tsOfDay;
                    daysInMonth[currentMonth] = 1;
                }
            }
            var months = results.Keys.ToList();
            // Take the average slice value for the month, not the total
            foreach (var mth in months)
            {
                var divisor = (daysInMonth[mth] * slicesPerDay * modelParams.MaxTempSuitability);
                results[mth] /= divisor;
            }
            // store the months for which we've generated outputs so the OutputDates property doesn't have to calculate 
            // them if it's called
            _iteratedOutputDates = results.Keys.ToArray();
            // output as float (we ran using doubles, but no need to maintain that as input / output tiffs will be float32)
            return results.Values.Select(d => (float)d).ToArray();
        }
    }
}
