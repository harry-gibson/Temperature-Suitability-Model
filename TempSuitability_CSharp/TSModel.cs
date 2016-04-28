using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace TempSuitability_CSharp
{
    /// <summary>
    /// Basic struct to pass TS model parameters relating to population dynamics and response
    /// </summary>
    struct PopulationParams
    {
        /// <summary>
        /// The maximum lifespan of a mosquito (they may die earlier according to temperature)
        /// </summary>
        public int MosquitoLifespanDays { get;  set; }
        /// <summary>
        /// The temporal resolution with which to run the model, e.g. 12 = 2 hour slices
        /// </summary>
        public int SlicesPerDay { get;  set; }
        /// <summary>
        /// The threshold temperature above which degree days start accumulating
        /// </summary>
        public double MinTempThreshold { get;  set; }
        /// <summary>
        /// The number of accumulated degree days above which a mosquito will be 
        /// deemed to become infectious. One degree day is a temperature of one degree 
        /// above the min temp threshold maintained for one full day.
        /// </summary>
        public double DegreeDayThreshold { get;  set; }
    }

    /// <summary>
    /// Basic struct to pass TS model parameters relating to model setup for a particular cell
    /// </summary>
    struct SpatioTemporalParams
    {
        public double Latitude { get;  set; }
        public double Longitude { get;  set; }
        public int ModelRuntimeDays { get;  set; }
        public int StartJulianDay { get;  set; }
    }

    
    /// <summary>
    /// Class to run the temperature suitability model for a given cell (location), for the entire time period. 
    /// </summary>
    class TSCellModel
    {
        private readonly PopulationParams modelParams;
        private readonly SpatioTemporalParams timeParams;
        private Population m_Pop;

        private double m_PreviousSunsetTemp;

        public TSCellModel(PopulationParams modelParams, SpatioTemporalParams timeParams)
        {
            this.modelParams = modelParams;
            this.timeParams = timeParams;
            double sliceLengthDays = 1.0 / modelParams.SlicesPerDay;
            int lifespanSlices = modelParams.SlicesPerDay * modelParams.MosquitoLifespanDays;
            m_Pop = new Population(lifespanSlices, sliceLengthDays, modelParams.MinTempThreshold, modelParams.DegreeDayThreshold);
        }

        /// <summary>
        /// Currently takes arrays with one value for each calendar day, for testing only.
        /// Will be modified to take one value for each source image (8-daily) and do the splining here.
        /// </summary>
        /// <param name="MinTemps">Array of LST Night temperatures, one per each calendar day</param>
        /// <param name="MaxTemps">Array of LST Day temperatires, one per each calendar day</param>
        /// <returns></returns>
        public double[] RunModel(double[] LST_Day_Temps, double[] LST_Night_Temps)
        {
            if (LST_Day_Temps.Length != LST_Night_Temps.Length || LST_Day_Temps.Length != timeParams.ModelRuntimeDays)
            {
                throw new ArgumentException("MinTemps and MaxTemps array must currently both be of same size as model runtime");

            }
            double sliceLengthDays = 1 / modelParams.SlicesPerDay;
            int sliceLengthHours = 24 / modelParams.SlicesPerDay;
            int lifespanSlices = modelParams.SlicesPerDay * modelParams.MosquitoLifespanDays;

            double[] tsOfDays = new double[timeParams.ModelRuntimeDays];
            double[] tSpinUpData = new double[lifespanSlices];

            double lstDayTemp = LST_Day_Temps[0];
            double lstNightTemp = LST_Night_Temps[0];
            double maxAirTemp, minAirTemp;
            ConvertTemperature(lstDayTemp, lstNightTemp, timeParams.StartJulianDay, out m_PreviousSunsetTemp, out minAirTemp);
            
            // "Spin up" the population model at this cell i.e. run it for an entire mosquito lifespan
            for (int runDay = 0; runDay<modelParams.MosquitoLifespanDays; runDay++)
            {
                int julianDay = (timeParams.StartJulianDay + runDay + 365) % 365;
                Tuple<double, double> sunriseSunset = GetSunriseSunsetTimes(julianDay);
                ConvertTemperature(LST_Day_Temps[runDay], LST_Night_Temps[runDay], julianDay, out maxAirTemp, out minAirTemp);
                for (int daySlice  = 0; daySlice < modelParams.SlicesPerDay; daySlice ++ )
                {
                    int hourOfDay = daySlice * sliceLengthHours;
                    double currentTemp = InterpolateHourlyTemperature(hourOfDay, sunriseSunset.Item1, sunriseSunset.Item2, minAirTemp, maxAirTemp);
                    tSpinUpData[runDay * modelParams.SlicesPerDay + daySlice] = currentTemp;
                }
            }
            m_Pop.Initialise(tSpinUpData);

            // Now, run the model for the remainder of the time period
            for (int runDay = modelParams.MosquitoLifespanDays; runDay < timeParams.ModelRuntimeDays; runDay++)
            {
                double tsOfDay = 0;

                int julianDay = (timeParams.StartJulianDay + runDay + 365) % 365;
                Tuple<double, double> sunriseSunset = GetSunriseSunsetTimes(julianDay);
                ConvertTemperature(LST_Day_Temps[runDay], LST_Night_Temps[runDay], julianDay, out maxAirTemp, out minAirTemp);

                // TODO handle case when sun did not rise or set
                for (int daySlice = 0; daySlice < modelParams.SlicesPerDay; daySlice ++)
                {
                    int hourOfDay = daySlice * sliceLengthHours;
                    double currentTemp = InterpolateHourlyTemperature(hourOfDay, 
                        sunriseSunset.Item1, sunriseSunset.Item2, 
                        minAirTemp, maxAirTemp);
                    tsOfDay += m_Pop.Iterate(currentTemp);
                }
                tsOfDays[runDay] = tsOfDay;
            }
            return tsOfDays;
        }
  
        /// <summary>
        /// Implementation of the algorithm at http://williams.best.vwh.net/sunrise_sunset_algorithm.htm
        /// </summary>
        /// <param name="JulianDay"></param>
        /// <returns>2-Tuple containing the sunrise and sunset times for this Julian day.
        /// A value of -9999 is returned if the sun never rises or sets on this day.</returns>
        private Tuple<double, double> GetSunriseSunsetTimes(int JulianDay)
        {
            double longitude_hr = timeParams.Longitude / 15;
            const double Zenith = 90.8333;
            double sunriseTime = -9999, sunsetTime = -9999;

            double sunriseTimeInit = JulianDay + ((6 - longitude_hr) / 24);
            double sunsetTimeInit = JulianDay + ((18 - longitude_hr) / 24);

            double meanAnom = (0.9856 * sunriseTimeInit) - 3.289;
            double sunLong = meanAnom +
                (1.916 * Math.Sin(meanAnom * Math.PI / 180)) +
                (0.020 * Math.Sin(2 * meanAnom * Math.PI / 180)) +
                282.634;
            sunLong = (sunLong + 360) % 360;

            double rightAscension = (180 / Math.PI) * 
                Math.Atan(0.91764 * Math.Tan(sunLong * Math.PI / 180));
            rightAscension = (rightAscension + 360) % 360;

            double l_quadrant = (Math.Floor(sunLong / 90)) * 90;
            double ra_quadrant = (Math.Floor(rightAscension / 90)) * 90;
            rightAscension += (l_quadrant - ra_quadrant);

            rightAscension /= 15;

            double sinDecl = 0.39782 * Math.Sin(sunLong * Math.PI / 180);
            double cosDecl = Math.Cos((Math.Asin(sinDecl) * 180 / Math.PI) * Math.PI / 180);

            double sunHourAngle = (Math.Cos(Zenith * Math.PI / 180) -
                (sinDecl * Math.Sin(timeParams.Latitude * Math.PI / 180))) / 
                (cosDecl * Math.Cos(timeParams.Latitude * Math.PI / 180));
                
            if (sunHourAngle <= 1)
            {
                sunriseTime = 360 - Math.Acos(sunHourAngle) * 180 / Math.PI;
                sunriseTime /= 15;
                sunriseTime = sunriseTime + rightAscension - (0.06571 * sunriseTimeInit) - 6.622;
                sunriseTime = sunriseTime - longitude_hr;
                sunriseTime = (sunriseTime + 24) % 24;
            }
            if (sunHourAngle >= -1)
            {
                sunsetTime = Math.Acos(sunHourAngle) * 180 / Math.PI;
                sunsetTime /= 15;
                sunsetTime = sunsetTime + rightAscension - (0.06571 * sunsetTimeInit) - 6.622;
                sunsetTime = sunsetTime - longitude_hr;
                sunsetTime = (sunsetTime + 24) % 24;
            }
            return new Tuple<double, double>(sunriseTime, sunsetTime);
            //if (sunriseTime != -9999 && sunsetTime != -9999)
            //{
            //    daylengths[i] = sunsetTime - sunriseTime;
            //}
            //else
            //{
            //    // todo, check what we want to say if it never rises or sets
            //    daylengths[i] = 0;
            //}
            //return daylengths;
        }

        /// <summary>
        /// Estimates max and min daily air temperature from the daytime and nighttime land surface temperatures, according 
        /// to the models published in Weiss et al. 2011
        /// </summary>
        /// <param name="LST_Day_Temp"></param>
        /// <param name="LST_Night_Temp"></param>
        /// <param name="JulianDay"></param>
        /// <param name="MaxTemp"></param>
        /// <param name="MinTemp"></param>
        private void ConvertTemperature(double LST_Day_Temp, double LST_Night_Temp, int JulianDay, out double MaxTemp, out double MinTemp)
        {
            double LST_Diff = LST_Day_Temp - LST_Night_Temp;
            MinTemp = 0.209087 * 0.970841 * LST_Night_Temp;
            var tSrSs = GetSunriseSunsetTimes(JulianDay);
            double dayHrs = tSrSs.Item2 - tSrSs.Item1;
            MaxTemp = 18.148887 + LST_Day_Temp* 0.949445 + (LST_Diff * -0.541052) +
                (dayHrs * -0.865620);
        }

        /// <summary>
        /// Estimates the temperature at a given hour of the day based on the daily max / min temperatures and the time 
        /// of sunrise and sunset, according to the model of Garske et al. 2013
        /// </summary>
        /// <param name="hourOfDay"></param>
        /// <param name="sunriseTime"></param>
        /// <param name="sunsetTime"></param>
        /// <param name="dayMinTemp"></param>
        /// <param name="dayMaxTemp"></param>
        /// <returns></returns>
        private double InterpolateHourlyTemperature(double hourOfDay, double sunriseTime, double sunsetTime, double dayMinTemp, double dayMaxTemp)
        {
            double daylightHrs = sunsetTime - sunriseTime;
            double hrTemp;
            if (hourOfDay > sunriseTime && hourOfDay < sunsetTime)
            {
                hrTemp = dayMinTemp + (dayMaxTemp - dayMinTemp) *
                    Math.Sin(Math.PI * (hourOfDay - sunriseTime) / (daylightHrs + 3.72));
                m_PreviousSunsetTemp = hrTemp;
            }
            else if (hourOfDay > sunsetTime)
            {
                hrTemp = dayMinTemp + (m_PreviousSunsetTemp - dayMinTemp) *
                    Math.Exp(-2.2 * ((hourOfDay - sunsetTime) / (24 - daylightHrs)));
            }
            else // (hourOfDay < sunriseSunset.Item2)
            {
                hrTemp = dayMinTemp + (m_PreviousSunsetTemp - dayMinTemp) *
                    Math.Exp(-2.2 * (((hourOfDay + 24) - sunsetTime) / (24 - daylightHrs)));
            }
            return hrTemp;
        }

        private double InterpolateDailyTemp(int JulianDay)
        {
            // PLACEHOLDER
            // Create a field of type CubicSpline and initialise it with input 8-daily temperatures.
            // Use that to interpolate temp at a given day

            // Or see code at https://gist.github.com/dreikanter/3526685
            return 0;
        }
    }
}
