using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TempSuitability_CSharp
{

    /// <summary>
    /// Basic struct to pass TS model parameters relating to model setup for a particular cell
    /// </summary>
    struct GeographicCellLocation
    {
        public double Latitude { get; set; }
        public double Longitude { get; set; }
        public double CellSizeDegrees { get; set; }
    }

    class GeographicCell
    {
        private readonly GeographicCellLocation locationParams;
        public GeographicCell(GeographicCellLocation CellLocation)
        {
            this.locationParams = CellLocation;
        }

        /// <summary>
        /// Calculates the approximate number of daylight hours on the given day of the year at the current location.
        /// Calculated according to the "CBD" model published in:
        /// Forsythe et al. (1995) A model comparison for daylength as a function of latitude and day of year.  
        /// Ecological Modeling 80(1) pp. 87-95
        /// </summary>
        /// <param name="JulianDay"></param>
        /// <returns></returns>
        public double CalcDaylightHrsForsyth(int JulianDay)
        {
            var lat = locationParams.Latitude;

            const double daylengthCoefficient = 0.8333;
            var theta = 0.2163108 + 2 * Math.Atan(0.9671396 * Math.Tan(0.00860 * (JulianDay - 186)));
            var phi = Math.Asin(0.39795 * Math.Cos(theta));
            var hrs = 24 - (24 / Math.PI) * Math.Acos(
                (Math.Sin(daylengthCoefficient * Math.PI / 180) + Math.Sin(lat * Math.PI / 180) * Math.Sin(phi))
                /
                (Math.Cos(lat * Math.PI / 180) * Math.Cos(phi)));
            return hrs;
        }

        public List<double> CalcDailyDaylightHrs()
        {
            var res = new List<double>(365);
            for (int i = 0; i< 365; i++)
            {
                res[i] = CalcDaylightHrsForsyth(i);
            }
            return res;
        }
        public List<Tuple<double, double>> CalcDailySunriseSunset()
        {
            var res = new List<Tuple<double, double>>(365);
            for (int i = 0; i < 365; i++)
            {
                res[i] = GetSunriseSunsetTimes(i);
            }
            return res;
        }
        public List<double> CalcMonthlyAvgDaylightHrs()
        {
            var res = new List<double>(12);
            var n = new List<int>(12);
            var sampleDate = new DateTime(2001, 1, 1);
            var oneDay = new TimeSpan(1, 0, 0, 0);
            for (int i = 0; i<365; i++)
            {
                var mth = sampleDate.Month;
                res[mth] += CalcDaylightHrsForsyth(i);
                n[mth] += 1;
            }
            for (int i = 0; i< 12; i++)
            {
                res[i] /= n[i]; 
            }
            return res;
        }

        private enum SunriseOrSunset
        {
            Sunrise,
            Sunset
        }

        /// <summary>
        /// Implementation of the algorithm at http://williams.best.vwh.net/sunrise_sunset_algorithm.htm
        /// </summary>
        /// <param name="JulianDay"></param>
        /// <param name="which"></param>
        /// <returns>double representing the hour of sunrise or sunset, according to which was requested.
        /// A value of -9999 is returned if the sun never rises or sets on this day.</returns>
        private double GetSunriseOrSunsetTime(int JulianDay, SunriseOrSunset which, double? lonDegrees, double? latDegrees)
        {
            double lat, lon;
            lat = latDegrees.HasValue ? latDegrees.Value : locationParams.Latitude;
            lon = lonDegrees.HasValue ? lonDegrees.Value : locationParams.Longitude;

            double longitude_hr = lon / 15;
            const double Zenith = 90.8333;
            double degConv = Math.PI / 180;

            double timeInit;
            switch (which)
            {
                case SunriseOrSunset.Sunrise:
                    timeInit = JulianDay + ((6 - longitude_hr) / 24);
                    break;
                case SunriseOrSunset.Sunset:
                    timeInit = JulianDay + ((18 - longitude_hr) / 24);
                    break;
                default:
                    throw new ArgumentException();
            }
            double meanAnom = (0.9856 * timeInit) - 3.289;

            double sunLong = meanAnom +
                (1.916 * Math.Sin(meanAnom * degConv)) +
                (0.020 * Math.Sin(2 * meanAnom * degConv)) +
                282.634;
            sunLong = (sunLong + 360) % 360;

            double rightAscension = (180 / Math.PI) *
                Math.Atan(0.91764 * Math.Tan(sunLong * degConv));
            rightAscension = (rightAscension + 360) % 360;
            double l_quadrant = (Math.Floor(sunLong / 90)) * 90;
            double ra_quadrant = (Math.Floor(rightAscension / 90)) * 90;
            rightAscension += (l_quadrant - ra_quadrant);
            rightAscension /= 15;

            double sinDecl = 0.39782 * Math.Sin(sunLong * degConv);

            double cosDecl = Math.Cos((Math.Asin(sinDecl) * 180 / Math.PI) * degConv);

            double sunHourAngle = (Math.Cos(Zenith * degConv) -
                (sinDecl * Math.Sin(lat * degConv))) /
                (cosDecl * Math.Cos(lat * degConv));

            if (sunHourAngle > 1 || sunHourAngle < -1)
            {
                return -9999;
            }
            double time;
            switch (which)
            {
                case SunriseOrSunset.Sunrise:
                    time = 360 - Math.Acos(sunHourAngle) * 180 / Math.PI;
                    break;
                case SunriseOrSunset.Sunset:
                    time = Math.Acos(sunHourAngle) * 180 / Math.PI;
                    break;
                default:
                    throw new ArgumentException();
            }
            time /= 15;
            time = time + rightAscension - (0.06571 * timeInit) - 6.622;
            // this would put the time into UTC
            // time = time - longitude_hr;
            // but we don't want to do that as we will iterate throug hours of the day 0-24 i.e. local time 
            if (time > 24)
            {
                time -= 24;
            }
            else if (time < 0)
            {
                time += 24;
            }
            return time;
        }

        /// <summary>
        /// Calculates sunrise and sunset times at the current location for a given day of the year
        /// </summary>
        /// <param name="JulianDay"></param>
        /// <returns>2-Tuple containing the sunrise and sunset times for this Julian day.
        /// A value of -9999 is returned if the sun never rises or sets on this day.</returns>
        public Tuple<double, double> GetSunriseSunsetTimes(int JulianDay)
        {
            return new Tuple<double, double>(
                GetSunriseOrSunsetTime(JulianDay, SunriseOrSunset.Sunrise, null, null),
                GetSunriseOrSunsetTime(JulianDay, SunriseOrSunset.Sunset, null, null));
        }


    }
}
