using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace TempSuitability_CSharp
{
    /// <summary>
    /// A class to parse the date encoded in filenames in the MODIS 8-daily source data format, which begin with 
    /// AYYYYJJJ..... where the letter A is constant and is followed by a four-digit year and a three-digit julian day
    /// </summary>
    class FilenameDateParser_MODIS8DayRaw : IFilenameDateParser
    {
        public DateTime? TryParseFilenameDate(string Filename)
        {
            try
            {
                string pat = @"^A(\d+)";
                Regex r = new Regex(pat);
                string basename = System.IO.Path.GetFileName(Filename);
                string datedigits = r.Match(basename).Groups[1].ToString();
                int yearnum = Convert.ToInt32(datedigits.Substring(0, 4));
                int juliannum = Convert.ToInt32(datedigits.Substring(4));
                DateTime newYrDay = new DateTime(yearnum, 1, 1);
                DateTime newDate = newYrDay.AddDays(juliannum - 1);
                return newDate;
            }
            catch
            {
                return null;
            }
        }
    }

    /// <summary>
    /// A class to parse the date encoded in filenames complying to the mastergrid naming convention 
    /// </summary>
    class FilenameDateParser_Mastergrid : IFilenameDateParser
    {
        private int _forceDayOfMonth;

        public FilenameDateParser_Mastergrid(int? overrideOrImplicitDayOfMonth = 15)
        {
            if (!overrideOrImplicitDayOfMonth.HasValue)
            {
                _forceDayOfMonth = -1;
            }
            else if (overrideOrImplicitDayOfMonth.Value < 1 || overrideOrImplicitDayOfMonth.Value > 31)
            {
                throw new ArgumentOutOfRangeException("months have between 1 and 31 days fool");
            }
            else
            {
                _forceDayOfMonth = overrideOrImplicitDayOfMonth.Value;
            }
        }

        /// <summary>
        /// Attempts to parse and return the date of a file names according to the mastergrid naming convention. 
        /// Synoptic files will return null as these do not have a date associated. Monthly files (mathcing the pattern YYYY.MM) 
        /// will return a date with the day set to whatever was configured in the constructor: 15 by default. 
        /// Otherwise the date can be specified as YYYY.JJJ where JJJ is the julian day-of-year, or YYYY.MMDD.
        /// </summary>
        /// <param name="Filename"></param>
        /// <returns></returns>
        public DateTime? TryParseFilenameDate(string Filename)
        {
            string basename = System.IO.Path.GetFileName(Filename);
            var parts = basename.Split('.');
            if (parts.Length < 6)
            {
                throw new ArgumentException("Does not appear to be a valid mastergrid filename, which should have 6 dot-delimited tokens: "
                    + Filename);
            }
            int yearnum;
            if (parts[1] == "Synoptic")
            {
                return null;
            }
            else
            {
                yearnum = Convert.ToInt32(parts[1]);
            }

            var monthOrDayOrDate = parts[2];
            int monthNum;
            int dayNum;
            DateTime outDate;
            if (monthOrDayOrDate == "Overall")
            {
                throw new ArgumentException("3rd token can only be 'Overall' if 2nd token is 'Synoptic' (and then no date can be parsed)");
            }
            if (monthOrDayOrDate.Length == 2)
            {
                // it should be a month, set day to 15 unless overridden on instantiation; 
                // invalid values e.g. month not 1-12 will throw error later
                monthNum = Convert.ToInt32(monthOrDayOrDate);
                dayNum = _forceDayOfMonth == -1 ? 15 : _forceDayOfMonth;
            }
            else if (monthOrDayOrDate.Length == 4)
            {
                // it should be a month and day i.e. MMDD
                monthNum = Convert.ToInt32(monthOrDayOrDate.Substring(0, 2));
                dayNum = Convert.ToInt32(monthOrDayOrDate.Substring(2, 2));
            }
            else if (monthOrDayOrDate.Length == 3)
            {
                var juliannum = Convert.ToInt32(monthOrDayOrDate);
                var newYrDay = new DateTime(yearnum, 1, 1);
                outDate = newYrDay.AddDays(juliannum - 1);
                return outDate;
            }
            else
            {
                throw new ArgumentException("3rd token must have 2, 3, or 4 digits or be 'Overall'");
            }
            outDate = new DateTime(yearnum, monthNum, dayNum);
            return outDate;
        }
    }

    class FilenameDateParser_AutoDetect : IFilenameDateParser
    {
        public DateTime? TryParseFilenameDate(string Filename)
        {
            DateTime? parsedDate = null;
            try
            {
                IFilenameDateParser mgParser = new FilenameDateParser_MODIS8DayRaw();
                parsedDate = mgParser.TryParseFilenameDate(Filename);
            }
            catch (ArgumentException e)
            {

            }
            if (parsedDate != null)
            {
                return parsedDate;
            }
            IFilenameDateParser rawParser = new FilenameDateParser_Mastergrid();
            return rawParser.TryParseFilenameDate(Filename);
        }
    }
}
    
