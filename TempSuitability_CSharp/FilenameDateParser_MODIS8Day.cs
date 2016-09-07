using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace TempSuitability_CSharp
{
    class FilenameDateParser_MODIS8Day : IFilenameDateParser
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
}
