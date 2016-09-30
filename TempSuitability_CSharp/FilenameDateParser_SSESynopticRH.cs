using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Text.RegularExpressions;
namespace TempSuitability_CSharp
{
    class FilenameDateParser_SSESynopticRH : IFilenameDateParser
    {
        private static Dictionary<string, int> MonthNames = new Dictionary<string, int>
        {
            { "Jan",1 },
            { "Feb", 2 },
            { "Mar",3 },
            { "Apr", 4 },
            { "May",5 },
            { "Jun", 6 },
            { "Jul",7 },
            { "Aug", 8 },
            { "Sep",9 },
            { "Oct", 10 },
            { "Nov",11 },
            { "Dec", 12 },
        };

        public DateTime? TryParseFilenameDate(string Filename)
        {
            try
            {

                string pat = @"_(\w+)$";
                Regex r = new Regex(pat);
                string basename = System.IO.Path.GetFileNameWithoutExtension(Filename);
                string monthtxt = r.Match(basename).Groups[1].ToString();
                int monthnum = MonthNames[monthtxt];
                int yearnum = 2000;
                int daynum = 15;
                DateTime newDate = new DateTime(yearnum, monthnum, daynum);
                return newDate;
            }
            catch
            {
                return null;
            }
        }
    }
}
