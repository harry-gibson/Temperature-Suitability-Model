using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace TempSuitability_CSharp
{
    class TempSuitabilityModelTest
    {
        static double[,] RunForStripe(string FilePattern)
        {
            string fnPattern = "C:\\Users\\zool1301\\Documents\\Code_General\\temp-suitability\\Original\\c_code\\test5\\data\\{0}.{1}.dat";
            string latFile = String.Format(fnPattern, "lat", FilePattern);
            string lonFile = String.Format(fnPattern, "long", FilePattern);
            string minFile = String.Format(fnPattern, "tmin", FilePattern);
            string maxFile = String.Format(fnPattern, "tmax", FilePattern);
            double[] tLats = File.ReadAllLines(latFile).
                Select(l => double.Parse(l)).ToArray();
            double[] tLongs = File.ReadAllLines(lonFile).
                Select(l => double.Parse(l)).ToArray();
            double[] tMins = File.ReadAllLines(minFile).
                Select(l => double.Parse(l)).ToArray();
            double[] tMaxs = File.ReadAllLines(maxFile).
                Select(l => double.Parse(l)).ToArray();
            int numCells = tLats.Length;
            int nDays = tMins.Length / tLats.Length;
            PopulationParams popParams = new PopulationParams();
            popParams.DegreeDayThreshold = 111;
            popParams.MinTempThreshold = 11;
            popParams.MosquitoLifespanDays = 31;
            popParams.SlicesPerDay = 12;
            double[,] tsOut = new double[numCells, nDays];
            //int c = 0;
            Parallel.For(0, numCells, c =>
            {
                double lat = tLats[c];
                double lon = tLongs[c];
                int startIdx = c * nDays;
                int endIdx = startIdx + nDays;
                double[] tMinsForCell = new double[nDays];
                double[] tMaxsForCell = new double[nDays];
                Array.Copy(tMins, startIdx, tMinsForCell, 0, nDays);
                Array.Copy(tMaxs, startIdx, tMaxsForCell, 0, nDays);

                SpatioTemporalParams geogParams = new SpatioTemporalParams();
                geogParams.Latitude = tLats[c];
                geogParams.Longitude = tLongs[c];
                geogParams.ModelRuntimeDays = nDays;
                geogParams.StartJulianDay = 0;

                TSCellModel tsModel = new TSCellModel(popParams, geogParams);
                double[] tsCell = tsModel.RunModel(tMinsForCell, tMaxsForCell);
                for (int i = 0; i < nDays; i++)
                {
                    tsOut[c, i] = tsCell[i];
                }
            });
            return tsOut;
        }
        static void Main(string[] args)
        {
            
            for (int i = 0; i <= 45; i += 5)
            {
                var res = RunForStripe(i.ToString());
               /* string outPath = String.Format("C:\\temp\\{0}.TestOutput.dat", i);
                for (int i = 0; i < res.GetLength(0); i++)
                {

                }*/
            }
            System.Console.WriteLine("Done");
        }
    }
}
