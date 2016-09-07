using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;

namespace TempSuitability_CSharp
{
    class TempSuitabilityModelTest
    {
        static double[][] RunForStripe(string FilePattern, PopulationTypes ModelType)
        {
            string fnPattern = "C:\\Users\\zool1301.NDPH\\Documents\\Code_General\\temp-suitability\\Original\\c_code\\test5\\data\\{0}.{1}.dat";
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
            popParams.MosquitoLifespanDays = new TimeSpan(31, 0, 0, 0);
            popParams.SliceLength = new TimeSpan(2, 0, 0);
            double[][] tsOut = new double[numCells][];
            //int c = 0;
            Parallel.For(0, numCells, c =>
            {
                tsOut[c] = new double[nDays];
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
                
                TSCellModel tsModel = new TSCellModel(popParams, geogParams, ModelType);
                //double[] tsCell = tsModel.RunModel(tMinsForCell, tMaxsForCell);
                //for (int i = 0; i < nDays; i++)
                //{
                //    tsOut[c][i] = tsCell[i];
                //}
            });
            return tsOut;
        }

        static double[][] RunForArray(
            float[] minTemps, 
            float[] maxTemps,
            float[] latitudes, 
            float[] longitudes, 
            float[] timestamps,
            PopulationTypes ModelType)
        {
            int numCells = latitudes.Length * longitudes.Length;
            int numRows = latitudes.Length;
            int numCols = longitudes.Length;

            Debug.Assert(minTemps.Length == maxTemps.Length);
            Debug.Assert(minTemps.Length % numCells == 0);
            Debug.Assert(minTemps.Length / numCells == timestamps.Length);

            PopulationParams popParams = new PopulationParams();
            popParams.DegreeDayThreshold = 111;
            popParams.MinTempThreshold = 11;
            popParams.MosquitoLifespanDays = new TimeSpan(31, 0, 0, 0); ;
            popParams.SliceLength = new TimeSpan(2, 0, 0);
            
            double[][] tsOut = new double[numCells][];

            Parallel.For(0, numCells, c =>
            {
                int rownum = c / numCols;
                int colnum = c % numCols;
                SpatioTemporalParams geogParams = new SpatioTemporalParams();
                geogParams.Latitude = latitudes[rownum];
                geogParams.Longitude = longitudes[colnum];

                //geogParams.ModelRuntimeDays = nDays;
                //geogParams.StartJulianDay = 0;

            });
            return null;
        }
        static double[][] RunForCoarseStripe(string FilePattern, PopulationTypes ModelType)
        {
            string fnPattern = "C:\\Users\\zool1301.NDPH\\Documents\\Code_General\\temp-suitability\\Original\\c_code\\test5\\data\\{0}.{1}.dat";
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
            popParams.MosquitoLifespanDays = new TimeSpan(31,0,0,0);
            popParams.SliceLength = new TimeSpan(2, 0, 0);
            double[][] tsOut = new double[numCells][];
            //int c = 0;
            int numBlocks = 200;
            IEnumerable < int > blocks = Enumerable.Range(0, numBlocks);
            int cellsPerBlock = numCells / numBlocks;
            Parallel.ForEach(blocks, blockID =>
            {
                int blockCellStartIdx = blockID * cellsPerBlock;
                for (int c = blockCellStartIdx; (c < (blockCellStartIdx + cellsPerBlock) && c < numCells); c++)
                {
                    tsOut[c] = new double[nDays];
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
                
                    TSCellModel tsModel = new TSCellModel(popParams, geogParams, ModelType);
                    //double[] tsCell = tsModel.RunModel(tMinsForCell, tMaxsForCell);
                    //for (int i = 0; i < nDays; i++)
                    //{
                    //    tsOut[c][i] = tsCell[i];
                    //}
                }

                
            });
            return tsOut;
        }
        static void Main(string[] args)
        {

            Stopwatch sw = new Stopwatch();
            sw.Start();
            for (int i = 0; i <= 45; i += 5)
            {
                var res = RunForStripe(i.ToString(), PopulationTypes.Pointers);
               /* string outPath = String.Format("C:\\temp\\{0}.TestOutput.dat", i);
                for (int i = 0; i < res.GetLength(0); i++)
                {

                }*/
            }

            sw.Stop();
            Console.WriteLine("Time elapsed running array-based model = {0}", sw.Elapsed);
            sw.Restart();
            for (int i = 0; i<=45; i+= 5)
            {
                var res = RunForCoarseStripe(i.ToString(), PopulationTypes.Pointers);
            }
            sw.Stop();
            Console.WriteLine("Time elapse running on coarse stripes = {0}", sw.Elapsed);
            Console.ReadLine();
            //sw = new Stopwatch();
            //sw.Start();
            //for (int i = 0; i <= 45; i += 5)
            //{
            //    var res = RunForStripe(i.ToString(), PopulationTypes.OOCohorts);
            //    /* string outPath = String.Format("C:\\temp\\{0}.TestOutput.dat", i);
            //     for (int i = 0; i < res.GetLength(0); i++)
            //     {

            //     }*/
            //}
            //sw.Stop();
            //Console.WriteLine("Time elapsed running cohort-based model = {0}", sw.Elapsed);
            
            System.Console.WriteLine("Done");
        }
    }
}
