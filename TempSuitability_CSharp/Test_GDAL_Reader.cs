using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace TempSuitability_CSharp
{
    /// <summary>
    /// this is a test class only for figuring out the fastest way to read a block from a large number of tiff files and 
    /// reformat the resulting data into vertical drill-down slices
    /// </summary>
    class Test_GDAL_Reader
    {
        static List<Tuple<string, DateTime>> GetFilenamesAndDates(string filenamePath, IFilenameDateParser filenameParser)
        {
            string _dir = System.IO.Path.GetDirectoryName(filenamePath);
            string _fnPattern = System.IO.Path.GetFileName(filenamePath);
            // add in new parser classes for different file formats
            
            List<string> _files = System.IO.Directory.GetFiles(_dir, _fnPattern).ToList();
            List<Tuple<string, DateTime>> _fileDates = 
                _files.Select(fn => new Tuple<string, DateTime?>(fn, filenameParser.TryParseFilenameDate(fn)))
                .Where(tpl => tpl.Item2.HasValue)
                .Select(tpl => new Tuple<string, DateTime>(tpl.Item1, tpl.Item2.Value))
                .ToList();
            // sort the list of filename/date pairs by date
            _fileDates.Sort((x, y) => x.Item2.CompareTo(y.Item2));
            return _fileDates;        
        }
    
        static float[] TestReadTileAcrossTime(int column, int row, int xSize = 512, int ySize = 512)
        {
            int xOff = column * xSize;
            int yOff = row * ySize;
            
            string fileWildCard = "F:\\MOD11A2_Gapfilled_Output\\LST_Day\\Output_Final_30k_2030pc\\*Data.tif";
            IFilenameDateParser modisFileParse = new FilenameDateParser_MODIS8DayRaw();
            var details = GetFilenamesAndDates(fileWildCard, modisFileParse);
            string firstFileName = details[0].Item1;
            double[] overallGT = GDAL_Operations.GetGeoTransform(firstFileName);
            var shape = GDAL_Operations.GetRasterShape(firstFileName);
            
            if (yOff > shape.Item1 || xOff > shape.Item2)
            {
                throw new ArgumentException("you specified a column or row greater than the number of tiles available");
            }
            if (yOff + ySize > shape.Item1)
            {
                ySize = shape.Item1 - yOff;
            }
            if (xOff + xSize > shape.Item2)
            {
                xSize = shape.Item2 - xOff;
            }
            int nPix = xSize * ySize;
            float[] tileData = new float[nPix * details.Count];
            
            for (int t = 0; t<details.Count; t++)
            {
                int pxStart = nPix * t;
                string filename = details[t].Item1;
                DateTime filedate = details[t].Item2;
                var newshape = GDAL_Operations.GetRasterShape(filename);
                if (newshape.Item1 != shape.Item1 || newshape.Item2 != shape.Item2)
                {
                    throw new ArgumentException("Raster shapes don't match");
                }
                var bandarr = GDAL_Operations.ReadGDALRasterBandsToFlatArray(filename, xSize, ySize, xOff, yOff, 1);
                Array.Copy(bandarr, 0, tileData, pxStart, bandarr.Length);
               // tileData[pxStart : pxStart+nPix] = arr;
            }
            return tileData;
        }

        static void Main(string [] args)
        {
            
            Environment.SetEnvironmentVariable("PATH", Environment.GetEnvironmentVariable("PATH")
         + ";C:\\Users\\zool1301.NDPH\\Documents\\Code_General\\temp-suitability\\TempSuitability_CSharp\\packages\\GDAL.Native.1.11.1\\gdal\\x64");

            //string testFile = "G:\\DataPrep\\ts_global\\TempSuitability.Pf.AnnualInfectiousDays.1k.2010.global.tif";
            string testFile = "F:\\MOD11A2_Gapfilled_Output\\LST_Day\\Output_Final_30k_2030pc\\LST_Day_All.vrt";

            int xsize = 512;
            int ysize = 512;
            int nbands = 727;

            Stopwatch sw = new Stopwatch();
            sw.Start();
            //var arr = GDAL_Operations.ReadGDALRasterBandsToFlatArray(testFile, xsize, ysize, 20480, 9216, null);
            var arr = TestReadTileAcrossTime(40, 18);
            sw.Stop();
            Console.WriteLine("Time elapsed reading tif data into flat array via pointers: {0}", sw.Elapsed);
            sw.Restart();
            int nPx = xsize * ysize;

            var cellArr = new float[nPx][];
            for (int y = 0; y < ysize; y++)
            {
                for (int x = 0; x < xsize; x++)
                {
                    int cellNum = y * xsize + x;
                    cellArr[cellNum] = new float[nbands];
                    for (int z = 0; z < nbands; z++)
                    {
                        cellArr[cellNum][z] = arr[nPx * z + cellNum];
                    }
                }
            }
            Console.WriteLine("Time elapsed reformatting flat data into cell-band order: {0}", sw.Elapsed);
            arr = null;

            sw.Restart();
            var arr3 = GDAL_Operations.ReadGDALRasterBandsToFlatArray(testFile, xsize, ysize, 20480, 9216, null);
            sw.Stop();
            Console.WriteLine("Time elapsed reading vrt data into flat array via pointers: {0}", sw.Elapsed);
            sw.Restart();
            var cellArr3 = new float[nPx][];
            for (int y = 0; y < ysize; y++)
            {
                for (int x = 0; x < xsize; x++)
                {
                    int cellNum = y * xsize + x;
                    cellArr3[cellNum] = new float[nbands];
                    for (int z = 0; z < nbands; z++)
                    {
                        cellArr3[cellNum][z] = arr3[nPx * z + cellNum];
                    }
                }
            }
            Console.WriteLine("Time elapsed reformatting flat vrt data into cell-band order: {0}", sw.Elapsed);
            arr3 = null;

            sw.Restart();
            var arr2 = GDAL_Operations.ReadGDALRasterBandsToJaggedArray(testFile, xsize, ysize, 20480, 9216, null);
            sw.Stop();
            Console.WriteLine("Time elapsed reading data into jagged array via pointers: {0}", sw.Elapsed);
            sw.Restart();
            var cellArr2 = new float[xsize * ysize][];
            for (int y = 0; y < ysize; y++)
            {
                for (int x = 0; x < xsize; x++)
                {
                    int cellNum = y * xsize + x;
                    cellArr2[cellNum] = new float[nbands];
                    for (int z = 0; z < nbands; z++)
                    {
                        cellArr2[cellNum][z] = arr2[z][cellNum];
                    }
                }
            }
            Console.WriteLine("Time elapsed reformatting jagged data into cell-band order: {0}", sw.Elapsed);
          
            System.Console.ReadKey();
        }
    }
}
/*512 * 512 * 727
44 secs to read
0.7 secs to mung
*/