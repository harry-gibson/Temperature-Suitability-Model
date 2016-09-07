﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace TempSuitability_CSharp
{
    class TSModelRunner
    {
        private static readonly object _lockobj = new object();
        private readonly TiffCubeReader _maxReader;
        private readonly TiffCubeReader _minReader;
        private readonly IFilenameDateParser _fnParser;
        private readonly string _maskPath;
        private readonly string _outDir;
        private readonly float _maskValidValue;
        private readonly float _maskNDV;
        private DateTime[] _outputDates;
        static void Main(string[] args)
        {
            Environment.SetEnvironmentVariable("PATH", Environment.GetEnvironmentVariable("PATH")
      + ";C:\\Users\\zool1301.NDPH\\Documents\\Code_General\\temp-suitability\\TempSuitability_CSharp\\packages\\GDAL.Native.1.11.1\\gdal\\x64");
            string maskPath, dayPath, nightPath, outDir;
            int maskValidValue;
            try {
                maskPath = args[0];
                dayPath = args[1];
                nightPath = args[2];
                outDir = args[3];
                maskValidValue = Convert.ToInt32(args[4]);
                /* "\\map-fs1.ndph.ox.ac.uk\map_data\mastergrids\Global_Masks\Land_Sea_Masks\CoastGlobal_5k.tif" 
                "\\map-fs1.ndph.ox.ac.uk\map_data\mastergrids\MODIS_Global\MOD11A2_LST\LST_Day\5km\8-Daily\*Max.tif" 
                "\\map-fs1.ndph.ox.ac.uk\map_data\mastergrids\MODIS_Global\MOD11A2_LST\LST_Night\5km\8-Daily\*Min.tif"
                "C:\temp\tsmodel"
                */
            }
            catch
            {
                maskPath = "\\\\map-fs1.ndph.ox.ac.uk\\map_data\\mastergrids\\Global_Masks\\Land_Sea_Masks\\CoastGlobal.tiff";
                dayPath = "F:\\MOD11A2_Gapfilled_Output\\LST_Day\\Output_Final_30k_2030pc\\*Data.tif";
                nightPath = "F:\\MOD11A2_Gapfilled_Output\\LST_Night\\Output_Final_30k_2030pc\\*Data.tif";
                outDir = "C:\\Temp\\TSModel";
                maskValidValue = 1;
            }

            TSModelRunner runner = new TSModelRunner(new FilenameDateParser_MODIS8Day(), maskPath, dayPath, nightPath, outDir, maskValidValue);

            Stopwatch sw = new Stopwatch();
            sw.Start();
            //runner.RunAllTiles(-18, 52, 38, -35, 512);
            runner.RunAllTiles(8,14,6,0,512);
            sw.Stop();
            Console.WriteLine("Time elapsed running model = {0}", sw.Elapsed);
            Console.ReadKey();
          
        }
        private TSModelRunner(IFilenameDateParser _parser, string maskPath, string dayPath, string nightPath, string outDir, int maskValidValue)
        {
            _fnParser = _parser;
            _maxReader = new TiffCubeReader(dayPath, _fnParser);
            _minReader = new TiffCubeReader(nightPath, _fnParser);
            if (!_maxReader.GeoTransform.SequenceEqual(_minReader.GeoTransform))
            {
                throw new ArgumentException("Day and night image transforms do not match!");
            }
            if (!GDAL_Operations.GetGeoTransform(maskPath).SequenceEqual(_maxReader.GeoTransform))
            {
                throw new ArgumentException("Land-sea mask doesn't match images!");
            }
            _maskPath = maskPath;
            var ndv = GDAL_Operations.GetNoDataValue(maskPath);
            if(ndv.HasValue)
            {
                _maskNDV = ndv.Value;
            }
            else
            {
                _maskNDV = float.MinValue;
            }
            if (!System.IO.Directory.Exists(outDir))
            {
                System.IO.Directory.CreateDirectory(outDir);
            }
            _outDir = outDir;
            _maskValidValue = maskValidValue;
        }
        public void RunAllTiles(double WestDegrees, double EastDegrees, double NorthDegrees, double SouthDegrees, int TileSize)
        {
            var globalGT = _maxReader.GeoTransform;
            var pxOverall = GDAL_Operations.CalculatePixelCoordsOfBlock(globalGT, WestDegrees, EastDegrees, NorthDegrees, SouthDegrees);
            var latsToRun = _maxReader.GetSubsetLatitudeCoords((int)pxOverall.NorthPixelCoord, (int)(pxOverall.SouthPixelCoord - pxOverall.NorthPixelCoord));
            var lonsToRun = _maxReader.GetSubsetLongitudeCoords((int)pxOverall.WestPixelCoord, (int)(pxOverall.EastPixelCoord - pxOverall.WestPixelCoord));
            var nTilesX = (int)Math.Ceiling((double)lonsToRun.Length / TileSize);
            var nTilesY = (int)Math.Ceiling((double)latsToRun.Length / TileSize);
            string runDir = WestDegrees.ToString() + "W-"
                            + EastDegrees.ToString() + "E-"
                            + NorthDegrees.ToString() + "N-"
                            + SouthDegrees.ToString() + "S";
            //if (!System.IO.Directory.Exists(System.IO.Path.Combine(_outDir, runDir)))
            //{
            System.IO.Directory.CreateDirectory(System.IO.Path.Combine(_outDir, runDir));
            //}
            for (int tileRow = 0; tileRow < nTilesY; tileRow++)
            {
                int yOff = (int) pxOverall.NorthPixelCoord + tileRow * TileSize;
                int yEnd = yOff + TileSize;
                int thisTileYSize = TileSize;
                if (yEnd > pxOverall.SouthPixelCoord)
                {
                    thisTileYSize = (int) pxOverall.SouthPixelCoord - yOff;
                }
                for (int tileCol = 0; tileCol < nTilesX; tileCol++)
                {
                    int xOff = (int) pxOverall.WestPixelCoord + tileCol * TileSize;
                    int xEnd = xOff + TileSize;
                    int thisTileXSize = TileSize;
                    if (xEnd > pxOverall.EastPixelCoord)
                    {
                        thisTileXSize = (int)pxOverall.EastPixelCoord - xOff;
                    }

                    var tileCoords = new PixelLims(xOff, xOff + thisTileXSize, yOff, yOff + thisTileYSize);
                    var cellRes = RunTile(xOff, yOff, thisTileXSize, thisTileYSize);
                    if (cellRes.Length == 0)
                    {
                        // whole tile was in masked / sea area. Do not bother to write output.
                        continue;
                    }
                    var tileRes = TransposeCellData(
                        cellRes,
                        thisTileXSize, thisTileYSize);
                    var tileGT = _maxReader.GetSubsetGeoTransform(tileCoords);
                    var tileProj = _maxReader.Projection;
                    // write each tile to a separate tiff file - we can mosaic them later.
                    var tileLims = new PixelLims(0, thisTileXSize, 0, thisTileYSize);
                    for (int t = 0; t<tileRes.Length; t++)
                    {
                        var tData = tileRes[t];
                        string fn = _outputDates[t].Date.ToString("yyyyMMdd") + "_r" + tileRow.ToString("D3") + "_c" + tileCol.ToString("D3") + ".tif";
                        string outFile = System.IO.Path.Combine(_outDir, runDir, fn);
                        GDAL_Operations.WriteWholeTiff(outFile, tData, tileGT, tileProj, tileLims, true, _maxReader.NoDataValue);
                    }
                }
            }
        }
        
        /// <summary>
        /// Runs a temperature suitability model for all pixels in a region specified by pixel limits.
        /// Output is a jagged array with one value for each cell starting at top left and then row by 
        /// row to bottom right, EXCEPT if no pixel in the tile is in a data area in which case the output 
        /// is an array with length zero.
        /// Otherwise, each value is an array with one TS value for each month of the run period, 
        /// EXCEPT if the cell is in the sea / masked area, in which case it is an array of length 0.
        /// </summary>
        /// <param name="xOff"></param>
        /// <param name="yOff"></param>
        /// <param name="xSize"></param>
        /// <param name="ySize"></param>
        /// <returns></returns>
        public float[][] RunTile(int xOff, int yOff, int xSize, int ySize)
        {
            // test area: int xOff=20480, int yOff=9216, int xSize=512, int ySize=512
            var lsMask = GDAL_Operations.ReadGDALRasterBandsToFlatArray(
                _maskPath,
                xSize, ySize, xOff, yOff, 1);
            if (!lsMask.Any(v => v == _maskValidValue))
            {
                // this whole tile is in the sea, no need to run, return as a special case a zero length cell array
                return new float[0][];
            }
            var dayData = _maxReader.ReadCellData(xOff, yOff, xSize, ySize);
            var lats = _maxReader.GetSubsetLatitudeCoords(yOff, ySize);
            var lons = _maxReader.GetSubsetLongitudeCoords(xOff, xSize);
            var dates = _maxReader.Filedates.ToArray();
            var nFiles = dates.Length;
            var nightData = _minReader.ReadCellData(xOff, yOff, xSize, ySize);
            int numCells = xSize * ySize;
            Debug.Assert(nightData.Length == dayData.Length);
            Debug.Assert(dayData.Length == numCells);
            
            PopulationParams popParams = new PopulationParams();
            popParams.DegreeDayThreshold = 111;
            popParams.MinTempThreshold = 16;
            popParams.MosquitoDeathTemperature = 40.0;
            popParams.MosquitoLifespanDays = new TimeSpan(31, 0, 0, 0); ;
            popParams.SliceLength = new TimeSpan(2, 0, 0);
            popParams.MaxTempSuitability = 34.2467; // 33.89401 was in Weiss code, not sure how generated. I got this using iterative solver.

            float[][] tsOut = new float[numCells][];
            DateTime[] outputDates = null;
            ParallelOptions b = new ParallelOptions();
          //  b.MaxDegreeOfParallelism = 1;

            Parallel.For(0, numCells, b, c =>
            {
                if (lsMask[c] != 1)
                {
                    // this cell is in the sea, no need to run, return as a special case a zero length result array
                    tsOut[c] = new float[0];
                }
                int rownum = c / xSize;
                int colnum = c % xSize;
                SpatioTemporalParams geogParams = new SpatioTemporalParams();
                geogParams.Latitude = lats[rownum];
                geogParams.Longitude = lons[colnum];

                //geogParams.ModelRuntimeDays = nDays;
                //geogParams.StartJulianDay = 0;
                if (dayData[c].Count(v => v != _maxReader.NoDataValue) < nFiles / 2)
                {
                    tsOut[c] = new float[0];
                }
                else
                {
                    TSCellModel tsModel = new TSCellModel(popParams, geogParams, PopulationTypes.Pointers);
                    if (!tsModel.SetData(dayData[c], nightData[c], dates, _maxReader.NoDataValue.Value))
                    {
                        throw new ApplicationException("Pop goes the weasel");
                    };
                    // run the entire ts model for this location
                    float[] tsCell = tsModel.RunModel();
                    int nOutputPoints = tsCell.Length;
                    tsOut[c] = new float[nOutputPoints];
                    for (int i = 0; i < nOutputPoints; i++)
                    {
                        tsOut[c][i] = tsCell[i];
                    }
                    lock (_lockobj) // ensure only 1 thread makes this check 
                    {
                        // calculate and record the dates of the outputs for an arbitrary one of the cell models
                        if (outputDates == null)
                        {
                            outputDates = tsModel.OutputDates;
                        }
                    }
                }
            }
            );
            if (_outputDates == null)
            {
                _outputDates = outputDates;
            }
            else if (!_outputDates.SequenceEqual(outputDates))
            {
                throw new ApplicationException("Dates have changed between tiles somehow, this shouldn't happen...");
            }
            return tsOut;
   //         Console.ReadKey();
        }

        /// <summary>
        /// Transposes a jagged array - from cell * time, to time * cell.
        /// Each value of the output therefore represents the cell values across a tile 
        /// at a particular time point, and can be written to a tiff file of that time.
        /// </summary>
        /// <param name="CellData"></param>
        /// <param name="xSize"></param>
        /// <param name="ySize"></param>
        /// <returns></returns>
        private float[][] TransposeCellData(float[][] CellData, int xSize, int ySize)
        {
            var nCells = xSize * ySize;
            if (nCells != CellData.Length)
            {
                throw new ArgumentException("Specified tile shape does not match the given number of cells");
            }
            var nTiles = CellData.First(cd => cd.Length != 0).Length;
            var tileArr = new float[nTiles][];
            var ndv = _maxReader.NoDataValue.Value;
            for (int b = 0; b < nTiles; b++)
            {
                tileArr[b] = new float[nCells];
            }
            int cellNum;
            for (int y = 0; y < ySize; y++)
            {
                for (int x = 0; x < xSize; x++)
                {
                    cellNum = y * xSize + x;
                    if (CellData[cellNum].Length == 0)
                    {
                        // this cell was in the sea, indicated by the special case result that the cell data is of length 0
                        // Write nodata for this location into all tiles
                        for (int b = 0; b < nTiles; b++)
                        {
                            tileArr[b][cellNum] = ndv;
                        }
                    }
                    else
                    {
                        for (int b = 0; b < nTiles; b++)
                        {
                            tileArr[b][cellNum] = CellData[cellNum][b];
                        }
                    }
                }
            }
            return tileArr;
        }
    }
}