using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Configuration;
namespace TempSuitability_CSharp
{
    class TSModelRunner
    {
        private static readonly object _lockobj = new object();
        private readonly TiffCubeReader _maxReader;
        private readonly TiffCubeReader _minReader;
        private readonly TiffCubeReader _humReader;
        private bool _humIsSynoptic;
        private readonly string _maskPath;
        private readonly string _outDir;
        private readonly float _maskValidValue;
        private readonly float _maskNDV;
        private DateTime[] _outputDates;
        private PopulationParams _DefaultPopParams;

        public TSModelRunner(IFilenameDateParser modisParser, string maskPath, string dayPath, string nightPath, string outDir, int maskValidValue)
        {
            //System.Console.WriteLine("Looking for files in " + m_FileWildCard);
            //var d = new System.Diagnostics.DefaultTraceListener();
            _maxReader = new TiffCubeReader(dayPath, modisParser);
            var nMax = _maxReader.Filenames.Count;
            System.Console.WriteLine("Looking for LST Day files in " + dayPath + 
                " - found "+ nMax.ToString());
            _minReader = new TiffCubeReader(nightPath, modisParser);
            var nMin = _minReader.Filenames.Count;
            System.Console.WriteLine("Looking for LST Night files in " + nightPath +
                " - found " + nMin.ToString());

            if (nMax == 0 || nMin == 0)
            {
                throw new ArgumentException("Can't continue without data");
                // I don't see any reason not to continue if the numbers for day and night aren't actually equal
            }
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
            _DefaultPopParams = ParseDefaultPopulationSettings();
        }

        public TSModelRunner(IFilenameDateParser modisParser, IFilenameDateParser humParser, 
            string maskPath, string dayPath, string nightPath, string humPath,
            string outDir, int maskValidValue)
            :this(modisParser, maskPath, dayPath, nightPath, outDir, maskValidValue)
        {
            if (_DefaultPopParams.RequiresHumidity)
            {
                // _humParser = humParser;
                _humReader = new TiffCubeReader(humPath, humParser);
                var nHum = _humReader.Filenames.Count;
                System.Console.WriteLine("Looking for humidity files in " + humPath +
                   " - found " + nHum.ToString());
                if (nHum == 0)
                {
                    _humReader = null;
                    _DefaultPopParams.SurvivalFunction = MossieMethods.SurvivalFunctions.Martens2;
                    System.Console.WriteLine("Humidity-based model was requested but no humidity data loaded. Using Martens2 survival function!");
                }
                else
                {
                    if (_humReader.Filedates.First().Year == _humReader.Filedates.Last().Year)
                    {
                        // there's only one year indicating that it's synoptic data. 
                        _humIsSynoptic = true;
                    }
                    else
                    {
                        _humIsSynoptic = false;
                    }
                }
            }
        }

        /// <summary>
        /// Runs TS model for all cells of all tiles required to cover the given bounding box at the given tile size.
        /// One output file will be generated for each tile * each output timestamp (e.g. month). The tiles will be 
        /// identified with a name relative to this run e.g. a suffix like *_r001_c002.tif.
        /// Any tile wholly in the sea will be skipped (no output files generated).
        /// <param name="WestDegrees"></param>
        /// <param name="EastDegrees"></param>
        /// <param name="NorthDegrees"></param>
        /// <param name="SouthDegrees"></param>
        /// <param name="TileSize"></param>
        public void RunAllTiles(ProjectedLims LatLonLims, int TileSize)
        //public void RunAllTiles(double WestDegrees, double EastDegrees, double NorthDegrees, double SouthDegrees, int TileSize)
        {
            var WestDegrees = LatLonLims.WestPixelCoord;
            var EastDegrees = LatLonLims.EastPixelCoord;
            var NorthDegrees = LatLonLims.NorthPixelCoord;
            var SouthDegrees = LatLonLims.SouthPixelCoord;

            string runDir = WestDegrees.ToString() + "W-"
                           + EastDegrees.ToString() + "E-"
                           + NorthDegrees.ToString() + "N-"
                           + SouthDegrees.ToString() + "S-"
                           + _DefaultPopParams.SurvivalFunction.ToString() +"-"
                           + TileSize.ToString() + "px";
            System.IO.Directory.CreateDirectory(System.IO.Path.Combine(_outDir, runDir));

            var globalTempGT = _maxReader.GeoTransform;
            var pixelCoordsOfWholeRun_Temp = GDAL_Operations.CalculatePixelCoordsOfBlock(globalTempGT, WestDegrees, EastDegrees, NorthDegrees, SouthDegrees);
            
            // get the real-world coordinates of each pixel along each axis, so each model will know where it is geographically
            var tempPixelLats = _maxReader.GetSubsetLatitudeCoords(pixelCoordsOfWholeRun_Temp);
            var tempPixelLons = _maxReader.GetSubsetLongitudeCoords(pixelCoordsOfWholeRun_Temp);
            var nTilesX = (int)Math.Ceiling((double)tempPixelLons.Length / TileSize);
            var nTilesY = (int)Math.Ceiling((double)tempPixelLats.Length / TileSize);

            PixelLims pixelCoordsOfWholeRun_Hum = null;
            int TileSizeHum = TileSize;
            if (_DefaultPopParams.RequiresHumidity)
            {
                var humResolution = _humReader.GeoTransform[1];
                var tempResolution = _maxReader.GeoTransform[1];
                if (humResolution % tempResolution > 0.00001*tempResolution)
                {
                    throw new ArgumentException("Temperature and humidity datasets must have resolutions that are full multiples of one another");
                }
                var humAgg = (int)(humResolution / tempResolution);

                pixelCoordsOfWholeRun_Hum = GDAL_Operations.CalculatePixelCoordsOfBlock(_humReader.GeoTransform, WestDegrees, EastDegrees, NorthDegrees, SouthDegrees);
                var humPixelLats = _humReader.GetSubsetLatitudeCoords(pixelCoordsOfWholeRun_Hum);
                var humPixelLons = _humReader.GetSubsetLongitudeCoords(pixelCoordsOfWholeRun_Hum);
                TileSizeHum = TileSize / humAgg;
            }

            System.Console.WriteLine("Initiating run of  " + (nTilesX * nTilesY).ToString() + " tiles");

            for (int tileRow = 0; tileRow < nTilesY; tileRow++)
            {
                int yOff = (int) pixelCoordsOfWholeRun_Temp.NorthPixelCoord + tileRow * TileSize;
                int yEnd = yOff + TileSize;
                int thisTileYSize = TileSize;
           
                if (yEnd > pixelCoordsOfWholeRun_Temp.SouthPixelCoord)
                {
                    // the last row might need to be less high than the others if we've reached the southern limit
                    thisTileYSize = (int) pixelCoordsOfWholeRun_Temp.SouthPixelCoord - yOff;
                }
                for (int tileCol = 0; tileCol < nTilesX; tileCol++)
                {
                    var tileNum = tileRow * nTilesX + tileCol + 1;
                    // we will save the tiles with a relative row / column number in the filename rather than specific geographic info,
                    // because each run (each given set of tiles) will be in a uniquely-named folder
                    var tilenameLocPart = "_r" + tileRow.ToString("D3") + "_c" + tileCol.ToString("D3") + ".tif";
                    if (System.IO.Directory.EnumerateFiles(System.IO.Path.Combine(_outDir, runDir), "*" + tilenameLocPart).Count() != 0)
                    {
                        System.Console.WriteLine("Tile " + tileNum.ToString() + " appears to be already done, skipping");
                        continue;
                    }

                    // Get the pixel limits of the tile relative to the temperature input files origin
                    int xOff = (int) pixelCoordsOfWholeRun_Temp.WestPixelCoord + tileCol * TileSize;
                    int xEnd = xOff + TileSize;
                    int thisTileXSize = TileSize;
                    if (xEnd > pixelCoordsOfWholeRun_Temp.EastPixelCoord)
                    {
                        // the last column might need to be less wide than the others if we've reached the eastern limit
                        thisTileXSize = (int)pixelCoordsOfWholeRun_Temp.EastPixelCoord - xOff;
                    }
                    var pixelCoordsOfTile_temp = new PixelLims(xOff, xOff + thisTileXSize, yOff, yOff + thisTileYSize);

                    // If we're using humidity data, get the pixel coordinates of the tile relative to the humidity input files origin
                    // Note that these might have slightly different geographic offset to the temp data if the tile does not coincide 
                    // with the edge of a (larger) humidity pixel
                    PixelLims pixelCoordsOfTile_hum = null;
                    if (_DefaultPopParams.RequiresHumidity)
                    {
                        int yOffHum = (int)pixelCoordsOfWholeRun_Hum.NorthPixelCoord + tileRow * TileSizeHum;
                        int yEndHum = yOffHum + TileSizeHum;
                        int thisHumTileYSize = TileSizeHum;
                        if (yEndHum > pixelCoordsOfWholeRun_Hum.SouthPixelCoord)
                        {
                            thisHumTileYSize = (int)pixelCoordsOfWholeRun_Hum.SouthPixelCoord - yOffHum;
                        }

                        int xOffHum = (int)pixelCoordsOfWholeRun_Hum.WestPixelCoord + tileCol * TileSizeHum;
                        int xEndHum = xOffHum + TileSizeHum;
                        int thisHumTileXSize = TileSizeHum;
                        if (xEndHum > pixelCoordsOfWholeRun_Hum.EastPixelCoord)
                        {
                            thisHumTileXSize = (int)pixelCoordsOfWholeRun_Hum.EastPixelCoord - xOffHum;
                        }
                        pixelCoordsOfTile_hum = new PixelLims(xOffHum, xOffHum + thisHumTileXSize, yOffHum, yOffHum + thisHumTileYSize);
                    }

                    float[][] cellRes;
                    if (!_DefaultPopParams.RequiresHumidity)
                    {
                        cellRes = RunTile(pixelCoordsOfTile_temp);
                    }
                    else
                    {
                        cellRes = RunTile(pixelCoordsOfTile_temp, pixelCoordsOfTile_hum);
                    }

                    if (cellRes.Length == 0)
                    {
                        // whole tile was in masked / sea area. Do not bother to write output.
                        System.Console.WriteLine("Tile "+tileNum.ToString() +" was wholly in sea - skipped");
                        continue;
                    }
                    
                    System.Console.WriteLine("Tile computation completed, writing output");
                    var tileRes = TransposeCellData(
                        cellRes,
                        thisTileXSize, thisTileYSize);
                    var tileGT = _maxReader.GetSubsetGeoTransform(pixelCoordsOfTile_temp);
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
                    System.Console.WriteLine("Tile " + tileNum.ToString() + " finished - wrote "+tileRes.Length.ToString()+" files");
                }
            }
        }

        private PopulationParams ParseDefaultPopulationSettings()
        {
            var set = Properties.Settings.Default;
            // get the model parameters from the default settings file
            PopulationParams popParams = new PopulationParams();
            popParams.DegreeDayThreshold = set.ParamSporogDegreeDays;
            popParams.MinTempThreshold = set.ParamMinTempThreshCelsius;
            popParams.MosquitoDeathTemperature = set.ParamDeathTempCelsius; ;
            popParams.MosquitoLifespanDays = new TimeSpan((int)set.ParamLifespanDays, 0, 0, 0); ;
            popParams.SliceLength = new TimeSpan((int)set.ParamSliceLengthHours, 0, 0);
            // the Weiss code had 33.89401 for max ts, not sure how generated. 
            // I got 34.2467 using iterative solver in excel.
            popParams.MaxTempSuitability = set.ParamMaxTSNormaliser;

            var sfStr = set.SurvivalFunction;
            if (Enum.IsDefined(typeof(MossieMethods.SurvivalFunctions), sfStr))
            {
                popParams.SurvivalFunction = (MossieMethods.SurvivalFunctions)
                    Enum.Parse(typeof(MossieMethods.SurvivalFunctions), sfStr);
            }
            else
            {
                popParams.SurvivalFunction = MossieMethods.SurvivalFunctions.Martens2;
            }
            return popParams;
        }

        public float[][] RunTile(PixelLims TemperaturePixelCoordinates)
        {
            return RunTile(TemperaturePixelCoordinates, null);
        }

        /// <summary>
        /// Runs a temperature suitability model for all pixels in a region specified by pixel limits.
        /// Output is a jagged array with one value for each cell starting at top left and then row by 
        /// row to bottom right, EXCEPT if no pixel in the tile is in a data area in which case the output 
        /// is an array with length zero.
        /// Otherwise, each value is an array with one TS value for each month of the run period, 
        /// EXCEPT if the cell is in the sea / masked area, in which case it is an array of length 0.
        /// Each cell location is done independently and this is therefore multithreaded.
        /// </summary>
        /// <param name="xOff"></param>
        /// <param name="yOff"></param>
        /// <param name="xSize"></param>
        /// <param name="ySize"></param>
        /// <returns></returns>
        public float[][] RunTile(PixelLims TemperaturePixelCoordinates, PixelLims HumidityPixelCoords = null)
        //public float[][] RunTile(int xOff, int yOff, int xSize, int ySize)
        {
            var xSize = TemperaturePixelCoordinates.XSize;
            var ySize = TemperaturePixelCoordinates.YSize;
            var xOff = (int)TemperaturePixelCoordinates.WestPixelCoord;
            var yOff = (int)TemperaturePixelCoordinates.NorthPixelCoord;

            var lsMask = GDAL_Operations.ReadGDALRasterBandsToFlatArray(
                _maskPath,
                xSize, ySize, xOff, yOff, 1);
            if (!lsMask.Any(v => v == _maskValidValue))
            {
                // this whole tile is in the sea, no need to run, return as a special case a zero length cell array
                return new float[0][];
            }
            var isHumidity = HumidityPixelCoords != null;
            
            // read athe temperature data for the tile
            var dayData = _maxReader.ReadCellData(xOff, yOff, xSize, ySize);
            var nightData = _minReader.ReadCellData(xOff, yOff, xSize, ySize);

            // get an array of the latitudes of each pixel in the y dimension of this tile
            var lats = _maxReader.GetSubsetLatitudeCoords(yOff, ySize);
            // and same for longitudes / x dimension
            var lons = _maxReader.GetSubsetLongitudeCoords(xOff, xSize);
            var dates = _maxReader.Filedates.ToArray();
            var nFiles = dates.Length;
          
            int numCells = xSize * ySize;
            Debug.Assert(nightData.Length == dayData.Length);
            Debug.Assert(dayData.Length == numCells);

            int agg = 1;
            float[][] humData = null;
            var temperatureTileGT = _maxReader.GetSubsetGeoTransform(TemperaturePixelCoordinates);
            double[] humTileGT = null;
            if (isHumidity)
            {
                var xSizeHum = HumidityPixelCoords.XSize;
                var ySizeHum = HumidityPixelCoords.YSize;
                var xOffHum = (int)HumidityPixelCoords.WestPixelCoord;
                var yOffHum = (int)HumidityPixelCoords.NorthPixelCoord;
                humData = _humReader.ReadCellData(xOffHum, yOffHum, xSizeHum, ySizeHum);
             //   var globalHumGT = _humReader.GeoTransform;
              //  var pxHumOverall = GDAL_Operations.CalculatePixelCoordsOfBlock(globalHumGT, WestDegrees, EastDegrees, NorthDegrees, SouthDegrees);
                var humLats = _humReader.GetSubsetLatitudeCoords(yOffHum, ySizeHum);
                var humLons = _humReader.GetSubsetLongitudeCoords(xOffHum, xSizeHum);
                // compare the resolutions, we will only bother looking at the x size as all our data are square pixels...
                var tempResolution = _maxReader.GeoTransform[1];
                var humResolution = _humReader.GeoTransform[1];
                if (humResolution % tempResolution > 0.00001 * tempResolution)
                {
                    throw new ArgumentException("Temperature and humidity datasets must have resolutions that are full multiples of one another");
                }
                agg = (int)(humResolution / tempResolution);
                humTileGT = _humReader.GetSubsetGeoTransform(HumidityPixelCoords);
            }

            System.Console.WriteLine("Tile data loaded, computation beginning");

            float[][] tsOut = new float[numCells][];
            DateTime[] outputDates = null;
            ParallelOptions b = new ParallelOptions();
            if (Properties.Settings.Default.MaxThreads != 0)
            {
                // set threads to 1 for easier step-through debugging or some other number to not hog the whole machine
                b.MaxDegreeOfParallelism = (int)Properties.Settings.Default.MaxThreads;
            }

            // begin parallel processing of each cell location in the tile. Each location is entirely independent 
            // so just basic automated approach is fine. 
            Parallel.For(0, numCells, b, c =>
            {
                if (lsMask[c] != 1)
                {
                    // this cell is in the sea, no need to run, return as a special case a zero length result array
                    tsOut[c] = new float[0];
                }
                else
                {
                    int rownum = c / xSize;
                    int colnum = c % xSize;
                    
                    // cut
                    // This doesn't work accurately because the "first" temperature may be midway "into" a humidity 
                    // pixel so for an aggregation factor of say 10 we may only need to go 7 temperature pixels before we 
                    // get into the second humidity pixel. 
                    //int humRow = (int)Math.Floor((double)rownum / agg);
                    //int humCol = (int)Math.Floor((double)colnum / agg);
                    //int humCell = HumidityPixelCoords.XSize * humRow + humCol;
                    // .. or do it geographically:
                    var tempLonLat = GDAL_Operations.GetLocationOfPixelCoords(temperatureTileGT, colnum, rownum);
                    var humXY = GDAL_Operations.GetPixelCoordsOfLocation(humTileGT, tempLonLat[0], tempLonLat[1]);
                    int humCell = HumidityPixelCoords.XSize * humXY[1] + humXY[0];

                    // end cut

                    GeographicCellLocation geogParams = new GeographicCellLocation();
                    geogParams.Latitude = lats[rownum];
                    geogParams.Longitude = lons[colnum];
                    
                    // run only if we've got at least 50% of the data and 100+ points
                    var nValid = Math.Min(
                        dayData[c].Count(v => v != _maxReader.NoDataValue),
                        nightData[c].Count(v => v != _minReader.NoDataValue));
                    if (nValid < nFiles / 2 || nValid < 100)
                    {
                        tsOut[c] = new float[0];
                    }
                    else
                    {
                        THSCellModel thsModel = new THSCellModel(_DefaultPopParams, geogParams);
                        //TSCellModel tsModel = new TSCellModel(popParams, geogParams, PopulationTypes.Pointers);
                        
                        if (!thsModel.SetTempData(dayData[c], nightData[c], dates, _maxReader.NoDataValue.Value, true))
                        {
                            throw new ApplicationException("Pop goes the weasel");
                        };
                        if (!thsModel.SetSynopticHumidityData(humData[humCell], _humReader.NoDataValue.Value))
                        {
                            throw new ApplicationException("Pop goes the weasel");
                        };

                        // run the entire ts model for this location
                        float[] tsCell = thsModel.RunModel();
                        int nOutputPoints = tsCell.Length;
                        tsOut[c] = tsCell;
                        lock (_lockobj) // ensure only 1 thread makes this check 
                        {
                            // calculate and record the dates of the outputs for an arbitrary one of the cell models
                            if (outputDates == null)
                            {
                                outputDates = thsModel.OutputDates;
                            }
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
