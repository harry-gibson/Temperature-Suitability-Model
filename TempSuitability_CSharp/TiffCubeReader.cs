using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TempSuitability_CSharp
{
    class TiffCubeReader
    {
        private string m_FileWildCard;
        private IFilenameDateParser m_FilenameDateParser;
        private SortedList<DateTime, string> m_FileDates;
        private int m_readsizeX, m_readsizeY;
        private DateTime? m_MinDateToRead, m_MaxDateToRead;
        public IList<string> Filenames
        {
            get
            {
                return m_FileDates.Values;
            }
        }
        public IList<DateTime> Filedates
        {
            get
            {
                return m_FileDates.Keys;
            }
        }
        public double[] GlobalLatCoords
        {
            get; private set;
        }
        public double[] GlobalLonCoords
        {
            get;private set;
        }

        
        /// <summary>
        /// The geotransform is set from the first image found, and an error will occur on read if any 
        /// subsequent images do not match this.
        /// </summary>
        public double[] GeoTransform { get; private set; }
        /// <summary>
        /// The shape is set from the first image found, and an error will occur on read if any 
        /// subsequent images do not match this.
        /// </summary>
        public Tuple<int, int > Shape { get; private set; }
        /// <summary>
        /// The projection is set from the first image found, and an error will occur on read if any 
        /// subsequent images do not match this.
        /// </summary>
        public string Projection { get; private set; }

        /// <summary>
        /// The no-data value is set (or set to null) from the first image found, and an error will
        /// occur on read if any subsequent images do not match this.
        /// </summary>
        public float? NoDataValue { get; private set; }

        public TiffCubeReader(string FilepathWildcard, IFilenameDateParser FilenameDateParserObject, DateTime? MinDate, DateTime? MaxDate)
        {
            m_FileWildCard = FilepathWildcard;
            m_FilenameDateParser = FilenameDateParserObject;
            m_MinDateToRead = MinDate;
            m_MaxDateToRead = MaxDate;
        
            PopulateFilenamesAndDates();
            if (m_FileDates.Count > 0)
            {
                PopulateRasterProperties();
            }
        }
        
        public double[] GetSubsetLongitudeCoords(int xOffset, int xSize)
        {
            if (xOffset > Shape.Item2)
            {
                throw new ArgumentException("An X offset was specified that's larger than the file width");
            }
            if (xOffset + xSize > Shape.Item2)
            {
                throw new ArgumentException("An X offset + size was specified that's larger than the file width");
                //xSize = Shape.Item2 - xOffset;
            }
            var outArr = new double[xSize];
            Array.Copy(GlobalLonCoords, xOffset, outArr, 0, xSize);
            return outArr;
        }
        public double[] GetSubsetLatitudeCoords(int yOffset, int ySize)
        {
            if (yOffset > Shape.Item1)
            {
                throw new ArgumentException("A Y offset was specified that's larger than the file width");
            }
            if (yOffset + ySize > Shape.Item1)
            {
                throw new ArgumentException("A Y offset + size was specified that's larger than the file width");
                //xSize = Shape.Item2 - xOffset;
            }
            var outArr = new double[ySize];
            Array.Copy(GlobalLatCoords, yOffset, outArr, 0, ySize);
            return outArr;
        }
        public double[] GetSubsetGeoTransform(PixelLims TileCoords)
        {
            string firstFilename = m_FileDates.First().Value;
            var res = GDAL_Operations.GetClippedGeoTransform(firstFilename, TileCoords);
            return res;
        }
        /// <summary>
        /// Read the data for all cells within the specified pixel offsets, across the first band of all the rasters.
        /// Returns a jagged array with one element for each cell, each of which contains the all of the values for 
        /// that cell (location) across time. The cells are in the order given by iterating through the values of 
        /// Longitude within the values of Latitude (i.e. latitude on outer loop, northwestern cell first, through 
        /// all cells of that latitude, then moving down until south eastern cell is reached).
        /// The values at each cell correspond in order to Filedates.
        /// </summary>
        /// <param name="xOffset"></param>
        /// <param name="yOffset"></param>
        /// <param name="xSize"></param>
        /// <param name="ySize"></param>
        /// <returns></returns>
        public float[][] ReadCellData(int xOffset, int yOffset, int xSize, int ySize)
        {
            if (xOffset > Shape.Item2 || yOffset > Shape.Item1)
            {
                throw new ArgumentException("An X or Y offset was specified that's larger than the file size");
            }
            if (yOffset + ySize > Shape.Item1)
            {
                //ySize = Shape.Item1 - yOffset;
                throw new ArgumentException("A Y offset + size was specified that's larger than the file size");
            }
            if (xOffset + xSize > Shape.Item2)
            {
                throw new ArgumentException("An X offset + size was specified that's larger than the file size");
                //xSize = Shape.Item2 - xOffset;
            }
            m_readsizeX = xSize;
            m_readsizeY = ySize;
            
            int npx = m_readsizeX * m_readsizeY;
            int nbands = m_FileDates.Count;

            var cellArr = new float[npx][];

            if (((long)npx * nbands * 4) < (2 ^ 31))
            {
                float[] tileFlatArr = ReadRegionAcrossFiles_Flat(xOffset, yOffset, xSize, ySize);
                int cellNum;
                for (int y = 0; y < ySize; y++)
                {
                    for (int x = 0; x < xSize; x++)
                    {
                        cellNum = y * xSize + x;
                        cellArr[cellNum] = new float[nbands];
                        for (int z = 0; z < nbands; z++)
                        {
                            cellArr[cellNum][z] = tileFlatArr[npx * z + cellNum];
                        }
                    }
                }

            }
            else
            {
                float[][] tileFlatArr = ReadRegionAcrossFiles_MP(xOffset, yOffset, xSize, ySize);
                int cellNum;
                for (int y = 0; y < ySize; y++)
                {
                    for (int x = 0; x < xSize; x++)
                    {
                        cellNum = y * xSize + x;
                        cellArr[cellNum] = new float[nbands];
                        for (int z = 0; z < nbands; z++)
                        {
                            cellArr[cellNum][z] = tileFlatArr[z][cellNum];
                        }
                    }
                }

            }

            return cellArr;
        }

        /// <summary>
        /// Read the values within the specified pixel offsets from the first band of all this TiffCubeReader's rasters,
        /// returning them as a 1-D array in X,Y,Z order (i.e. leftmost pixel of first row of first file to rightmost 
        /// pixel of last row of last file)
        /// </summary>
        /// <param name="xOffset"></param>
        /// <param name="yOffset"></param>
        /// <param name="xSize"></param>
        /// <param name="ySize"></param>
        /// <returns></returns>
        public float[] ReadRegionAcrossFiles_Flat(int xOffset, int yOffset, int xSize, int ySize)
        {
           // var f = m_FileDates.Values;
            int nPixPerTile = xSize * ySize;
            float[] tileData = new float[nPixPerTile * m_FileDates.Count];
            int fileNum = 0;
            foreach (var t in m_FileDates)
            {
                int pxStart = nPixPerTile * fileNum;
                var newshape = GDAL_Operations.GetRasterShape(t.Value);
                if (newshape.Item1 != Shape.Item1 || newshape.Item2 != Shape.Item2)
                {
                    throw new ArgumentException("Raster shapes don't match");
                }
                var newGT = GDAL_Operations.GetGeoTransform(t.Value);
                if (!GeoTransform.SequenceEqual(newGT))
                {
                    throw new ArgumentException("Raster geotransforms don't match");
                }
                var newNDV = GDAL_Operations.GetNoDataValue(t.Value);
                if (newNDV != NoDataValue)
                {
                    throw new ArgumentException("Raster nodata values don't match");
                }
                var newProj = GDAL_Operations.GetProjection(t.Value);
                if (newProj != Projection)
                {
                    throw new ArgumentException("Raster projections don't match");
                }
                var bandarr = GDAL_Operations.ReadGDALRasterBandsToFlatArray(
                    t.Value, xSize, ySize, xOffset, yOffset, 1);
                Array.Copy(bandarr, 0, tileData, pxStart, bandarr.Length);
                fileNum += 1;
            }
            return tileData;
        }

        public float[][] ReadRegionAcrossFiles_MP(int xOffset, int yOffset, int xSize, int ySize)
        {
            int nPixPerTile = xSize * ySize;
            float[][] tileData = new float[m_FileDates.Count][];
            ParallelOptions b = new ParallelOptions();
            b.MaxDegreeOfParallelism = 6;
            var keys = m_FileDates.Keys;
            Parallel.For(0, keys.Count, b, c =>
            {
                var fDate = keys[c];
                var fName = m_FileDates[fDate];
                var newshape = GDAL_Operations.GetRasterShape(fName);
                if (newshape.Item1 != Shape.Item1 || newshape.Item2 != Shape.Item2)
                {
                    throw new ArgumentException("Raster shapes don't match");
                }
                var newGT = GDAL_Operations.GetGeoTransform(fName);
                if (!GeoTransform.SequenceEqual(newGT))
                {
                    throw new ArgumentException("Raster geotransforms don't match");
                }
                var newNDV = GDAL_Operations.GetNoDataValue(fName);
                if (newNDV != NoDataValue)
                {
                    throw new ArgumentException("Raster nodata values don't match");
                }
                var newProj = GDAL_Operations.GetProjection(fName);
                if (newProj != Projection)
                {
                    throw new ArgumentException("Raster projections don't match");
                }
                tileData[c] = GDAL_Operations.ReadGDALRasterBandsToFlatArray(
                    fName, xSize, ySize, xOffset, yOffset, 1);

            });
            return tileData;
        }

        public float[][]ReadRegionAcrossFiles(int xOffset, int yOffset, int xSize, int ySize)
        {
            int nPixPerTile = xSize * ySize;
            float[][] tileData = new float[m_FileDates.Count][];
            int fileNum = 0;
            foreach (var t in m_FileDates)
            {
                var newshape = GDAL_Operations.GetRasterShape(t.Value);
                if (newshape.Item1 != Shape.Item1 || newshape.Item2 != Shape.Item2)
                {
                    throw new ArgumentException("Raster shapes don't match");
                }
                var newGT = GDAL_Operations.GetGeoTransform(t.Value);
                if (!GeoTransform.SequenceEqual(newGT))
                {
                    throw new ArgumentException("Raster geotransforms don't match");
                }
                var newNDV = GDAL_Operations.GetNoDataValue(t.Value);
                if (newNDV != NoDataValue)
                {
                    throw new ArgumentException("Raster nodata values don't match");
                }
                var newProj = GDAL_Operations.GetProjection(t.Value);
                if (newProj != Projection)
                {
                    throw new ArgumentException("Raster projections don't match");
                }
                tileData[fileNum] = GDAL_Operations.ReadGDALRasterBandsToFlatArray(
                    t.Value, xSize, ySize, xOffset, yOffset, 1);
                fileNum += 1;
            }
            return tileData;
        }

        /// <summary>
        /// Parse the dates for all files matching the given wildcard, filter to only include those within the 
        ///  requested date range if one has been set on this reader, and store the dates against the 
        /// filenames in a sorted (by date) list of key value pairs date:filename
        /// </summary>
        private void PopulateFilenamesAndDates()
        {
            string _dir = System.IO.Path.GetDirectoryName(m_FileWildCard);
            string _fnPattern = System.IO.Path.GetFileName(m_FileWildCard);
         
            string[] _files = System.IO.Directory.GetFiles(_dir, _fnPattern);
            m_FileDates = new SortedList<DateTime,string> (_files
                .Select(fn => new KeyValuePair<DateTime?, string>(m_FilenameDateParser.TryParseFilenameDate(fn), fn))
                .Where(kvp => kvp.Key.HasValue)
                .Select(kvp => new KeyValuePair<DateTime, string>(kvp.Key.Value, kvp.Value))
                .Where(kvp => !m_MinDateToRead.HasValue || kvp.Key >= m_MinDateToRead)
                .Where(kvp => !m_MaxDateToRead.HasValue || kvp.Key <= m_MaxDateToRead)
                .ToDictionary(kvp => kvp.Key, kvp=> kvp.Value)
                );
        }

        /// <summary>
        /// Parse key information about the first raster: the geotransform, shape (n pixels X and Y), and 
        /// nodata value (of the first band). Also calculate the latitude of each row and the longitude of 
        /// each column.
        /// </summary>
        private void PopulateRasterProperties()
        {
            string firstFilename = m_FileDates.First().Value;
            GeoTransform = GDAL_Operations.GetGeoTransform(firstFilename);
            Shape = GDAL_Operations.GetRasterShape(firstFilename);
            NoDataValue = GDAL_Operations.GetNoDataValue(firstFilename);
            Projection = GDAL_Operations.GetProjection(firstFilename);
            double[] latcoords = new double[Shape.Item1];
            double[] loncoords = new double[Shape.Item2];
            double originX = GeoTransform[0];
            double cellsizeX = GeoTransform[1];
            double originY = GeoTransform[3];
            double cellsizeY = GeoTransform[5];
            for (int i = 0; i < Shape.Item1; i++)
            {
                latcoords[i] = originY + i * cellsizeY;
            }
            for (int i = 0; i < Shape.Item2; i++)
            {
                loncoords[i] = originX + i * cellsizeX;
            }
            GlobalLatCoords = latcoords;
            GlobalLonCoords = loncoords;
        }
    }
}
