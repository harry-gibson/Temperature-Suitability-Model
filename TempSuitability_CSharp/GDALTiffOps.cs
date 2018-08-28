using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OSGeo.GDAL;

namespace TempSuitability_CSharp
{
    static class GDAL_Operations
    {
        static GDAL_Operations()
        {
            Gdal.AllRegister();
        }

        public static double[] GetGeoTransform(string gdalDataset)
        {
            if (!System.IO.File.Exists(gdalDataset))
            {
                throw new ArgumentException("File does not exist!");
            }
            Dataset ds = Gdal.Open(gdalDataset, Access.GA_ReadOnly);
            if (ds == null)
            {
                throw new ArgumentException("Can't open " + gdalDataset);
            }
            double[] inGT = new double[6];
            ds.GetGeoTransform(inGT);
            ds.Dispose();
            return inGT;
        }
        public static string GetProjection(string gdalDataset)
        {
            if (!System.IO.File.Exists(gdalDataset))
            {
                throw new ArgumentException("File does not exist!");
            }
            Dataset ds = Gdal.Open(gdalDataset, Access.GA_ReadOnly);
            if (ds == null)
            {
                throw new ArgumentException("Can't open " + gdalDataset);
            }
            var proj = ds.GetProjection();
            ds.Dispose();
            return proj;
        }
        public static float? GetNoDataValue(string gdalDataset)
        {
            if (!System.IO.File.Exists(gdalDataset))
            {
                throw new ArgumentException("File does not exist!");
            }
            Dataset ds = Gdal.Open(gdalDataset, Access.GA_ReadOnly);
            if (ds == null)
            {
                throw new ArgumentException("Can't open " + gdalDataset);
            }
            var band = ds.GetRasterBand(1);
            double ndv;
            int hasndv;
            band.GetNoDataValue(out ndv, out hasndv);
            if (hasndv == 1)
            {
                return (float)ndv;
            }
            return null;
        }

        public static double[] GetClippedGeoTransform(string gdalDataset, PixelLims subsetCoords)
        {
            if (!System.IO.File.Exists(gdalDataset))
            {
                throw new ArgumentException("File does not exist!");
            }
            Dataset ds = Gdal.Open(gdalDataset, Access.GA_ReadOnly);
            if (ds == null)
            {
                throw new ArgumentException("Can't open " + gdalDataset);
            }
            double[] inGT = new double[6];
            ds.GetGeoTransform(inGT);
            var topLeftLongIn = inGT[0];
            var topLeftLatIn = inGT[3];
            var resX = inGT[1];
            var resY = inGT[5];
            if (inGT[2] != 0.0 || inGT[4] != 0)
            {
                throw new InvalidOperationException("Only datasets with zero skew parameters (i.e. those aligned north-south) are supported");
            }
            var topLeftLongOut = topLeftLongIn + subsetCoords.WestPixelCoord * resX;
            var topLeftLatOut = topLeftLatIn + subsetCoords.NorthPixelCoord * resY;
            var clippedGT = new double[] { topLeftLongOut, resX, 0.0, topLeftLatOut, 0.0, resY };
            return clippedGT;
        }

        public static PixelLims CalculatePixelCoordsOfBlock(string gdalDataset, double westDegrees, double eastDegrees, double northDegrees, double southDegrees)
        {
            Dataset ds = Gdal.Open(gdalDataset, Access.GA_ReadOnly);
            if (ds == null)
            {
                throw new ArgumentException("Can't open " + gdalDataset);
            }
            double[] inGT = new double[6];
            ds.GetGeoTransform(inGT);
            var westLimit = inGT[0];
            var northLimit = inGT[3];
            return CalculatePixelCoordsOfBlock(inGT, westDegrees, eastDegrees, northDegrees, southDegrees, westLimit, northLimit);
        }

        /// <summary>
        /// Returns the pixel coordinates of a block of the requested lat / long coordinates, within an image of the given geotransform. 
        /// If the relativeTo*Limit parameters are provided then the coordinates will be relative to this origin rather than the image's true 
        /// origin - use this if wanting to output an image with a different overall extent to the input
        /// </summary>
        /// <param name="geoTransform"></param>
        /// <param name="westDegrees"></param>
        /// <param name="eastDegrees"></param>
        /// <param name="northDegrees"></param>
        /// <param name="southDegrees"></param>
        /// <param name="relativeToWesternLimit"></param>
        /// <param name="relativeToNorthernLimit"></param>
        /// <returns></returns>
        public static PixelLims CalculatePixelCoordsOfBlock(double[] geoTransform, double westDegrees, double eastDegrees, double northDegrees, double southDegrees,
            double? relativeToWesternLimit = null, double? relativeToNorthernLimit = null)
        {
            var resX = geoTransform[1];
            var resY = geoTransform[5];
            if (eastDegrees < westDegrees || northDegrees < southDegrees)
            {
                throw new ArgumentException("East coord must be greater than West, and North greater than South");
            }
            if (!relativeToWesternLimit.HasValue)
            {
                relativeToWesternLimit = geoTransform[0];
            }
            if (!relativeToNorthernLimit.HasValue)
            {
                relativeToNorthernLimit = geoTransform[3];
            }
            var x0 = (uint)((westDegrees - relativeToWesternLimit) / resX);
            var x1 = (uint)(((eastDegrees - relativeToWesternLimit) / resX) + 0.5);
            var y0 = (uint)((relativeToNorthernLimit - northDegrees) / (-resY));
            var y1 = (uint)(((relativeToNorthernLimit - southDegrees) / (-resY)) + 0.5);
            return new PixelLims(x0, x1, y0, y1);
        }
        public static Tuple<int,int> GetRasterShape(string gdalDataset)
        {
            if (!System.IO.File.Exists(gdalDataset))
            {
                throw new ArgumentException("File does not exist!");
            }
            Dataset ds = Gdal.Open(gdalDataset, Access.GA_ReadOnly);
            if (ds == null)
            {
                throw new ArgumentException("Can't open " + gdalDataset);
            }
            Tuple<int, int> shape = new Tuple<int, int>(ds.RasterYSize, ds.RasterXSize);
            ds.Dispose();
            return shape;
        }

        /// <summary>
        /// Read a GDAL float32 raster dataset into a one-dimensional array which contains the pixel values in order (i.e. col, row, band).
        /// Either reads one band, or all bands, depending on whether bandNum is set. Returns a one-d array of dimensions[ySize*xSize*nBands]
        /// 
        /// </summary>
        /// <param name="gdalDataset"></param>
        /// <param name="xSize"></param>
        /// <param name="ySize"></param>
        /// <param name="xOff"></param>
        /// <param name="yOff"></param>
        /// <param name="bandNum"></param>
        /// <returns></returns>
        public static unsafe float[] ReadGDALRasterBandsToFlatArray(string gdalDataset, int? xSize, int? ySize, int? xOff, int? yOff, int? bandNum)
        {
            if (!System.IO.File.Exists(gdalDataset))
            {
                throw new ArgumentException("File does not exist!");
            }
            Dataset ds = Gdal.Open(gdalDataset, Access.GA_ReadOnly);
            if (ds == null)
            {
                throw new ArgumentException("Can't open " + gdalDataset);
            }
           
            if (xSize < 0 || ySize < 0 || xOff < 0 || yOff < 0 || bandNum < 1)
            {
                throw new ArgumentOutOfRangeException("Size, offset and band num must all be positive or null - bands are numbered from 1");
            }

            if (xOff == null ) {
                xOff = 0;
            }
            if (xSize == null) {
                xSize = ds.RasterXSize - xOff;
            }
            else if (xSize + xOff > ds.RasterXSize)
            {
                throw new ArgumentOutOfRangeException("XSize + XOffset cannot be larger than dataset XSize of " + ds.RasterXSize.ToString());
            }

            if (yOff == null)
            {
                yOff = 0;
            }
            if (ySize == null) {
                ySize = ds.RasterYSize - yOff;
            }
            else if (ySize + yOff > ds.RasterYSize)
            {
                throw new ArgumentOutOfRangeException("YSize + YOffset cannot be larger than dataset YSize of " + ds.RasterYSize.ToString());
            }
            int bandFrom = 1, bandTo = 1;
            if (bandNum > ds.RasterCount)
            {
                throw new ArgumentOutOfRangeException("Band num, if specified, must not be larger than the number of bands in the dataset");
            }

            float[] buffer;
            if (bandNum == null)
            {
                // read all bands
                buffer = new float[(int)xSize * (int)ySize * ds.RasterCount];
                bandFrom = 1;
                bandTo = ds.RasterCount;
            }
            else
            {
                buffer = new float[(int)xSize * (int)ySize];
                bandFrom = bandNum.Value;
                bandTo = bandNum.Value;
            }
            
            fixed(float* bufPtr = buffer)
            {
                IntPtr bufPtrManaged = (IntPtr)bufPtr;
                for (int b = bandFrom; b <= bandTo; b++)
                {
                    Band bnd = ds.GetRasterBand(b);
                    bnd.ReadRaster((int)xOff, (int)yOff, (int)xSize, (int)ySize, bufPtrManaged,
                        (int)xSize, (int)ySize, DataType.GDT_Float32, 0, 0);
                    bufPtrManaged += (int)xSize * (int)ySize * sizeof(float);
                }
            }
            ds.Dispose();
            ds = null;
            System.GC.Collect();
            return buffer;
        }

        /// <summary>
        /// returns a jagged array of dimensions[nBands][ySize*xSize]
        /// </summary>
        /// <param name="gdalDataset"></param>
        /// <param name="xSize"></param>
        /// <param name="ySize"></param>
        /// <param name="xOff"></param>
        /// <param name="yOff"></param>
        /// <param name="bandNum"></param>
        /// <returns></returns>
        public static unsafe float[][] ReadGDALRasterBandsToJaggedArray(string gdalDataset, int? xSize, int? ySize, int? xOff, int? yOff, int? bandNum)
        {
            if (!System.IO.File.Exists(gdalDataset))
            {
                throw new ArgumentException("File does not exist!");
            }
            Dataset ds = Gdal.Open(gdalDataset, Access.GA_ReadOnly);
            if (ds == null)
            {
                throw new ArgumentException("Can't open " + gdalDataset);
            }
            if (xSize < 0 || ySize < 0 || xOff < 0 || yOff < 0 || bandNum < 1)
            {
                throw new ArgumentOutOfRangeException("Size, offset and band num must all be positive or null");
            }

            if (xOff == null)
            {
                xOff = 0;
            }
            if (xSize == null)
            {
                xSize = ds.RasterXSize - xOff;
            }
            else if (xSize + xOff > ds.RasterXSize)
            {
                throw new ArgumentOutOfRangeException("XSize + XOffset cannot be larger than dataset XSize of " + ds.RasterXSize.ToString());
            }

            if (yOff == null)
            {
                yOff = 0;
            }
            if (ySize == null)
            {
                ySize = ds.RasterYSize - yOff;
            }
            else if (ySize + yOff > ds.RasterYSize)
            {
                throw new ArgumentOutOfRangeException("YSize + YOffset cannot be larger than dataset YSize of " + ds.RasterYSize.ToString());
            }
            int bandFrom = 1, bandTo = 1;
            if (bandNum > ds.RasterCount)
            {
                throw new ArgumentOutOfRangeException("Band num, if specified, must not be larger than the number of bands in the dataset");
            }

            float[][] buffer;
            if (bandNum == null)
            {
                buffer = new float[ds.RasterCount][];
                bandFrom = 1;
                bandTo = ds.RasterCount;
            }
            else
            {
                buffer = new float[1][]; //[(int)xSize * (int)ySize];
                bandFrom = bandNum.Value;
                bandTo = bandNum.Value;
            }

            for (int b = bandFrom; b<= bandTo; b++)
            {
                buffer[b-1] = new float[(int)xSize * (int)ySize];
                Band bnd = ds.GetRasterBand(b);
                fixed (float* bufPtr = buffer[b-1])
                {
                    IntPtr bufPtrManaged = (IntPtr)bufPtr;
                    bnd.ReadRaster((int)xOff, (int)yOff, (int)xSize, (int)ySize, bufPtrManaged,
                            (int)xSize, (int)ySize, DataType.GDT_Float32, 0, 0);
                }
            }
            
            ds.Dispose();
            return buffer;
        }
        
        public static bool WriteWholeTiff(string Filename, float[] Data, double[] GeoTransform, string Projection, PixelLims shape, bool CreateCompressed, float? NDV)
        {
            // we will create a tiff with tiles (blocks) but let the size of them be specified automatically
            var ds = CreateTiff(Filename, GeoTransform, Projection, shape, CreateCompressed, NDV, null);
            var creationOpts = new string[0];
            if (CreateCompressed)
            {
                creationOpts = new string[] { "COMPRESS=DEFLATE", "ZLEVEL=9", "TILED=YES", "SPARSE_OK=FALSE", "BIGTIFF=YES" };
            }
            var drv = Gdal.GetDriverByName("GTiff");
            var shapeX = shape.EastPixelCoord - shape.WestPixelCoord;
            var shapeY = shape.SouthPixelCoord - shape.NorthPixelCoord;
            if (shapeX * shapeY != Data.Length) { 
                return false;
            }
            var band = ds.GetRasterBand(1);

            band.WriteRaster(0, 0, (int)shapeX, (int)shapeY, Data, (int)shapeX, (int)shapeY, 0, 0);
            band.FlushCache();
            ds.FlushCache();
            return true;
        }
        public static bool WritePartTiff(string Filename, float[] Data, PixelLims WriteShape)
        {
            var ds = Gdal.Open(Filename, Access.GA_Update);
            var band = ds.GetRasterBand(1);
            var fileShapeX = ds.RasterXSize;
            var fileShapeY = ds.RasterYSize;
            var writeSizeX = WriteShape.EastPixelCoord - WriteShape.WestPixelCoord;
            var writeSizeY = WriteShape.SouthPixelCoord - WriteShape.NorthPixelCoord;
            if (WriteShape.SouthPixelCoord > fileShapeY || WriteShape.EastPixelCoord > fileShapeX
                || writeSizeX * writeSizeY != Data.Length)
            {
                return false;
            }
            int dataBlockSizeX, dataBlockSizeY;
            band.GetBlockSize(out dataBlockSizeX, out dataBlockSizeY);
            if (dataBlockSizeX != writeSizeX 
                || dataBlockSizeY != writeSizeY
                || WriteShape.WestPixelCoord % dataBlockSizeX != 0
                || WriteShape.NorthPixelCoord % dataBlockSizeY != 0)
            {
                // only for testing, we will allow non-whole block writing eventually
                return false;
            }
            band.WriteRaster(
                (int)WriteShape.WestPixelCoord, (int)WriteShape.NorthPixelCoord, 
                (int)writeSizeX, (int)writeSizeY,
                Data,
                (int)writeSizeX, (int)writeSizeY, 
                0, 0);
            band.FlushCache();
            ds.FlushCache();
            return true;
        }
        private static Dataset CreateTiff(string Filename, double[] GeoTransform, string Projection, PixelLims shape, bool CreateCompressed, float? NDV, int? BlockSize)
        {
            var creationOpts = new List<string>();
            if (CreateCompressed)
            {
                creationOpts = new List<string>(){ "COMPRESS=DEFLATE", "ZLEVEL=9", "TILED=YES", "SPARSE_OK=FALSE", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS" };
            }
            if (BlockSize.HasValue)
            {
                creationOpts.Add("BLOCKXSIZE=" + BlockSize.Value.ToString());
                creationOpts.Add("BLOCKYSIZE=" + BlockSize.Value.ToString());
            }
            var drv = Gdal.GetDriverByName("GTiff");
            var shapeX = shape.EastPixelCoord - shape.WestPixelCoord;
            var shapeY = shape.SouthPixelCoord - shape.NorthPixelCoord;

            var ds = drv.Create(Filename, (int)shapeX, (int)shapeY, 1, DataType.GDT_Float32, creationOpts.ToArray());
            ds.SetGeoTransform(GeoTransform);
            ds.SetProjection(Projection);
            var band = ds.GetRasterBand(1);
            if (NDV.HasValue)
            {
                band.SetNoDataValue(NDV.Value);
            }
            return ds;
        }
    }
    public class PixelLims
    {
        public uint WestPixelCoord { get; }
        public uint EastPixelCoord { get; }
        public uint NorthPixelCoord { get; }
        public uint SouthPixelCoord { get; }
        public PixelLims(uint West, uint East, uint North, uint South)
        {
            this.WestPixelCoord= West;
            this.EastPixelCoord= East;
            this.NorthPixelCoord= North;
            this.SouthPixelCoord = South;
        }
        public PixelLims(int West, int East, int North, int South)
        {
            if (East < 0 || West < 0 || North<0 || South < 0)
            {
                throw new ArgumentException("Values must be positive");
            }
            this.WestPixelCoord = (uint)West;
            this.EastPixelCoord = (uint)East;
            this.NorthPixelCoord = (uint)North;
            this.SouthPixelCoord = (uint)South;
        }
    }
    
}
