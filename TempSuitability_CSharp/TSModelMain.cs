using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace TempSuitability_CSharp
{
    class TSModelMain
    {
        static void Main(string[] args)
        {
            Environment.SetEnvironmentVariable("PATH", Environment.GetEnvironmentVariable("PATH")
      + ";C:\\Users\\zool1301.NDPH\\Documents\\Code_General\\temp-suitability\\TempSuitability_CSharp\\packages\\GDAL.Native.1.11.1\\gdal\\x64");
            string maskPath, dayPath, nightPath, outDir;
            int maskValidValue;
            try
            {
                maskPath = Properties.Settings.Default.LS_Mask_File;
                dayPath = Properties.Settings.Default.LST_Day_Files;
                nightPath = Properties.Settings.Default.LST_Night_Files;
                outDir = Properties.Settings.Default.OutputFolder;
                maskValidValue = Properties.Settings.Default.MaskValidValue;
                /* args:  "\\map-fs1.ndph.ox.ac.uk\map_data\mastergrids\Global_Masks\Land_Sea_Masks\CoastGlobal_5k.tif" 
                "\\map-fs1.ndph.ox.ac.uk\map_data\mastergrids\MODIS_Global\MOD11A2_LST\LST_Day\5km\8-Daily\*Max.tif" 
                "\\map-fs1.ndph.ox.ac.uk\map_data\mastergrids\MODIS_Global\MOD11A2_LST\LST_Night\5km\8-Daily\*Min.tif"
                "C:\temp\tsmodel"
                */
            }
            catch
            {
                //maskPath = "\\\\map-fs1.ndph.ox.ac.uk\\map_data\\mastergrids\\Global_Masks\\Land_Sea_Masks\\CoastGlobal.tiff";
                //dayPath = "F:\\MOD11A2_Gapfilled_Output\\LST_Day\\Output_Final_30k_2030pc\\*Data.tif";
                //nightPath = "F:\\MOD11A2_Gapfilled_Output\\LST_Night\\Output_Final_30k_2030pc\\*Data.tif";
                //outDir = "C:\\Temp\\TSModel";
                //maskValidValue = 1;

                maskPath = "\\\\map-fs1.ndph.ox.ac.uk\\map_data\\mastergrids\\Global_Masks\\Land_Sea_Masks\\CoastGlobal_5k.tif";
                dayPath = "G:\\tmp_local\\LST_Day_5km_8day\\*tif";
                nightPath = "G:\\tmp_local\\LST_Night_5km_8day\\*.tif";
                outDir = "C:\\Temp\\TSModel\\5km";
                maskValidValue = 1;

            }

            TSModelRunner runner = new TSModelRunner(new FilenameDateParser_MODIS8Day(), maskPath, dayPath, nightPath, outDir, maskValidValue);

            Stopwatch sw = new Stopwatch();
            sw.Start();
            //runner.RunAllTiles(-18, 52, 38, -35, 512);
            double e, w, n, s;
            int size;
            w = Properties.Settings.Default.WestLim;
            e = Properties.Settings.Default.EastLim;
            n = Properties.Settings.Default.NorthLim;
            s = Properties.Settings.Default.SouthLim;
            size = (int)Properties.Settings.Default.TileSizePx;
            runner.RunAllTiles(w, e, n, s, size);
            //runner.RunAllTiles(13,14,6,5,512);
            sw.Stop();
            Console.WriteLine("Time elapsed running model = {0}", sw.Elapsed);
            Console.ReadKey();


        }
    }
}
