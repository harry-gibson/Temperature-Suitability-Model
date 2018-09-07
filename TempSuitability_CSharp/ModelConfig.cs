using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Serialization;
using System.IO;

namespace TempSuitability_CSharp
{
    [XmlRoot("TSI_Model_Configuration")]
    public class TSIModelConfig
    {
        [XmlElement("Population_Model_Parameters")]
        public ModelConfig modelConfig { get; set; }
        [XmlElement("Data_Paths")]
        public DataPathConfig dataPathConfig { get; set; }
        [XmlElement("Spatial_Extent_Degrees")]
        public SpatialLimits spatialLimits { get; set; }
        [XmlElement("Execution_Parameters")]
        public ModelRunConfig modelRunConfig { get; set; }

        static public void SaveConfig(TSIModelConfig config, string pathFolder, string fileName)
        {
            string pathFileName;
            string root = Path.GetPathRoot(pathFolder);
            if (root != pathFolder)
            {
                pathFileName = @pathFolder + Path.DirectorySeparatorChar + fileName + ".xml";
            }
            else
            {
                pathFileName = @pathFolder + fileName + ".xml";
            }
            XmlSerializer serializer = new XmlSerializer(typeof(TSIModelConfig));

            using (TextWriter writer = new StreamWriter(pathFileName))
            {
                serializer.Serialize(writer, config);
            }
        }
        static public TSIModelConfig LoadOrCreateConfig(string path = null)
        {
            if (path == null)
            {
                return LoadConfigFromSettings();
            }
            try
            {
                if (File.Exists(@path))
                {
                    XmlSerializer serializer = new XmlSerializer(typeof(TSIModelConfig));
                    using (FileStream fileStream = new FileStream(@path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                    {
                        var modelConfig = (TSIModelConfig)serializer.Deserialize(fileStream);
                        return modelConfig;
                    }
                }
                else if (Directory.Exists(@path))
                {
                    Console.WriteLine("Folder was provided, attempting to create sample config from application settings");
                    var modelConfig = LoadConfigFromSettings();
                    SaveConfig(modelConfig, path, "Sample_Config");
                    Console.WriteLine("Sample_Config.xml has been created in the specified folder, "+
                        "please edit this as appropriate and re-run specifying this file");
                    return null;
                }
                else
                {
                    throw new FileNotFoundException(@path + " does not exist as a file or directory");
                }
                
            }
            catch (Exception e)
            {
                Console.WriteLine("Error occurred reading the specified config file - cannot continue");
                Console.WriteLine(e.Message);
                return null;
            }
        }

        static public TSIModelConfig LoadConfigFromSettings()
        {
            var set = Properties.Settings.Default;
            DataPathConfig dpc = new DataPathConfig();
            try
            {
                dpc.MaskFile = set.LS_Mask_File;
                dpc.MaxTempFiles = set.Max_Temp_Files;
                dpc.MinTempFiles = set.Min_Temp_Files;
                dpc.OutputFolder = set.OutputFolder;
            }
            catch
            {
                Console.WriteLine("Could not find paths specification in application config file");
                return null;
            }
            SpatialLimits spc = new SpatialLimits();
            try
            {
                spc.WestLimitDegrees = set.WestLim;
                spc.EastLimitDegrees = set.EastLim;
                spc.WestLimitDegrees = set.NorthLim;
                spc.SouthLimitDegrees = set.SouthLim;
            }
            catch
            {
                Console.WriteLine("Failed to read geographic limits from application config file");
                return null;
            }
            ModelConfig mcfg = new ModelConfig();
            try
            {
                mcfg.SliceLengthHours = set.ParamSliceLengthHours;
                mcfg.LifespanDays = set.ParamLifespanDays;
                mcfg.DeathTempCelsius = set.ParamDeathTempCelsius;
                mcfg.SporogenesisDegreeDays = set.ParamSporogDegreeDays;
                mcfg.MinTempThresholdCelsius = set.ParamMinTempThreshCelsius;
                mcfg.MaxTSNormaliser = set.ParamMaxTSNormaliser;
            }
            catch
            {
                Console.WriteLine("Failed to read temp suitability parameters from application config file");
                return null;
            }
            ModelRunConfig rcfg = new ModelRunConfig();
            try
            {
                rcfg.OutputFileTag = set.Output_File_Tag;
                rcfg.MaxTileSizePx = set.TileSizePx;
                rcfg.MaxThreads = set.MaxThreads;
                rcfg.MaxTempFilesAreLST = set.Max_Temp_Convert_From_LST;
                rcfg.MinTempFilesAreLST = set.Min_Temp_Convert_From_LST;
                rcfg.MinRequiredDataPoints = set.Min_Required_Data_Points;
                rcfg.MaskValidValue = set.MaskValidValue;
                rcfg.ReadFromDate = set.Read_From_Date;
                rcfg.ReadToDate = set.Read_To_Date;
            }
            catch
            {
                Console.WriteLine("Failed to read model execution parameters from application config file");
                return null;
            }
            var tsiCfg = new TSIModelConfig();
            tsiCfg.modelConfig = mcfg;
            tsiCfg.spatialLimits = spc;
            tsiCfg.dataPathConfig = dpc;
            tsiCfg.modelRunConfig = rcfg;
            return tsiCfg;
        }
    }

    // TODO this is redundant in PopulationParams so get rid
    public class ModelConfig
    {
        public int SliceLengthHours { get; set; }
        public int LifespanDays { get; set; }
        public double DeathTempCelsius { get; set; }
        public double SporogenesisDegreeDays { get; set; }
        public double MinTempThresholdCelsius { get; set; }
        public double MaxTSNormaliser { get; set; }
    }

    public class DataPathConfig
    {
        public string MaxTempFiles { get; set; }
        public string MinTempFiles { get; set; }
        public string MaskFile { get; set; }
        public string OutputFolder { get; set; }
    }

    public class SpatialLimits
    {
        public double WestLimitDegrees { get; set; }
        public double EastLimitDegrees { get; set; }
        public double NorthLimitDegrees { get; set; }
        public double SouthLimitDegrees { get; set; }
    }

    public class ModelRunConfig
    {
        public string OutputFileTag { get; set; }
        public ushort MaxTileSizePx { get; set; }
        public ushort MaxThreads { get; set; }
        public bool MaxTempFilesAreLST { get; set; }
        public bool MinTempFilesAreLST { get; set; }
        public int MinRequiredDataPoints { get; set; }
        public int MaskValidValue { get; set; }
        public DateTime ReadFromDate { get; set; }
        public DateTime ReadToDate { get; set; }
    }

}