using System;
using System.Configuration;
using System.Reflection;
using System.Linq;

namespace TempSuitability_CSharp
{
    //https://stackoverflow.com/a/6151688
    public abstract class AppConfig : IDisposable
    {
        public static AppConfig Change(string path)
        {
            try {
                return new ChangeAppConfig(path);
            }
            catch (NullReferenceException e)
            {
                // The dirty reflection-based method fails on mono, because the mono code has the private fields named differently 
                // e.g. configSystem instead of s_configSystem
                // https://github.com/mono/mono/blob/effa4c07ba850bedbe1ff54b2a5df281c058ebcb/mcs/class/System.Configuration/System.Configuration/ConfigurationManager.cs#L48
                // So try this approach for mono
                // https://stackoverflow.com/a/39394998/4150190
                AppDomain.CurrentDomain.SetData("APP_CONFIG_FILE", path);
                System.Configuration.Configuration newConfiguration = ConfigurationManager.OpenExeConfiguration(path);
                FieldInfo configSystemField = typeof(ConfigurationManager).GetField("configSystem", BindingFlags.NonPublic | BindingFlags.Static);
                object configSystem = configSystemField.GetValue(null);
                FieldInfo cfgField = configSystem.GetType().GetField("cfg", BindingFlags.Instance | BindingFlags.NonPublic);
                cfgField.SetValue(configSystem, newConfiguration);
                return null;
            }
        }
        public abstract void Dispose();

        private class ChangeAppConfig : AppConfig
        {
            private readonly string oldConfig =
                AppDomain.CurrentDomain.GetData("APP_CONFIG_FILE").ToString();

            private bool disposedValue;

            public ChangeAppConfig(string path)
            {
                AppDomain.CurrentDomain.SetData("APP_CONFIG_FILE", path);
                ResetConfigMechanism();
            }

            public override void Dispose()
            {
                if (!disposedValue)
                {
                    AppDomain.CurrentDomain.SetData("APP_CONFIG_FILE", oldConfig);
                    ResetConfigMechanism();
                    disposedValue = true;
                }
                GC.SuppressFinalize(this);
            }

            private static void ResetConfigMechanism()
            {
                typeof(System.Configuration.ConfigurationManager)
                    .GetField("s_initState", BindingFlags.NonPublic | BindingFlags.Static)
                    .SetValue(null, 0);
                
                typeof(ConfigurationManager)
                    .GetField("s_configSystem", BindingFlags.NonPublic |
                                           BindingFlags.Static)
                    .SetValue(null, null);
                
                typeof(ConfigurationManager)
                    .Assembly.GetTypes()
                    .Where(x => x.FullName ==
                                "System.Configuration.ClientConfigPaths")
                    .First()
                    .GetField("s_current", BindingFlags.NonPublic |
                                           BindingFlags.Static)
                    .SetValue(null, null);
            }
        }
    }
    


}
