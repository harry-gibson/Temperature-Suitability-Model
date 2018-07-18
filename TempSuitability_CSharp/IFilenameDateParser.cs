using System;

namespace TempSuitability_CSharp
{
    /// <summary>
    /// Specifies an interface for classes that can take a string filename and parse a DateTime from it, if one is present
    /// </summary>
    interface IFilenameDateParser
    {
        DateTime? TryParseFilenameDate(string Filename);
    }
}
