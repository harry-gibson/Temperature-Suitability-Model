using System;

namespace TempSuitability_CSharp
{
    interface IFilenameDateParser
    {
        DateTime? TryParseFilenameDate(string Filename);
    }
}
