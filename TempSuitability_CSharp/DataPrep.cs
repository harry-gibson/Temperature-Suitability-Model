using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Interpolation;

namespace TempSuitability_CSharp
{
    class DataPrep
    {
        private void GetDailyTemps()
        {
            double[] tDays = new double[500];
            double[] tTemps = new double[500];
            CubicSpline minInterpolator = CubicSpline.InterpolateNatural(tDays, tTemps);
            //minInterpolator.Interpolate(blah);
        }
    }
}
