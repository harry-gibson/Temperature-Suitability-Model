using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TempSuitability_CSharp
{
    /// <summary>
    /// Basic struct to pass TS model parameters relating to population dynamics and response
    /// </summary>
    struct PopulationParams
    {
        /// <summary>
        /// The maximum lifespan of a mosquito (they may die earlier according to temperature)
        /// </summary>
        public TimeSpan MosquitoLifespanDays { get; set; }
        /// <summary>
        /// The temporal resolution with which to run the model, e.g. 2 hour slices
        /// </summary>
        public TimeSpan SliceLength { get; set; }
        public int LifespanSlices
        {
            get
            {
                return (int)(MosquitoLifespanDays.TotalDays / SliceLength.TotalDays);
            }
        }
        /// <summary>
        /// The threshold temperature above which degree days start accumulating.
        /// If temperature never exceeds this then the suitability is by definition zero
        /// </summary>
        public double MinTempThreshold { get; set; }

        /// <summary>
        /// The maximum temperature that is deemed to be survivable. If temperature never drops 
        /// below this then the suitability is by defintion zero.
        /// </summary>
        public double MosquitoDeathTemperature { get; set; }

        /// <summary>
        /// The number of accumulated degree days above which a mosquito will be 
        /// deemed to become infectious. One degree day is a temperature of one degree 
        /// above the min temp threshold maintained for one full day.
        /// </summary>
        public double DegreeDayThreshold { get; set; }

        /// <summary>
        /// The maximum value that the temperature suitability calculation can produce
        /// </summary>
        public double MaxTempSuitability { get; set; }

        /// <summary>
        /// The minimum proportion of the data series that must be valid (non-NDV) to run the model
        /// </summary>
        public double ValidDataProportion { get; set; }

        /// <summary>
        /// Which of the published survival functions should be used to calculate what proportion of the existing 
        /// population survives each model timestep. Which function is chosen determines whether the model can be 
        /// run with only temperature data or if humidity data is also required.
        /// </summary>
        public MossieMethods.SurvivalFunctions SurvivalFunction { get; set; }


    }
}
