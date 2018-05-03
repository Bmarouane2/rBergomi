using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Meta.Numerics.Matrices;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics;
using Meta.Numerics.SignalProcessing;


namespace rBergomi
{
    public class OUApproxVolterra
    {
        public Grid grid { get; set; }

        //Gaussian Simulator
        private GaussianSimulator simulator = new GaussianSimulator();
        private double H;

        //Structor
        public OUApproxVolterra( Grid grid_, double H_ = 0.07)
        {
            grid = grid_;
            H = H_;
        }


    }
}

