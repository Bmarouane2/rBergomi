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
    class OUApproxVolterra
    {
        //test
        public Grid grid { get; set; }//grid of time Discretisation
        private GaussianSimulator simulator = new GaussianSimulator();
        private double alpha;
        private double r;
        private int N;//number of factors in the approximation
        //private int M; 

        private ColumnVector C;
        private ColumnVector Gamma;

        public OUApproxVolterra(Grid grid_, double H_, double r_, int N_)
        {
            grid = grid_;
            alpha = H_ - 1 / 2;
            r = r_;
            N = N_;
            //M = M_;
            Func<int,int, double> c_ = (j,i) => Math.Pow(r, ((i - N / 2 - 1) * (1 - alpha))) * (Math.Pow(r, 1 - alpha) - 1) / ((1 - alpha) * AdvancedMath.Gamma(1 - alpha));
            Func<int, int, double> gamma_ = (j, i) => Math.Pow(r, (i - N / 2 - 1)) *
                                            ((1 - alpha) * (Math.Pow(r, 2 - alpha) - 1))
                                            / ((2 - alpha) * (Math.Pow(r, 1 - alpha) - 1));
            C = new ColumnVector(N);
            Gamma = new ColumnVector(N);
            C.Fill(c_);
            Gamma.Fill(gamma_);
        }

        public void simulate(out ColumnVector dw, out double volterra)
        {
            volterra = 0;
            int M = grid.get_timeNmbrStep();
            dw = new ColumnVector(M - 1);
            ColumnVector X = new ColumnVector(M);
            X[0] = 0.0;
            for (int i = 0; i < M-1; i++)
            {
                dw[i] = Math.Sqrt(grid.get_Step()) * simulator.Next();
               // X[i] = ...
            }
        }

        private double OUEuler(double a, double b, double sigma, int N, double dw, double X_0)
        {
            double result = X_0;
            for(int i=1;i<=N;i++)
            {
                result = result + (a - b * result) * result + sigma * dw;
            }
            return result;
        }
    }
}
