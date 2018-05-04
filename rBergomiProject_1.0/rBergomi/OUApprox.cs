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

        //constuct
        public OUApproxVolterra(Grid grid_, double H_, double r_, int N_)
        {
            grid = grid_;
            alpha = H_ - 1 / 2;
            r = r_;
            N = N_;
            //M = M_;
            Func<int,int, double> c_ = (i,j) => Math.Pow(r, ((i+1 - N / 2 - 1) * (1 - alpha))) * (Math.Pow(r, 1 - alpha) - 1) / ((1 - alpha) * AdvancedMath.Gamma(1 - alpha));
            Func<int, int, double> gamma_ = (i,j) => Math.Pow(r, (i+1 - N / 2 - 1)) *
                                            ((1 - alpha) * (Math.Pow(r, 2 - alpha) - 1))
                                            / ((2 - alpha) * (Math.Pow(r, 1 - alpha) - 1));
            C = new ColumnVector(N);
            Gamma = new ColumnVector(N);
            C.Fill(c_);
            Gamma.Fill(gamma_);
        }

        public void simulate(out ColumnVector dw, out double volterra)
        {
            int M = grid.get_timeNmbrStep();
            dw = new ColumnVector(M - 1);
            ColumnVector X = new ColumnVector(N);
            for (int i = 0; i < M-1; i++)
            {
                dw[i] =  simulator.Next();
            }
            for (int i = 0; i < N; i++)
            {
                X[i] = OUEuler(0, Gamma[i], 1, M, dw, 0.0);
            }
            volterra= C.Transpose() * X;

        }
        
        private double OUEuler(double a, double b, double sigma, int M, ColumnVector dw, double X_0)
        {
            double result = X_0;
            for (int i = 0; i < M-1; i++)
            {
                result = result + (a - b * result) * grid.get_Step() + sigma * Math.Sqrt(grid.get_Step()) * dw[i];
            }
            return result;
        }
    }
}
