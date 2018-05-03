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
    public class HybridScheme
    {
        public int kappa { get; set; }
        public  SquareMatrix choleskyCorrel { get; set; }
        public  Grid grid { get; set; }

        //Gaussian Simulator
        private GaussianSimulator simulator = new GaussianSimulator();
        //private  RandGenNormal normalRng = new RandGenNormal(0, 1);
        private double H;


        private ColumnVector Gamma;

        double var_0;
        double var_1;
        double var_2;
        //Structor
        public HybridScheme(SquareMatrix choleskyCorrel_, Grid grid_, int kappa_ = 2, double H_ = 0.07)
        {
            kappa = kappa_;
            choleskyCorrel = choleskyCorrel_;
            grid = grid_;
            H = H_;
            int n_T = grid.get_timeNmbrStep();
            double n = 1 / grid.get_Step();
            Gamma = new ColumnVector(n_T);
            for (int i = kappa; i < n_T; i++)
            {
                double b_ = ((Math.Pow(i, H + 0.5) - Math.Pow(i - 1, H + 0.5)) / (H + 0.5));
                Gamma[i] = (b_ / Math.Pow(n, H - 0.5));
            }

            var_0 = (double)1 / n;
            var_1 = (double)1 / (Math.Pow(n, 2 * (H - 0.5) + 1) * (2 * (H - 0.5) + 1));
            var_2 = (double)Math.Pow(2, 2 * (H - 0.5) + 1) / (Math.Pow(n, 2 * (H - 0.5) + 1) * (2 * (H - 0.5) + 1));
        }

        public void simulate(out ColumnVector Z, out ColumnVector volterra)
        {
            

            int n_T = grid.get_timeNmbrStep();
            //Z_i defined in A.1
            Z = new ColumnVector(n_T);
            RectangularMatrix Z_k = new RectangularMatrix(2, n_T);

            SimulateCorrelated(Z, Z_k);//Simulate 3 Correlated Gaussians Z_i; Z_k_{1,i}; Z_k_{2,i}  i:= 0, ..., n_T
            //Calculation of Matrix M 
            SquareMatrix M = new SquareMatrix(n_T);
            for (int i = 1; i <= n_T; i++)
            {
                for (int j = 1; j <= Math.Min(kappa, i); j++)
                {
                    M[i - 1, j - 1] = Z_k[j - 1, i - j];
                }
                for (int j = kappa + 1; j <= i; j++)
                {
                    M[i - 1, j - 1] = Z[i - j];
                }
            }

            volterra = new ColumnVector(n_T);
            //Vector of b*
            ColumnVector b = new ColumnVector(n_T);
            for (int i = 1; i <= n_T; i++)
            {
                if (i <= kappa) b[i - 1] = 1;
                else
                {
                    double b_ = ((Math.Pow(i, H + 0.5) - Math.Pow(i - 1, H + 0.5)) / (H + 0.5));
                    b[i - 1] = b_ / Math.Pow(n_T, H - 0.5);
                }
            }
            // 1-simulation of volterra
            volterra = M * b;
        }

        public void simulateFFT(out ColumnVector Z, out ColumnVector volterra)
        {
            int n_T = grid.get_timeNmbrStep();

            //Z_i defined in A.1
            Z = new ColumnVector(n_T);
            RectangularMatrix Z_k = new RectangularMatrix(2,n_T);
            volterra = new ColumnVector(n_T);

            //DateTime starttestM = DateTime.Now;
            //for (int zz = 1; zz <= 1.0E4; zz++)
            //{
            //    SimulateCorrelated(Z, Z_k);
            //}
            //TimeSpan timeExecutiontestM = DateTime.Now - starttestM;

            SimulateCorrelated(Z, Z_k);//Simulate 3 Correlated Gaussians Z_i; Z_k_{1,i}; Z_k_{2,i}  i:= 0, ..., n_T



            //Gamma for Convolution
            double[] Zdouble = Z.ToArray();
            double[] convolution;
            DateTime starttest = DateTime.Now;
            //for (int zz = 1; zz <= 1.0E4; zz++)
            //    alglib.convr1d(Gamma.ToArray(), n_T, Z.ToArray(), n_T, out convolution2);
            //TimeSpan timeExecutiontest = DateTime.Now - starttest;
            alglib.convr1d(Gamma.ToArray(), n_T, Z.ToArray(), n_T, out convolution);

            //DateTime starttestM = DateTime.Now;
            //for (int zz = 1; zz <= 1.0E5; zz++)
            //{

            //}
            //TimeSpan timeExecutiontestM = DateTime.Now - starttestM;

            //var convolutionM = new Double1DConvolution(Z, Gamma.Count());
            //DoubleVector convolution = convolutionM.Convolve(Gamma);

            for (int i = 1; i <= volterra.Count(); i++)
            {
                double v = convolution[i - 1];
                for (int k = 1; k < Math.Min(i, kappa); k++)
                {
                    v += Z_k[k - 1, i - k];
                }
                volterra[i - 1] = v;
            }
        }
        //Simulate the three COrrelated Processus Z ,Z_1 and Z_2
        public  void SimulateCorrelated(ColumnVector Z, RectangularMatrix Z_k)
        {
            double n_T = grid.get_timeNmbrStep();

            for (int i = 0; i < n_T; i++)
            {
                ColumnVector vectorG = new ColumnVector(3);//

                vectorG[0] = simulator.Next();
                vectorG[1] = simulator.Next();
                vectorG[2] = simulator.Next();

                ColumnVector ZVector = choleskyCorrel * vectorG;
                //Adjusting the Variance 
                Z[i] = Math.Sqrt(var_0) * ZVector[0];
                Z_k[0, i] = Math.Sqrt(var_1) * ZVector[1];//Z_{i,1}
                Z_k[1, i] = Math.Sqrt(var_2) * ZVector[2];//Z_{i,2}
            }
        }
    }
}
