using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.IO;

using rBergomi;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;


using CenterSpace.NMath.Core;

namespace rBergomiTestProject
{
    [TestClass]
    public class BergomiTest
    {
        [TestMethod]
        public void TestCorrelatedSimulation()//for Hybrid Methode
        {
            StreamWriter Z_text = new StreamWriter("C:\\Users\\Marouane\\Desktop\\M2IF\\rough volatility\\FileZ.txt");
            StreamWriter Z1_text = new StreamWriter("C:\\Users\\Marouane\\Desktop\\M2IF\\rough volatility\\FileZ1.txt");
            StreamWriter Z2_text = new StreamWriter("C:\\Users\\Marouane\\Desktop\\M2IF\\rough volatility\\FileZ2.txt");

            rBergomiVIXfuture model = new rBergomiVIXfuture();
            int n = 500;
            double T = 0.5;
            Grid grid = new Grid(0, T, (int)Math.Abs(T * n));

            //Correl
            SymmetricMatrix correl = model.MakeCorrel(n,grid);

            SquareMatrix choleskyCorrel = correl.CholeskyDecomposition().SquareRootMatrix();

            //Gaussian Simulator
            var simulator = new GaussianSimulator();
            double mc = 1.0E5;

            ColumnVector Z_test = new ColumnVector((int)mc);
            ColumnVector Z_test_1 = new ColumnVector((int)mc);
            ColumnVector Z_test_2 = new ColumnVector((int)mc);
            ColumnVector G = new ColumnVector((int)mc);

            for (int i = 1; i <= mc; i++)
            {
                DoubleVector Z = new DoubleVector(grid.get_timeNmbrStep());
                RectangularMatrix Z_k = new RectangularMatrix(2, grid.get_timeNmbrStep());

                HybridScheme hybridscheme = new HybridScheme(choleskyCorrel, grid, 2, 0.07);

                hybridscheme.SimulateCorrelated(Z, Z_k);//Simulate 3 Correlated Gaussians Z_i; Z_k_{1,i}; Z_k_{2,i}  i:= 0, ..., n_T
                int mid = (int)grid.get_timeNmbrStep() / 2;

                G[i - 1] = simulator.Next();

                Z_test[i - 1] = Z[mid];
                Z_test_1[i - 1] = Z_k[0, mid];
                Z_test_2[i - 1] = Z_k[1, mid];

                string Ztext_ = "";
                string Ztext1_ = "";
                string Ztext2_ = "";
                for (int k = 0; k < Z.RowCount; k++)
                {
                    Ztext_ += (Z[k].ToString() + "; ");
                    Ztext1_ += (Z_k[0, k].ToString() + "; ");
                    Ztext2_ += (Z_k[1, k].ToString() + "; ");
                }
                Z_text.WriteLine(Ztext_);
                Z1_text.WriteLine(Ztext1_);
                Z2_text.WriteLine(Ztext2_);
            }

            Z_text.Close();
            Z1_text.Close();
            Z2_text.Close();

            BivariateSample mybivariate = new BivariateSample();
            mybivariate.Add(Z_test, Z_test_2);

            double correlcoef = mybivariate.CorrelationCoefficient;
            double correlreal = correl[2, 0];

            Sample mysample = new Sample(Z_test);
            double mean = mysample.Mean;//must be =0
            double variance = mysample.Variance;//must be =1/n
            var result = mysample.KolmogorovSmirnovTest(new NormalDistribution());

            Sample mysampleG = new Sample(G);
            var result2 = mysampleG.KolmogorovSmirnovTest(new NormalDistribution());

        }
    }
}
