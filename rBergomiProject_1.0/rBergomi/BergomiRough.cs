using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

using Meta.Numerics;
using Meta.Numerics.Matrices;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using System.IO;




namespace rBergomi
{
    public class rBergomiVIXfuture
    {
        //Model Inputs
        double H = 0.07;
        double Ch = 0.2893;
        double vi = 1.2287;// 1.9*Ch*sqrt(2 * H) / 2;
        double Delta = 1 / 12.0;
        double n_tH;

        public double VIXfuture_LogNormal(double T, Func<double,double> epsilon)
        {
            //For Integrale calculation
            EvaluationSettings settings = new EvaluationSettings();
            settings.AbsolutePrecision = 1.0E-5;
            settings.RelativePrecision = 0.0;

	        //Grid On [O,T]
	        Grid grid=new Grid(0, T, 1000);

            // Variance Calcul
            Func<double,double> f= s=>Math.Pow(Math.Pow(T-s+Delta,H+0.5)-Math.Pow(T-s,H+0.5),2);
            IntegrationResult integr_f = FunctionMath.Integrate(f, Meta.Numerics.Interval.FromEndpoints(0.0, T), settings);
            double var = (4 * Math.Pow(vi * Ch, 2)) / (Math.Pow(Delta * (H + 0.5), 2)) * integr_f.Value;
            

            IntegrationResult integr_epsilon = FunctionMath.Integrate(epsilon, Meta.Numerics.Interval.FromEndpoints(T, T+Delta), settings);
            //LogNormal Solution APproximation
            return Math.Pow(Delta, -0.5) * Math.Sqrt(integr_epsilon.Value) * Math.Exp(-var / 8);
        }

        public ColumnVector VIXfuture_HybridMethod(double T, Func<double,double> epsilon)
        {
            //Result
            ColumnVector VIXFutures;
            //For Exécution Time
            DateTime start = DateTime.Now;
            //Input:
            int kappa = 2;
           
            // the Gri t_0 ..... t_{100} = T
            int n = 500;
            Grid grid = new Grid(0, T, (int)Math.Abs(T * n));
            
            int n_T = grid.get_timeNmbrStep();
            n_tH = Math.Pow(n_T, H - 0.5); //optimizeMC

            //The second Grid  t_0=T, t_1 ..... t_N = T + Delta
            int Nbstep = 20;

            //Correlation :
            SymmetricMatrix correl = MakeCorrel(n, grid);
            SquareMatrix choleskyCorrel = correl.CholeskyDecomposition().SquareRootMatrix();

            int period = 100;
            VIXFutures = new ColumnVector(grid.get_timeNmbrStep() / period);

            
            //--------------------------------------------------------------------------------------------------------------------------------------------------------------
            //-----------------------------------MC----------------------------------------------------------------------------------------------------------------------
            //-----------------------------------------------------------------------------------------------------------------------------------------------------------
            var MCnbofSimulation = 3.0E4;
            HybridScheme hybridscheme = new HybridScheme(choleskyCorrel, grid, kappa, H);
            for (int mc = 1; mc <= MCnbofSimulation; mc++)
            {
                // 1-simulation of volterra Hybrid Scheme
                ColumnVector Z;
                ColumnVector volterra;
                ColumnVector Z_Brownian;

                hybridscheme.simulateFFT(out Z, out volterra);

                //2 - extract the path of the Brownian motion Z driving the Volterra process
                Z_Brownian = ExtractBrownian(kappa, n_T, Z, volterra);

                //Grid secondGrid2 = new Grid(grid.t(n_T), grid.t(n_T) + Delta, 20);
                //ColumnVector volterraT2 = EulerMEthode(grid, secondGrid2, n_T, volterra, Z_Brownian);
                //double VIX_future_i2 = VIX_Calculate(epsilon, secondGrid2, volterraT2);
                //VIXFutures[0] += VIX_future_i2;


                //DateTime starttestM = DateTime.Now;
                //for (int zz = 1; zz <= 1.0E4; zz++)
                //{
                //    hybridscheme.simulateFFT(out Z, out volterra);
                //}
                //TimeSpan timeExecutiontestM = DateTime.Now - starttestM;

                //DateTime starttest = DateTime.Now;
                //for (int zz = 1; zz <= 1.0E4; zz++)
                //{
                //    for (int i = 1; i <= volterra.Count() / period; i++)
                //    {
                //        int N_i2 = i * period;

                //        //3- approximate the continuous-time process V^T by the discrete-time version V^T~ defined via the  forward Euler scheme
                //        double[] volterraT2 = EulerMEthode(grid, secondGrid, N_i2, volterra, Z_Brownian);
                //        //double VIX_future_i2 = VIX_Calculate(epsilon, secondGrid, volterraT2);
                //        //VIXFutures[i - 1] += VIX_future_i2;
                //    }
                //}
                //TimeSpan timeExecutiontest = DateTime.Now - starttest;

                for (int i = 1; i <= volterra.Count() / period; i++)
                {
                    int N_i = i * period;
                    Grid secondGrid = new Grid(grid.t(N_i), grid.t(N_i) + Delta, Nbstep);
                    //3- approximate the continuous-time process V^T by the discrete-time version V^T~ defined via the  forward Euler scheme
                    double[] volterraT = EulerMEthode(grid, secondGrid, N_i, volterra, Z_Brownian);
                    double VIX_future_i = VIX_Calculate(epsilon, secondGrid, volterraT);
                    VIXFutures[i - 1] += VIX_future_i;
                }
            }// ENd of MC
            //MC Mean
            for (int i = 0; i < VIXFutures.Count(); i++)
            {
                VIXFutures[i] /= MCnbofSimulation;
            }
            TimeSpan timeExecution = DateTime.Now - start;
            return VIXFutures;
        }

        private ColumnVector EulerMEthode2(Grid grid, ColumnVector volterra, ColumnVector Z_Brownian)//Not Used !! 
        {
            int Nbstep = 20;
            Func<double, double> epsilon = t => Math.Pow(0.235, 2);

            var VIXFutures = new ColumnVector(grid.get_timeNmbrStep() / 50);
           

            for (int k = 1; k <= volterra.Count() / 50; k++)
            {
                int N_k = k * 50;
                Grid secondGrid = new Grid(grid.t(N_k), grid.t(N_k) + Delta, Nbstep);
                //3- approximate the continuous-time process V^T by the discrete-time version V^T~ defined via the  forward Euler scheme
                int N_tau = secondGrid.get_timeNmbrStep();
                // Volterra  V^T
                double[] volterraT = new double[N_tau + 1];

                volterraT[0] = volterra[N_k - 1];

                for (int j = 1; j <= N_tau; j++)
                {
                    double v = 0.0;
                    double[] v2= new double[volterra.Count() / 50];
                    for (int i = 1; i <= N_k; i++)
                    {
                        v += (Z_Brownian[i] - Z_Brownian[i - 1]) / Math.Pow((secondGrid.t(j) - grid.t(i - 1)), -H + 0.5);
                    }
                    volterraT[j] = (v);
                }

                double VIX_future_i = VIX_Calculate(epsilon, secondGrid, volterraT);
                VIXFutures[k - 1] += VIX_future_i;
            }
            return VIXFutures;
        }

        public double VIXfuture_TruncatedChlsky(double T, Func<double, double> epsilon,out double StddeviationMC)
        {
            DateTime start = DateTime.Now;
            //StreamWriter sw = new StreamWriter("C:\\Users\\Marouane\\Desktop\\M2IF\\rough volatility\\rBergomiFile.txt");
            StddeviationMC = 0.0;
            //Grid
            int N = 100;
            int S = 7;
            Grid grid = new Grid(T, T + Delta, N);

            #region  Correlation "Correl Matrix Construction"
            // Small Covariance Matrix t_0 , ... , t_7
            SymmetricMatrix covM_Small = new SymmetricMatrix(S);
            //the full Covariance Matrix
            SymmetricMatrix covM = new SymmetricMatrix(N+1);
            //Setting for integrale approximations
            EvaluationSettings settings = new EvaluationSettings();
            settings.AbsolutePrecision = 1.0E-7;
            settings.RelativePrecision = 0.0;

            // Small Covariance Calcul t_0,.....,t_7
            for (int i = 0; i < S; i++)
            {
                covM_Small[i, i] = (Math.Pow(grid.t(i), 2 * H) - Math.Pow(grid.t(i) -T, 2 * H) )/ (2 * H);
                for (int j = 0; j < i; j++)
                {
                    Func<double, double> covfunc = t => Math.Pow((grid.t(j) - t) * (grid.t(i) - t), H - 0.5);
                    IntegrationResult integresult = FunctionMath.Integrate(covfunc, Meta.Numerics.Interval.FromEndpoints(0.0, T), settings);
                    covM_Small[i, j] = integresult.Value;
                }
            }
            // full COrrelation Calcul
            for (int i = 0; i <= N; i++)
            {
                covM[i, i] = (Math.Pow(grid.t(i), 2 * H) - Math.Pow(grid.t(i) - T, 2 * H)) / (2 * H);
                for (int j = 0; j < i; j++)
                {
                    Func<double, double> covfunc = t => Math.Pow((grid.t(j) - t) * (grid.t(i) - t), H - 0.5);
                    IntegrationResult integresult = FunctionMath.Integrate(covfunc, Meta.Numerics.Interval.FromEndpoints(0.0, T), settings);
                    covM[i, j] = integresult.Value;
                }
            }

            Func<int, int, double> corrf_Small = (i, j) => covM_Small[i, j] / (Math.Sqrt(covM_Small[i, i] * covM_Small[j, j]));
            SymmetricMatrix Correl_Small = new SymmetricMatrix(S);
            Correl_Small.Fill(corrf_Small);

            Func<int, int, double> corrf = (i, j) => covM[i, j] / (Math.Sqrt(covM[i, i] * covM[j, j]));
            SymmetricMatrix Correl = new SymmetricMatrix(N+1);
            Correl.Fill(corrf);
            #endregion

            CholeskyDecomposition cholesky = Correl_Small.CholeskyDecomposition();
            SquareMatrix choleskyCorrel_Small = new SquareMatrix(S);
            choleskyCorrel_Small = cholesky.SquareRootMatrix();

            GaussianSimulator simulator = new GaussianSimulator();
            double VIX = 0.0;
            var MC = 1.0E5;
            for (int mc = 1; mc < MC; mc++)
            {
                ColumnVector GaussianVector = new ColumnVector(S);
                // Simulating Volterra at 8 first steps on [T, T+Delta]
                for (int i = 0; i < S; i++)
                {
                    GaussianVector[i] = simulator.Next();
                }
                ColumnVector Volterra_small = choleskyCorrel_Small * GaussianVector;
                //Adjusting the variance of Volterra Processus
                for (int i = 0; i < S; i++)
                {
                    Volterra_small[i] = Volterra_small[i] * Math.Sqrt((Math.Pow(grid.t(i), 2 * H) - Math.Pow(grid.t(i) - grid.t(0), 2 * H)) / (2 * H));
                }

                // Simulating VOlterra on t_8 ... t_N with the truncated Formula
                double[] Volterra = new double[N + 1];
                for (int i = 0; i <= N; i++)
                {
                    if (i < S) Volterra[i] = Volterra_small[i];
                    //Tranceted Formula
                    else
                        Volterra[i] = Math.Sqrt(covM[i, i]) * (Correl[i, i - 1] * Volterra[i - 1] / Math.Sqrt(covM[i - 1, i - 1]) + Math.Sqrt(1 - Math.Pow(Correl[i, i - 1], 2)) * simulator.Next());
                }
    
                double VIX_ = VIX_Calculate(epsilon, grid, Volterra);
                VIX += VIX_;
                StddeviationMC += Math.Pow(VIX_, 2);
            }
            //sw.Close();
            VIX /= MC;
            StddeviationMC = Math.Sqrt(StddeviationMC / MC - Math.Pow(VIX, 2))/Math.Sqrt(MC);
            TimeSpan dur = DateTime.Now - start;
            return VIX;
        }

        public double VIXfuture_TruncatedChlsky(double T, Func<double, double> epsilon)
        {
            DateTime start = DateTime.Now;
            //StreamWriter sw = new StreamWriter("C:\\Users\\Marouane\\Desktop\\M2IF\\rough volatility\\rBergomiFile.txt");
            //Grid
            int N = 100;
            int S = 7;
            Grid grid = new Grid(T, T + Delta, N);

            #region  Correlation "Correl Matrix Construction"
            // Small Covariance Matrix t_0 , ... , t_7
            SymmetricMatrix covM_Small = new SymmetricMatrix(S);
            //the full Covariance Matrix
            SymmetricMatrix covM = new SymmetricMatrix(N + 1);
            //Setting for integrale approximations
            EvaluationSettings settings = new EvaluationSettings();
            settings.AbsolutePrecision = 1.0E-7;
            settings.RelativePrecision = 0.0;

            // Small Covariance Calcul t_0,.....,t_7
            for (int i = 0; i < S; i++)
            {
                covM_Small[i, i] = (Math.Pow(grid.t(i), 2 * H) - Math.Pow(grid.t(i) - T, 2 * H)) / (2 * H);
                for (int j = 0; j < i; j++)
                {
                    Func<double, double> covfunc = t => Math.Pow((grid.t(j) - t) * (grid.t(i) - t), H - 0.5);
                    IntegrationResult integresult = FunctionMath.Integrate(covfunc, Meta.Numerics.Interval.FromEndpoints(0.0, T), settings);
                    covM_Small[i, j] = integresult.Value;
                }
            }
            // full COrrelation Calcul
            for (int i = 0; i <= N; i++)
            {
                covM[i, i] = (Math.Pow(grid.t(i), 2 * H) - Math.Pow(grid.t(i) - T, 2 * H)) / (2 * H);
                for (int j = 0; j < i; j++)
                {
                    Func<double, double> covfunc = t => Math.Pow((grid.t(j) - t) * (grid.t(i) - t), H - 0.5);
                    IntegrationResult integresult = FunctionMath.Integrate(covfunc, Meta.Numerics.Interval.FromEndpoints(0.0, T), settings);
                    covM[i, j] = integresult.Value;
                }
            }

            Func<int, int, double> corrf_Small = (i, j) => covM_Small[i, j] / (Math.Sqrt(covM_Small[i, i] * covM_Small[j, j]));
            SymmetricMatrix Correl_Small = new SymmetricMatrix(S);
            Correl_Small.Fill(corrf_Small);

            Func<int, int, double> corrf = (i, j) => covM[i, j] / (Math.Sqrt(covM[i, i] * covM[j, j]));
            SymmetricMatrix Correl = new SymmetricMatrix(N + 1);
            Correl.Fill(corrf);
            #endregion

            CholeskyDecomposition cholesky = Correl_Small.CholeskyDecomposition();
            SquareMatrix choleskyCorrel_Small = new SquareMatrix(S);
            choleskyCorrel_Small = cholesky.SquareRootMatrix();

            GaussianSimulator simulator = new GaussianSimulator();
            double VIX = 0.0;
            var MC = 1.0E5;
            for (int mc = 1; mc < MC; mc++)
            {
                ColumnVector GaussianVector = new ColumnVector(S);
                // Simulating Volterra at 8 first steps on [T, T+Delta]
                for (int i = 0; i < S; i++)
                {
                    GaussianVector[i] = simulator.Next();
                }
                ColumnVector Volterra_small = choleskyCorrel_Small * GaussianVector;
                //Adjusting the variance of Volterra Processus
                for (int i = 0; i < S; i++)
                {
                    Volterra_small[i] = Volterra_small[i] * Math.Sqrt((Math.Pow(grid.t(i), 2 * H) - Math.Pow(grid.t(i) - grid.t(0), 2 * H)) / (2 * H));
                }

                // Simulating VOlterra on t_8 ... t_N with the truncated Formula
                double[] Volterra = new double[N + 1];
                for (int i = 0; i <= N; i++)
                {
                    if (i < S) Volterra[i] = Volterra_small[i];
                    //Tranceted Formula
                    else
                        Volterra[i] = Math.Sqrt(covM[i, i]) * (Correl[i, i - 1] * Volterra[i - 1] / Math.Sqrt(covM[i - 1, i - 1]) + Math.Sqrt(1 - Math.Pow(Correl[i, i - 1], 2)) * simulator.Next());
                }

                double VIX_ = VIX_Calculate(epsilon, grid, Volterra);
                VIX += VIX_;
            }
            //sw.Close();
            VIX /= MC;
            TimeSpan dur = DateTime.Now - start;
            return VIX;
        }

        public double PriceVIXOption(VIXOption Option,Func<double, double> epsilon,double stock_0)
        {
            int n = 500;
            double T = Option.maturity;
            Grid grid = new Grid(0, T, (int)Math.Abs(T * n));
            int n_T =grid.get_timeNmbrStep();

            int kappa = 2;
            //Correlation :
            SymmetricMatrix correl = MakeCorrel(n, grid);
            SquareMatrix choleskyCorrel = correl.CholeskyDecomposition().SquareRootMatrix();

            double rho = 0.1;//correlation between the two Brownian
            Func<double, double> payoff;
            switch (Option.type)
            {
                case VIXOption.OptionType.Call:
                    payoff = S => Math.Max(S - Option.strike, 0);
                    break;
                case VIXOption.OptionType.Put:
                    payoff = S => Math.Max(Option.strike - S, 0);
                    break;
                default:
                    payoff = S => Math.Max(S - Option.strike, 0);
                    break;
            }

            double price = 0.0;
            double McNbSimulation = 1E5;
            for (int mc = 1; mc <= McNbSimulation;mc ++ )
            {
                // 1-simulation of volterra Hybrid Scheme
                ColumnVector Z;
                ColumnVector volterra;
                //ColumnVector Z_Brownian;

                HybridScheme hybridscheme = new HybridScheme(choleskyCorrel, grid, kappa, H);
                hybridscheme.simulate(out Z, out volterra);

                GaussianSimulator simulator = new GaussianSimulator();
                //Z_Brownian = ExtractBrownian(kappa, n_T, Z, volterra);

                ColumnVector Variance = new ColumnVector(n_T + 1);
                Variance[0] = epsilon(grid.t(0));
                for (int i = 1; i <= grid.get_timeNmbrStep(); i++)
                {
                    Variance[i] = epsilon(grid.t(i)) * Math.Exp(2 * vi * Ch * volterra[i - 1] - Math.Pow(vi * Ch, 2) * Math.Pow(grid.t(i), 2 * H));
                }
                double X = Math.Log(stock_0);
                for (int i = 0; i < n_T; i++)
                {
                    double dW = (rho * Z[i] + Math.Sqrt(1 - Math.Pow(rho, 2)) * Math.Sqrt(grid.get_Step()) * simulator.Next());
                    X = X - 0.5 * Variance[i] * grid.get_Step() + Math.Sqrt(Variance[i]) * dW;
                }

                double S = Math.Exp(X);
                price += payoff(S);
            }
            price /= McNbSimulation;
            return price;
        }
        //----------------------------------------------------------------------------------------------------------------------------------------------------------
        //Private Methodes:----------------------------------------------------------------------------------------------------------------------------------------------------------

        //Extract the path of the Brownian motion Z driving the Volterra process
        private ColumnVector ExtractBrownian(int kappa, int n_T, ColumnVector Z, ColumnVector volterra)
        {
            ColumnVector Z_Brownian = new ColumnVector(n_T +1);
            Z_Brownian[0] = 0.0;
            for (int i = 1; i <= n_T; i++)
            {
                if (i == 1) Z_Brownian[i] = (Z_Brownian[i - 1] + n_tH * (volterra[i - 1]));
                else
                {
                    if (i <= kappa) Z_Brownian[i] = (Z_Brownian[i - 1] + n_tH * (volterra[i - 1] - volterra[i - 2]));
                    else Z_Brownian[i] = (Z_Brownian[i - 1] + (Z[i - 2]));
                }
            }
            return Z_Brownian;
        }
        // approximate the continuous-time process V^T by the discrete-time version V^T~ defined via the  forward Euler scheme
        private double[] EulerMEthode(Grid grid, Grid secondGrid, int N_i, ColumnVector volterra, ColumnVector Z_Brownian)
        {
            int N_tau = secondGrid.get_timeNmbrStep();
            // Volterra  V^T
            double[] volterraT = new double[N_tau + 1];

            volterraT[0] = volterra[N_i - 1];

            for (int j = 1; j <= N_tau; j++)
            {
                double v = 0.0;
                for (int i = 1; i <= N_i; i++)
                {
                    v += (Z_Brownian[i]-Z_Brownian[i - 1]) / Math.Pow((secondGrid.t(j) - grid.t(i - 1)), -H + 0.5);
                }
                volterraT[j] = (v);
            }
            return volterraT;
        }
        // VIX future Calcul with Trapèze Approximation
        private double VIX_Calculate(Func<double, double> epsilon, Grid grid, double[] volterraT)
        {
            int N_tau = grid.get_timeNmbrStep();
            ColumnVector Q = new ColumnVector(N_tau + 1);
            for (int j = 0; j <= N_tau; j++)
            {
                Q[j] = (epsilon(grid.t(j)) * Math.Exp(2 * vi * Ch * volterraT[j]) * Math.Exp(Math.Pow(vi * Ch, 2) / H * (Math.Pow(grid.t(j) - grid.t(0), 2 * H) - Math.Pow(grid.t(j), 2 * H))));
            }

            double VIX_future_i = 0.0;
            for (int j = 0; j < N_tau; j++)
            {
                VIX_future_i += (Q[j] + Q[j + 1]) / 2 * (grid.get_Step());
            }

            VIX_future_i /= Delta;
            VIX_future_i = Math.Sqrt(VIX_future_i);
            return VIX_future_i;
        }
        // COrrel Making for Hybrid Method!
        public SymmetricMatrix MakeCorrel(int n, Grid grid)
        {
            //Correlation between Z , Z_1 and Z_2
            //integrationCalcul
            Func<double, double> func = t => (double)Math.Pow((1 * grid.get_Step() - t), (H - 0.5)) * Math.Pow((2 * grid.get_Step() - t), (H - 0.5));
            EvaluationSettings settings = new EvaluationSettings();
            settings.AbsolutePrecision = 1.0E-9;
            settings.RelativePrecision = 0.0;
            IntegrationResult integresult = FunctionMath.Integrate(func, Meta.Numerics.Interval.FromEndpoints(0.0, grid.get_Step()), settings);

            #region Correlation Calcul (correl Matrix)
            SymmetricMatrix correl = new SymmetricMatrix(3);
            for (int i = 0; i < 3; i++)
            {
                correl[i, i] = 1.0;
            }
            double var_0 = (double)1 / n;
            double var_1 = (double)1 / (Math.Pow(n, 2 * (H - 0.5) + 1) * (2 * (H - 0.5) + 1));
            double var_2 = (double)Math.Pow(2, 2 * (H - 0.5) + 1) / (Math.Pow(n, 2 * (H - 0.5) + 1) * (2 * (H - 0.5) + 1));
            correl[1, 0] = (double)1 / (Math.Pow(n, H + 0.5) * (H + 0.5))/ (Math.Pow(var_0 * var_1, 0.5));
            correl[2, 0] = (double)(Math.Pow(2, H + 0.5) - 1) / (Math.Pow(n, H + 0.5) * (H + 0.5)) / (Math.Pow(var_0 * var_2, 0.5));
            correl[2, 1] = integresult.Value/ (Math.Pow(var_2 * var_1, 0.5));
            #endregion
            return correl;
        }

    }
}
