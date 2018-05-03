using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using System.Data.OleDb;

using Meta.Numerics.Matrices;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics;


namespace rBergomi
{

    class Program
    {
        static void Main(string[] args)
        {

            //Ploting results Graph
            //Application.EnableVisualStyles();
            //Application.SetCompatibleTextRenderingDefault(false);
            //Application.Run(new Graph());


            rBergomiVIXfuture model = new rBergomiVIXfuture();
            Func<double, double> epsilon_0 = t => Math.Pow(0.235, 2);// *Math.Sqrt(1 + t);
            double T = 2.0;
            //// the Gri t_0 ..... t_{100} = T
            int n = 500;
            Grid grid = new Grid(0, T, (int)Math.Abs(T * n));

            //lognormal
            double vixfuture_LogNormal = model.VIXfuture_LogNormal(grid.t(35), epsilon_0);
           // //Mc with hybrid Scheme
            ColumnVector vixfuture_Hybrid = model.VIXfuture_HybridMethod(T, epsilon_0);
            ColumnVector newvixfutures = new ColumnVector(vixfuture_Hybrid.Count());
            ColumnVector newvixlognorm = new ColumnVector(vixfuture_Hybrid.Count());
            ColumnVector newvixtruncated = new ColumnVector(vixfuture_Hybrid.Count());

            ColumnVector ti = new ColumnVector(vixfuture_Hybrid.Count());
            int perdiod = 100;
            for (int i = 0; i < vixfuture_Hybrid.Count();i++)
            {
                int N_i = (i + 1) * perdiod;
                ti[i] = grid.t((N_i));
                newvixfutures[i] = vixfuture_Hybrid[i];
                newvixtruncated[i] = model.VIXfuture_TruncatedChlsky(ti[i], epsilon_0);
                newvixlognorm[i] = model.VIXfuture_LogNormal(ti[i], epsilon_0);
            }
            Application.Run(new Plot(ti, newvixlognorm, newvixfutures,newvixtruncated));
            //Application.Run(new Plot(ti, newvixfutures,newvixlognorm, "Volterra Simulation in hybrid Method", "T + Delta ", "V", "volterra process"));
            //Application.Run(new Plot(ti, newvixlognorm, "Volterra Simulation in hybrid Method", "T + Delta ", "V", "volterra process"));
            
           // //Mc with trancated Cholesky
            double TrancatedMethode_MCstdDeviation;
            double vixfuture_trancated = model.VIXfuture_TruncatedChlsky(T, epsilon_0,out  TrancatedMethode_MCstdDeviation);
           
            //S&P Option Pricing
            VIXOption option = new VIXOption(VIXOption.OptionType.Call, 1.0, 100);
            double optionPrice = model.PriceVIXOption(option, epsilon_0, 100);
        }//Main

    }//programm
}//namespace
