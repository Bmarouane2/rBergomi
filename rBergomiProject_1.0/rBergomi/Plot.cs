using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

using ZedGraph.Web;
using ZedGraph;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace rBergomi
{
    public partial class Plot : Form
    {

        public Plot(double[] x_i, double[] y_i,string title,string xName,string yName,string curveName)
        {
            InitializeComponent();

            GraphPane myPane = zedGraphControl1.GraphPane;
            myPane.Title.Text = title;
            myPane.XAxis.Title.Text = xName;
            myPane.YAxis.Title.Text = yName;

            PointPairList VIX_fut_Prices_LogN = new PointPairList();

            rBergomiVIXfuture model = new rBergomiVIXfuture();
            Func<double, double> epsilon_0 = t => Math.Pow(0.235, 2);// *Math.Pow(1 + t, 0.5);

            for (int i = 0; i < x_i.Count(); i++)
            {
                VIX_fut_Prices_LogN.Add(x_i[i],y_i[i]);
                //VIX_fut_Prices_TrCholesky.Add(x_i[i], model.TruncatedCholeskyBergomi(x_i[i], epsilon_0));
            }

            LineItem teamACurve = myPane.AddCurve(curveName,
                 VIX_fut_Prices_LogN, Color.Red, SymbolType.None);
            

            zedGraphControl1.AxisChange();

        }

        public Plot(ColumnVector x_i, ColumnVector y_i, string title, string xName, string yName, string curveName)
        {
            InitializeComponent();

            GraphPane myPane = zedGraphControl1.GraphPane;
            myPane.Title.Text = title;
            myPane.XAxis.Title.Text = xName;
            myPane.YAxis.Title.Text = yName;

            PointPairList VIX_fut_Prices_LogN = new PointPairList();

            rBergomiVIXfuture model = new rBergomiVIXfuture();
            Func<double, double> epsilon_0 = t => Math.Pow(0.235, 2);// *Math.Pow(1 + t, 0.5);

            for (int i = 0; i < x_i.Count(); i++)
            {
                VIX_fut_Prices_LogN.Add(x_i[i], y_i[i]);
                //VIX_fut_Prices_TrCholesky.Add(x_i[i], model.TruncatedCholeskyBergomi(x_i[i], epsilon_0));
            }

            LineItem teamACurve = myPane.AddCurve(curveName,
                 VIX_fut_Prices_LogN, Color.Red, SymbolType.None);


            zedGraphControl1.AxisChange();

        }

        public Plot(ColumnVector x_i, ColumnVector y_i, ColumnVector y_i2, string title, string xName, string yName, string curveName)
        {
            InitializeComponent();

            GraphPane myPane = zedGraphControl1.GraphPane;
            myPane.Title.Text = title;
            myPane.XAxis.Title.Text = xName;
            myPane.YAxis.Title.Text = yName;

            PointPairList VIX_fut_Prices_LogN = new PointPairList();
            PointPairList VIX_fut_Prices_LogN2 = new PointPairList();

            for (int i = 0; i < x_i.Count(); i++)
            {
                VIX_fut_Prices_LogN.Add(x_i[i], y_i[i]);
                VIX_fut_Prices_LogN2.Add(x_i[i], y_i2[i]);
                //VIX_fut_Prices_TrCholesky.Add(x_i[i], model.TruncatedCholeskyBergomi(x_i[i], epsilon_0));
            }

            LineItem teamACurve = myPane.AddCurve(curveName,
                 VIX_fut_Prices_LogN, Color.Red, SymbolType.None);
            LineItem teamACurve2 = myPane.AddCurve(curveName,
     VIX_fut_Prices_LogN2, Color.Blue, SymbolType.None);

            zedGraphControl1.AxisChange();

        }

        public Plot(ColumnVector x_i, ColumnVector y_i, ColumnVector y_i2, ColumnVector y_i3)
        {
            InitializeComponent();

            GraphPane myPane = zedGraphControl1.GraphPane;
            myPane.Title.Text = "VIX Futures";
            myPane.XAxis.Title.Text = "T";
            myPane.YAxis.Title.Text = "VIX futures";

            PointPairList VIX_fut_Prices_LogN = new PointPairList();
            PointPairList VIX_fut_Prices_LogN2 = new PointPairList();
            PointPairList VIX_fut_Prices_LogN3 = new PointPairList();

            for (int i = 0; i < x_i.Count(); i++)
            {
                VIX_fut_Prices_LogN.Add(x_i[i], y_i[i]);
                VIX_fut_Prices_LogN2.Add(x_i[i], y_i2[i]);
                VIX_fut_Prices_LogN3.Add(x_i[i], y_i3[i]);
                //VIX_fut_Prices_TrCholesky.Add(x_i[i], model.TruncatedCholeskyBergomi(x_i[i], epsilon_0));
            }

            LineItem teamACurve = myPane.AddCurve("Log Normal",
                 VIX_fut_Prices_LogN, Color.Red, SymbolType.None);

            LineItem teamACurve2 = myPane.AddCurve("Hybrid + Euler",
     VIX_fut_Prices_LogN2, Color.Blue, SymbolType.None);

            LineItem teamACurve3 = myPane.AddCurve("Truncated Cholesky",
VIX_fut_Prices_LogN3, Color.Green, SymbolType.None);

            zedGraphControl1.AxisChange();

        }

        private void zedGraphControl1_Load(object sender, EventArgs e)
        {

        }
    }
}
