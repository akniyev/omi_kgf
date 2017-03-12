using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using GraphBuilders;
using DiscreteFunctionsPlots;
using static System.Math;
using DiscreteFunctions;

namespace omi_kgf
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>    

    public partial class MainWindow : Window
    {
        GraphBuilder2DForm graphBuilder;
        Plot2D functionPlot;

        public MainWindow()
        {
            InitializeComponent();
            graphBuilder = new GraphBuilder2DForm();
            functionPlot = new Plot2D();
            graphBuilder.GraphBuilder.DrawPlot(functionPlot);
            graphBuilder.Show();
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            double[] xs = new double[100];
            double[] ys = new double[100];

            for (var i = 0; i < 100; i++)
            {
                xs[i] = 2 * PI / 100.0 * i;
                ys[i] = Sin(xs[i]);
            }

            var discreteFunction = new DiscreteFunction2D(xs, ys);

            functionPlot.DiscreteFunction = discreteFunction;

            functionPlot.Refresh();
        }
    }
}
