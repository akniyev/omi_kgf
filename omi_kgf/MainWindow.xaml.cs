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

    class FunctionApproximator
    {
        //Public
        public int NumberOfNodes
        {
            get
            {
                return _numberOfNodes;
            }

            set
            {
                if (!isPowerOfTwo(value))
                {
                    throw new Exception();
                }
                _numberOfNodes = value;
            }
        }

        public int NumberOfTriangleNodes
        {
            get
            {
                return _numberOfTriangleNodes;
            }
            set
            {
                if (value >= NumberOfNodes)
                {
                    throw new Exception();
                }
                _numberOfTriangleNodes = value;
            }
        }

        public double phi(int k, double x)
        {
            if (k == 0)
                return Sqrt(1.0 / NumberOfNodes);
            else
                return Sqrt(2.0 / NumberOfNodes) * Cos(k * x);
        }

        public double node(int i)
        {
            return (2 * i + 1.0) * PI / (2.0 * NumberOfNodes);
        }

        public void setNodesNumberAsPowerOfTwo(int power)
        {
            int result = 1;
            for (var i = 0; i < power; i++)
            {
                result *= 2;
            }
            NumberOfNodes = result;
        }

        //Helper

        public double[] dicretizeFunction1D(Func<double, double> function)
        {
            double[] discreteFunction = new double[NumberOfNodes];

            for (var i = 0; i < NumberOfNodes; i++)
            {
                discreteFunction[i] = function(node(i));
            }

            return discreteFunction;
        }

        public double[,] discretizeFunction2D(Func<double,double,double> function)
        {
            double[,] discreteFunction = new double[NumberOfNodes, NumberOfNodes];

            for (var i = 0; i < NumberOfNodes; i++)
            {
                for (var j = 0; j < NumberOfNodes; j++)
                {
                    discreteFunction[i, j] = function(node(i), node(j));
                }
            }

            return discreteFunction;
        }

        public double[] generateGrid1D()
        {
            double[] grid = new double[NumberOfNodes];

            for (int i = 0; i < NumberOfNodes; i++)
            {
                grid[i] = node(i);
            }

            return grid;
        }

        //1D Fourier transform

        public double[] getFourierTransform1D(double[] valuesOnGrid)
        {
            return getFourierTransform1DSlow(valuesOnGrid);
        }

        public double[] getInverseFourierTransform1D(double[] coefficients)
        {
            return getInverseFourierTransform1DSlow(coefficients);
        }

        //2D Fourier transform

        public double[,] getFourierTransform2D(double[,] valuesOnGrid)
        {
            return getFourierTransform2DSlow(valuesOnGrid);
        }

        public double[,] getInverseFourierTransform2D(double[,] coefficients)
        {
            return getInverseFourierTransform2DSlow(coefficients);
        }

        //Restoring function g's coefficients from II's problem

        public double[] getCoefficientsOfG(double[,] valuesOfK, double[] valuesOfF)
        {
            double[,] coefficientsOfK = getFourierTransform2D(valuesOfK);
            double[] coefficientsOfF = getFourierTransform1D(valuesOfF);

            double[] coefficientsOfG = new double[NumberOfTriangleNodes];

            coefficientsOfG[0] = coefficientsOfF[NumberOfTriangleNodes - 1] / coefficientsOfK[NumberOfTriangleNodes - 1, 0];
            for (var j = 1; j < NumberOfTriangleNodes; j++)
            {
                double sum = 0;
                for (var i = 0; i < j; i++)
                {
                    sum += coefficientsOfK[NumberOfTriangleNodes - 1 - j, i] * coefficientsOfG[i];
                }
                coefficientsOfG[j] = (coefficientsOfF[NumberOfTriangleNodes - 1 - j] - sum) / coefficientsOfK[NumberOfTriangleNodes - 1 - j, j];
            }

            return coefficientsOfG;
        }

        // Private

        private int _numberOfNodes = 16;

        private int _numberOfTriangleNodes = 8;

        private bool isPowerOfTwo(int x)
        {
            x = Abs(x);
            while (x > 1)
            {
                if (x % 2 == 0)
                {
                    x /= 2;
                }
                else
                {
                    return false;
                }
            }
            return true;
        }

        //One-dimensional SLOW Fourier transform

        private double[] getFourierTransform1DSlow(double[] discreteFunctionValuesOnGrid)
        {
            if (discreteFunctionValuesOnGrid.Count() != NumberOfNodes)
            {
                throw new Exception();
            }

            var coefficients = new double[NumberOfNodes];

            for (var coefficientIndex = 0; coefficientIndex < NumberOfNodes; coefficientIndex++)
            {
                coefficients[coefficientIndex] = computeFourierCoefficientSlow(coefficientIndex, discreteFunctionValuesOnGrid);
            }

            return coefficients;
        }

        private double[] getInverseFourierTransform1DSlow(double[] coefficients)
        {
            if (coefficients.Count() != NumberOfNodes)
            {
                throw new Exception();
            }

            var discreteFunctionValuesOnGrid = new double[NumberOfNodes];

            for (var nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++)
            {
                discreteFunctionValuesOnGrid[nodeIndex] = computeDiscreteFunctionValueSlow(nodeIndex, coefficients);
            }

            return discreteFunctionValuesOnGrid;
        }

        private double computeFourierCoefficientSlow(int coefficientIndex, double[] discreteFunctionValuesOnGrid)
        {
            if (coefficientIndex < 0 || coefficientIndex >= NumberOfNodes)
            {
                throw new Exception();
            }

            double coefficient = 0;

            for (var nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++)
            {
                coefficient += discreteFunctionValuesOnGrid[nodeIndex] * phi(coefficientIndex, node(nodeIndex));
            }

            return coefficient;
        }

        private double computeDiscreteFunctionValueSlow(int nodeIndex, double[] coefficients)
        {
            if (coefficients.Count() != NumberOfNodes)
            {
                throw new Exception();
            }

            double functionValue = 0;

            for (var coefficientIndex = 0; coefficientIndex < NumberOfNodes; coefficientIndex++)
            {
                functionValue += coefficients[coefficientIndex] * phi(coefficientIndex, node(nodeIndex));
            }

            return functionValue;
        }

        //Two-dimensional SLOW Fourier transform

        private double[,] getFourierTransform2DSlow(double[,] discreteFunctionValues)
        {
            if (discreteFunctionValues.GetLength(0) != NumberOfNodes || discreteFunctionValues.GetLength(1) != NumberOfNodes)
            {
                throw new Exception();
            }

            double[,] coefficients = new double[NumberOfNodes, NumberOfNodes];

            for (var coefficientIndex1 = 0; coefficientIndex1 < NumberOfNodes; coefficientIndex1++)
            {
                for (var coefficientIndex2 = 0; coefficientIndex2 < NumberOfNodes; coefficientIndex2++)
                {
                    coefficients[coefficientIndex1, coefficientIndex2] = compute2DFourierCoefficientSlow(coefficientIndex1, coefficientIndex2, discreteFunctionValues);
                }
            }

            return coefficients;
        }

        private double[,] getInverseFourierTransform2DSlow(double[,] coefficients)
        {
            if (coefficients.GetLength(0) != NumberOfNodes || coefficients.GetLength(1) != NumberOfNodes)
            {
                throw new Exception();
            }

            double[,] f = new double[NumberOfNodes, NumberOfNodes];

            for (int functionIndex1 = 0; functionIndex1 < NumberOfNodes; functionIndex1++)
            {
                for (int functionIndex2 = 0; functionIndex2 < NumberOfNodes; functionIndex2++)
                {
                    f[functionIndex1, functionIndex2] = compute2DDiscreteFunctionValueSlow(functionIndex1, functionIndex2, coefficients);
                }
            }

            return f;
        }

        private double compute2DFourierCoefficientSlow(int coefficientIndex1, int coefficientIndex2, double[,] discreteFunctionValuesOnGrid)
        {
            if (discreteFunctionValuesOnGrid.GetLength(0) != NumberOfNodes || discreteFunctionValuesOnGrid.GetLength(1) != NumberOfNodes)
            {
                throw new Exception();
            }

            double coefficient = 0;

            for (var functionIndex1 = 0; functionIndex1 < NumberOfNodes; functionIndex1++)
            {
                for (var functionIndex2 = 0; functionIndex2 < NumberOfNodes; functionIndex2++)
                {
                    coefficient += discreteFunctionValuesOnGrid[functionIndex1, functionIndex2] * phi(coefficientIndex1, node(functionIndex1)) * phi(coefficientIndex2, node(functionIndex2));
                }
            }

            return coefficient;
        }

        private double compute2DDiscreteFunctionValueSlow(int functionIndex1, int functionIndex2, double[,] coefficients)
        {
            double functionValue = 0;

            for (var n = 0; n < NumberOfNodes; n++)
            {
                for (var m = 0; m < NumberOfNodes; m++)
                {
                    functionValue += coefficients[n, m] * phi(n, node(functionIndex1)) * phi(m, node(functionIndex2));
                }
            }

            return functionValue;
        }
    }

    public partial class MainWindow : Window
    {
        GraphBuilder2DForm graphBuilder3dForm;
        Plot3D functionPlot3d;
        Plot3D coefficientsPlot3d;
        Plot3D restoredFunctionPlot3d;

        GraphBuilder2DForm graphBuilder2dForm;
        Plot2D plot1;
        Plot2D plot2;
        Plot2D plot3;

        public MainWindow()
        {
            InitializeComponent();

            // 3D plot
            graphBuilder3dForm = new GraphBuilder2DForm();

            functionPlot3d = new Plot3D();
            coefficientsPlot3d = new Plot3D();
            restoredFunctionPlot3d = new Plot3D();

            graphBuilder3dForm.GraphBuilder.Set3DMode();
            graphBuilder3dForm.GraphBuilder.ChartPanel.Aspect.Chart3DPercent = 100;

            graphBuilder3dForm.GraphBuilder.DrawPlot(functionPlot3d);
            graphBuilder3dForm.GraphBuilder.DrawPlot(coefficientsPlot3d);
            graphBuilder3dForm.GraphBuilder.DrawPlot(restoredFunctionPlot3d);
            graphBuilder3dForm.Show();

            // 2D plot
            graphBuilder2dForm = new GraphBuilder2DForm();

            plot1 = new Plot2D();
            plot2 = new Plot2D();
            plot3 = new Plot2D();

            graphBuilder2dForm.GraphBuilder.DrawPlot(plot1);
            graphBuilder2dForm.GraphBuilder.DrawPlot(plot2);
            graphBuilder2dForm.GraphBuilder.DrawPlot(plot3);
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            var functionApproximator = new FunctionApproximator();

            functionApproximator.NumberOfNodes = 32;
            functionApproximator.NumberOfTriangleNodes = 8;

            //double[] xs = functionApproximator.generateGrid1D();
            //double[] ys = functionApproximator.generateGrid1D();
            //double[,] functionValues = functionApproximator.discretizeFunction2D((x, y) => x);
            //double[,] coefficients = functionApproximator.getFourierTransform2D(functionValues);
            //double[,] restoredFunctionValues = functionApproximator.getInverseFourierTransform2D(coefficients);

            //var dfFunction = new DiscreteFunction3D(xs, ys, functionValues);
            //var dfCoefficients = new DiscreteFunction3D(xs, ys, coefficients);
            //var dfRestoredFunction = new DiscreteFunction3D(xs, ys, restoredFunctionValues);

            //functionPlot3d.DiscreteFunction = dfFunction;
            //coefficientsPlot3d.DiscreteFunction = dfCoefficients;
            //restoredFunctionPlot3d.DiscreteFunction = dfRestoredFunction;

            Func<double, double, double> K = (x, y) => 1;
            Func<double, double> f = x => x;

            var valuesOfK = functionApproximator.discretizeFunction2D(K);
            var valuesOfF = functionApproximator.dicretizeFunction1D(f);
            var coefficientsOfG = functionApproximator.getCoefficientsOfG(valuesOfK, valuesOfF);
            
            // 3D
            functionPlot3d.Refresh();
            coefficientsPlot3d.Refresh();
            restoredFunctionPlot3d.Refresh();

            // 2D
            plot1.Refresh();
            plot2.Refresh();
            plot3.Refresh();
        }
    }
}
