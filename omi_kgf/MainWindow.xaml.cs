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

        public double[] generateGrid1D(int NumberOfNodes = -1)
        {
            if (NumberOfNodes == -1) NumberOfNodes = this.NumberOfNodes;

            double[] grid = new double[NumberOfNodes];

            for (int i = 0; i < NumberOfNodes; i++)
            {
                grid[i] = node(i);
            }

            return grid;
        }

        public double[] getArrayWithAmendedZeros(double[] shortArray, int desiredLength)
        {
            if (shortArray.Length > desiredLength)
            {
                throw new Exception();
            }
            double[] resultArray = new double[desiredLength];

            for (int i = 0; i < shortArray.Length; i++)
            {
                resultArray[i] = shortArray[i];
            }

            for (int i = shortArray.Length; i < desiredLength; i++)
            {
                resultArray[i] = 0;
            }

            return resultArray;
        }

        //1D Fourier transform
        public double[] getFourierTransform1D(double[] valuesOnGrid, int NumberOfNodes = -1)
        {
            if (NumberOfNodes == -1) NumberOfNodes = this.NumberOfNodes;
            return getFourierTransform1DSlow(valuesOnGrid, NumberOfNodes);
        }

        public double[] getInverseFourierTransform1D(double[] coefficients, int NumberOfNodes = -1)
        {
            if (NumberOfNodes == -1) NumberOfNodes = this.NumberOfNodes;
            return getInverseFourierTransform1DSlow(coefficients, NumberOfNodes);
        }

        //2D Fourier transform
        public double[,] getFourierTransform2D(double[,] valuesOnGrid, int NumberOfNodes = -1)
        {
            if (NumberOfNodes == -1) NumberOfNodes = this.NumberOfNodes;
            return getFourierTransform2DSlow(valuesOnGrid, NumberOfNodes);
        }

        public double[,] getInverseFourierTransform2D(double[,] coefficients, int NumberOfNodes = -1)
        {
            if (NumberOfNodes == -1) NumberOfNodes = this.NumberOfNodes;
            return getInverseFourierTransform2DSlow(coefficients, NumberOfNodes);
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

        private double[] getFourierTransform1DSlow(double[] discreteFunctionValuesOnGrid, int NumberOfNodes)
        {
            if (discreteFunctionValuesOnGrid.Count() != NumberOfNodes)
            {
                throw new Exception();
            }

            var coefficients = new double[NumberOfNodes];

            for (var coefficientIndex = 0; coefficientIndex < NumberOfNodes; coefficientIndex++)
            {
                coefficients[coefficientIndex] = computeFourierCoefficientSlow(coefficientIndex, discreteFunctionValuesOnGrid, NumberOfNodes);
            }

            return coefficients;
        }

        private double[] getInverseFourierTransform1DSlow(double[] coefficients, int NumberOfNodes)
        {
            if (coefficients.Count() != NumberOfNodes)
            {
                throw new Exception();
            }

            var discreteFunctionValuesOnGrid = new double[NumberOfNodes];

            for (var nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++)
            {
                discreteFunctionValuesOnGrid[nodeIndex] = computeDiscreteFunctionValueSlow(nodeIndex, coefficients, NumberOfNodes);
            }

            return discreteFunctionValuesOnGrid;
        }

        private double computeFourierCoefficientSlow(int coefficientIndex, double[] discreteFunctionValuesOnGrid, int NumberOfNodes)
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

        private double computeDiscreteFunctionValueSlow(int nodeIndex, double[] coefficients, int NumberOfNodes)
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

        private double[,] getFourierTransform2DSlow(double[,] discreteFunctionValues, int NumberOfNodes)
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
                    coefficients[coefficientIndex1, coefficientIndex2] = compute2DFourierCoefficientSlow(coefficientIndex1, coefficientIndex2, discreteFunctionValues, NumberOfNodes);
                }
            }

            return coefficients;
        }

        private double[,] getInverseFourierTransform2DSlow(double[,] coefficients, int NumberOfNodes)
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
                    f[functionIndex1, functionIndex2] = compute2DDiscreteFunctionValueSlow(functionIndex1, functionIndex2, coefficients, NumberOfNodes);
                }
            }

            return f;
        }

        private double compute2DFourierCoefficientSlow(int coefficientIndex1, int coefficientIndex2, double[,] discreteFunctionValuesOnGrid, int NumberOfNodes)
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

        private double compute2DDiscreteFunctionValueSlow(int functionIndex1, int functionIndex2, double[,] coefficients, int NumberOfNodes)
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
        Plot2D originalFPlot;
        Plot2D coefficientsOfFPlot;
        Plot2D restoredFPlot;
        Plot2D guessedGPlot;
        Plot2D guessedGCoefficients;
        Plot2D guessedFFromGPlot;

        public MainWindow()
        {
            InitializeComponent();

            // 3D plot
            graphBuilder3dForm = new GraphBuilder2DForm();

            functionPlot3d = new Plot3D("functionPlot3d");
            coefficientsPlot3d = new Plot3D("coefficientsPlot3d");
            restoredFunctionPlot3d = new Plot3D("restoredFunctionPlot3d");

            graphBuilder3dForm.GraphBuilder.Set3DMode();
            graphBuilder3dForm.GraphBuilder.ChartPanel.Aspect.Chart3DPercent = 100;

            graphBuilder3dForm.GraphBuilder.DrawPlot(functionPlot3d);
            graphBuilder3dForm.GraphBuilder.DrawPlot(coefficientsPlot3d);
            graphBuilder3dForm.GraphBuilder.DrawPlot(restoredFunctionPlot3d);
            graphBuilder3dForm.Show();

            // 2D plot
            graphBuilder2dForm = new GraphBuilder2DForm();

            originalFPlot = new Plot2D("original F");
            coefficientsOfFPlot = new Plot2D("coefficients of F");
            restoredFPlot = new Plot2D("restored F");
            guessedGPlot = new Plot2D("guessed G");
            guessedFFromGPlot = new Plot2D("guessed F from G");
            guessedGCoefficients = new Plot2D("guessed coefficients of G");

            graphBuilder2dForm.GraphBuilder.DrawPlot(originalFPlot);
            graphBuilder2dForm.GraphBuilder.DrawPlot(coefficientsOfFPlot);
            graphBuilder2dForm.GraphBuilder.DrawPlot(restoredFPlot);
            graphBuilder2dForm.GraphBuilder.DrawPlot(guessedGPlot);
            graphBuilder2dForm.GraphBuilder.DrawPlot(guessedFFromGPlot);
            graphBuilder2dForm.GraphBuilder.DrawPlot(guessedGCoefficients);
            graphBuilder2dForm.Show();
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            var functionApproximator = new FunctionApproximator();

            functionApproximator.NumberOfNodes = 64;
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

            Func<double, double, double> K = (x, y) => 1.0 / (Abs(Sin(x)) + Abs(Cos(x)) + Abs(Sin(y)) + Abs(Cos(y)));
            Func<double, double> f = x => Sin(x);

            var grid = functionApproximator.generateGrid1D();

            var valuesOfK = functionApproximator.discretizeFunction2D(K);
            var coefficientsOfK = functionApproximator.getFourierTransform2D(valuesOfK);

            functionPlot3d.DiscreteFunction = new DiscreteFunction3D(grid, grid, valuesOfK);
            coefficientsPlot3d.DiscreteFunction = new DiscreteFunction3D(grid, grid, coefficientsOfK);

            restoredFunctionPlot3d.DiscreteFunction = new DiscreteFunction3D(grid, grid, functionApproximator.getInverseFourierTransform2D(coefficientsOfK));

            var valuesOfF = functionApproximator.dicretizeFunction1D(f);
            var coefficientsOfF = functionApproximator.getFourierTransform1D(valuesOfF);

            var coefficientsOfG = functionApproximator.getCoefficientsOfG(valuesOfK, valuesOfF);
            var coefficientsOfGWithZeros = functionApproximator.getArrayWithAmendedZeros(coefficientsOfG, functionApproximator.NumberOfNodes);

            var restoredG = functionApproximator.getInverseFourierTransform1D(coefficientsOfGWithZeros);

            guessedGPlot.DiscreteFunction = new DiscreteFunction2D(grid, restoredG);
            

            originalFPlot.DiscreteFunction = new DiscreteFunction2D(grid, valuesOfF);
            coefficientsOfFPlot.DiscreteFunction = new DiscreteFunction2D(grid, coefficientsOfF);
            restoredFPlot.DiscreteFunction = new DiscreteFunction2D(grid, functionApproximator.getInverseFourierTransform1D(coefficientsOfF));

            var restoredF = new double[functionApproximator.NumberOfNodes];
            for (int i = 0; i < functionApproximator.NumberOfNodes; i++)
            {
                double f_i = 0;
                for (int j = 0; j < functionApproximator.NumberOfNodes; j++)
                {
                    f_i += valuesOfK[i, j] * restoredG[j];
                }
                restoredF[i] = f_i;
            }

            guessedFFromGPlot.DiscreteFunction = new DiscreteFunction2D(grid, restoredF);

            // 3D
            functionPlot3d.Refresh();
            coefficientsPlot3d.Refresh();
            restoredFunctionPlot3d.Refresh();

            // 2D
            originalFPlot.Refresh();
            coefficientsOfFPlot.Refresh();
            restoredFPlot.Refresh();
            guessedGPlot.Refresh();
            //guessedFFromGPlot.Refresh();
            //guessedGCoefficients.Refresh();
        }
    }
}
