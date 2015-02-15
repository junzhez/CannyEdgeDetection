using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Windows.Forms;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace CannyEdgeDetection
{
    class Program
    {
        [STAThread]
        static void Main(string[] args)
        {
            string bitmapFilePath = @"E:\Github\ImageRead\lena.gif";
            Bitmap c = new Bitmap(bitmapFilePath);
            
            int Nx1 = 10;
            double Sigmax1 = 2;
            int Nx2 = 10;
            double Sigmax2 = 2;
            double Theta1 = Math.PI/2;

            int Ny1 = 10;
            double Sigmay1 = 2;
            int Ny2 = 10;
            double Sigmay2 = 2;
            double Theta2 = 0;

            double alpha = 0.1;

            Bitmap d = Helper.Ind2Gray(c);

            var mat = new DenseMatrix(d.Height, d.Width);

            for (int i = 0; i < d.Height; i++)
            {
                for (int j = 0; j < d.Width; j++)
                {
                    mat.At(i, j, d.GetPixel(j, i).R);
                }
            }

            DenseMatrix filterx = Helper.D2dGauss(Nx1, Sigmax1, Nx2, Sigmax2, Theta1);;

            DenseMatrix matx = Helper.Convolute(mat, filterx);

            DenseMatrix filtery = Helper.D2dGauss(Ny1, Sigmay1, Ny2, Sigmay2, Theta2); ;

            DenseMatrix maty = Helper.Convolute(mat, filtery);

            DenseMatrix matx2 = (DenseMatrix) matx.PointwisePower(2);

            DenseMatrix maty2 = (DenseMatrix) maty.PointwisePower(2);

            DenseMatrix matg = matx2 + maty2;

            DenseMatrix R = new DenseMatrix(matg.RowCount, matg.ColumnCount);

            double max = Helper.Max(matg);

            double min = Helper.Min(matg);

            double level = alpha * (max - min) + min;

            for (int i = 1; i < matg.RowCount - 1; i++)
            {
                for (int j = 1; j < matg.ColumnCount - 1; j++)
                {
                    if (matg[i, j] <= level)
                    {
                        R[i, j] = 255;
                    }
                    else
                    {
                        double[] xi =
                        {
                            matx[i, j] / matg[i, j], -matx[i, j] / matg[i, j]
                        };

                        double[] yi = 
                        {
                            maty[i, j] / matg[i, j], -maty[i, j] / matg[i, j]
                        };

                        double[] zi = { 0, 0 };

                        for(int ii = 0; ii < 2; ii++)
                        {
                            if (xi[ii] >= -1 && xi[ii] < 0)
                            {
                                if (yi[ii] >= -1 && yi[ii] < 0)
                                {
                                    double[] X = {-1, 0};
                                    double[] Y = {-1, 0};
                                    double[,] Values = 
                                        {
                                            {matg[i - 1, j - 1], matg[i - 1, j]}, 
                                            {matg[i, j - 1], matg[i, j]}
                                        };

                                    zi[ii] = Helper.Interp2(X, Y, Values, xi[ii], yi[ii]);
                                }
                                else if (yi[ii] >= 0 && yi[ii] <= 1)
                                {
                                    double[] X = { -1, 0 };
                                    double[] Y = { 0, 1 };
                                    double[,] Values = 
                                        {
                                            {matg[i - 1, j], matg[i - 1, j + 1]}, 
                                            {matg[i, j], matg[i, j + 1]}
                                        };

                                    zi[ii] = Helper.Interp2(X, Y, Values, xi[ii], yi[ii]);
                                }
                            }
                            else if (xi[ii] >= 0 && xi[ii] <= 1)
                            {
                                if (yi[ii] >= -1 && yi[ii] < 0)
                                {
                                    double[] X = { 0, 1 };
                                    double[] Y = { -1, 0 };
                                    double[,] Values = 
                                        {
                                            {matg[i, j - 1], matg[i, j]}, 
                                            {matg[i + 1, j - 1], matg[i + 1, j]}
                                        };

                                    zi[ii] = Helper.Interp2(X, Y, Values, xi[ii], yi[ii]);
                                }
                                else if (yi[ii] >= 0 && yi[ii] <= 1)
                                {
                                    double[] X = { 0, 1 };
                                    double[] Y = { 0, 1 };
                                    double[,] Values = 
                                        {
                                            {matg[i, j], matg[i, j + 1]}, 
                                            {matg[i + 1, j], matg[i + 1, j + 1]}
                                        };

                                    zi[ii] = Helper.Interp2(X, Y, Values, xi[ii], yi[ii]);
                                }
                            }

                            if (matg[i, j] >= zi[0] && matg[i, j] >= zi[1])
                            {
                                R[i, j] = 0;
                            }
                            else
                            {
                                R[i, j] = 255;
                            }
                        }
                    }
                }
            }

            Helper.MatrixShow(R);
        }
    }
}
