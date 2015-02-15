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
    public class Helper
    {
        public static double Interp2(
            double[] X, double[] Y,
            double[,] Values,
            double x, double y)
        {
            double q11 = Values[0, 0] * (X[1] - x) * (Y[1] - y) / ((X[1] - X[0]) * (Y[1] - Y[0]));
            double q21 = Values[1, 0] * (x - X[0]) * (Y[1] - y) / ((X[1] - X[0]) * (Y[1] - Y[0]));
            double q12 = Values[0, 1] * (X[1] - x) * (y - Y[0]) / ((X[1] - X[0]) * (Y[1] - Y[0]));
            double q22 = Values[1, 1] * (x - X[0]) * (y - Y[0]) / ((X[1] - X[0]) * (Y[1] - Y[0]));

            return q11 + q21 + q12 + q22;
        }

        public static double Max(DenseMatrix Mat)
        {
            double result = Double.MinValue;

            for (int i = 0; i < Mat.RowCount; i++)
            {
                for (int j = 0; j < Mat.ColumnCount; j++)
                {
                    if (Mat.At(i, j) > result)
                    {
                        result = Mat.At(i, j);
                    }
                }
            }

            return result;
        }

        public static double Min(DenseMatrix Mat)
        {
            double result = Double.MaxValue;

            for (int i = 0; i < Mat.RowCount; i++)
            {
                for (int j = 0; j < Mat.ColumnCount; j++)
                {
                    if (Mat.At(i, j) < result)
                    {
                        result = Mat.At(i, j);
                    }
                }
            }

            return result;
        }

        public static void MatrixShow(DenseMatrix Mat)
        {
            Bitmap img = new Bitmap(Mat.ColumnCount, Mat.RowCount);

            for (int i = 0; i < img.Width; i++)
            {
                for (int j = 0; j < img.Height; j++)
                {
                    int val = (int)Mat.At(j, i);
                    img.SetPixel(i, j, Color.FromArgb(val, val, val));
                }
            }

            Helper.ImageShow(img);
        }

        public static DenseMatrix Convolute(DenseMatrix Mat, DenseMatrix Kernel)
        {
            int w = Mat.ColumnCount;
            int h = Mat.RowCount;

            var result = new DenseMatrix(h, w);

            int halfKernelWidth = Kernel.ColumnCount / 2;
            int halfKernelHeight = Kernel.RowCount / 2;

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    double sum = 0;
                    for (int x = -halfKernelHeight; x < halfKernelHeight; x++)
                    {
                        for (int y = -halfKernelWidth; y < halfKernelWidth; y++)
                        {
                            int ki = x + halfKernelHeight;
                            int kj = y + halfKernelWidth;

                            int hi = i - x;
                            int hj = j - y;

                            double hval = 0;

                            if ((hi >= 0 && hi < h)
                                && (hj >= 0 && hj < w))
                            {
                                hval = Mat.At(hi, hj);
                            }

                            sum += hval * Kernel.At(ki, kj);
                        }
                    }

                    result.At(i, j, sum);
                }
            }

            return result;
        }

        public static void ImageShow(Bitmap Img)
        {
            PictureBox P = new PictureBox
            {
                Size = Img.Size,
            };

            Form form = new Form
            {
                Size = Img.Size,
            };

            P.Image = Img;

            form.Controls.Add(P);
            form.Show();

            Application.Run(form);
        }

        public static Bitmap Ind2Gray(Bitmap Img)
        {
            int height = Img.Height;
            int width = Img.Width;

            Bitmap d = new Bitmap(width, height);

            // Loop through the images pixels to reset color.
            for (int x = 0; x < Img.Width; x++)
            {
                for (int y = 0; y < Img.Height; y++)
                {
                    Color pixelColor = Img.GetPixel(x, y);
                    int grayScale = (int)((pixelColor.R * 0.3) + (pixelColor.G * 0.59) + (pixelColor.B * 0.11));
                    Color newColor = Color.FromArgb(pixelColor.A, grayScale, grayScale, grayScale);
                    d.SetPixel(x, y, newColor); // Now greyscale
                }
            }

            return d;
        }

        public static DenseMatrix D2dGauss(int N1, double Sigma1, int N2, double Sigma2, double Theta)
        {
            var r = new DenseMatrix(2, 2);

            r.At(0, 0, Math.Cos(Theta));
            r.At(0, 1, -Math.Sin(Theta));
            r.At(1, 0, Math.Sin(Theta));
            r.At(1, 1, Math.Cos(Theta));

            var h = new DenseMatrix(N2, N1);

            for (int i = 0; i < N2; i++)
            {
                for (int j = 0; j < N1; j++)
                {
                    var tmp = new DenseMatrix(2, 1);
                    tmp.At(0, 0, j - (N1 + 1) / 2);
                    tmp.At(1, 0, i - (N2 + 1) / 2);
                    var u = r.Multiply(tmp);

                    h.At(i, j, Helper.Gauss(u.At(0, 0), Sigma1) * Helper.Gauss(u.At(1, 0), Sigma2));
                }
            }

            return h;
        }

        public static double Gauss(double x, double std)
        {
            double y = Math.Exp(-Math.Pow(x, 2) / (2 * Math.Pow(std, 2))) / (std * Math.Sqrt(2 * Math.PI));

            return y;
        }

        public static double Dgauss(double x, double std)
        {
            double y = -x * Helper.Gauss(x, std) / Math.Pow(std, 2);

            return y;
        }
    }
}
