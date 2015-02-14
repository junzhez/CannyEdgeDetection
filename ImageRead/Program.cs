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

namespace ImageRead
{
    class Program
    {
        [STAThread]
        static void Main(string[] args)
        {
            string bitmapFilePath = @"J:\Github\ImageRead\lena.gif";
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

            Bitmap d = Program.Ind2Gray(c);

            var mat = new DenseMatrix(d.Height, d.Width);

            for (int i = 0; i < d.Height; i++)
            {
                for (int j = 0; j < d.Width; j++)
                {
                    mat.At(i, j, d.GetPixel(j, i).R);
                }
            }

            DenseMatrix filter = Program.D2dGauss(Nx1, Sigmax1, Nx2, Sigmax2, Theta1);


        }

        static Convolute(DenseMatrix Mat, DenseMatrix Kernel)
        {
            int w = Mat.ColumnCount;
            int h = Mat.RowCount;

            int halfKernelWidth = kernelWidth / 2;
            int halfKernelHeight = kernelHeight / 2;
        }

        static void ImageShow(Bitmap Img)
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

        static Bitmap Ind2Gray(Bitmap Img)
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

        static DenseMatrix D2dGauss(int N1, double Sigma1, int N2, double Sigma2, double Theta)
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

                    h.At(i, j, Program.Gauss(u.At(0, 0), Sigma1) * Program.Gauss(u.At(1, 0), Sigma2));
                }
            }

            return h;
        }

        static double Gauss(double x, double std)
        {
            double y = Math.Exp(-Math.Pow(x, 2) / (2 * Math.Pow(std, 2))) / (std * Math.Sqrt(2 * Math.PI));

            return y;
        }

        static double Dgauss(double x, double std)
        {
            double y = -x * Program.Gauss(x, std) / Math.Pow(std, 2);

            return y;
        }
    }
}
