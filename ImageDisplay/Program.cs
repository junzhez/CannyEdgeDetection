using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Drawing;
using System.Windows.Forms;

namespace ImageDisplay
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);

            string bitmapFilePath = @"J:\Github\ImageRead\lena.gif";
            Bitmap c = new Bitmap(bitmapFilePath);
            int hight = c.Height;
            int width = c.Width;

            Bitmap d = new Bitmap(width, hight);

            // Loop through the images pixels to reset color.
            for (int x = 0; x < c.Width; x++)
            {
                for (int y = 0; y < c.Height; y++)
                {
                    Color pixelColor = c.GetPixel(x, y);
                    int grayScale = (int)((pixelColor.R * 0.3) + (pixelColor.G * 0.59) + (pixelColor.B * 0.11));
                    Color newColor = Color.FromArgb(pixelColor.A, grayScale, grayScale, grayScale);
                    d.SetPixel(x, y, newColor); // Now greyscale
                }
            }

            PictureBox P = new PictureBox();

            P.Image = c;
        }
    }
}
