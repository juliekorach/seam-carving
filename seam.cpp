#include <cassert>
#include <iostream>
#include <limits>
#include <tgmath.h>
#include <vector>

#include "seam.h"

using namespace std;

// ***********************************
// TASK 1: COLOR
// ***********************************

// Returns red component (in the scale 0.0-1.0) from given RGB color.
// 110100001 10101111 10001110 01111010

double get_red(int rgb)
{
    return ((rgb >> 16) & 0xff) / 255.0;
}

// Returns green component (in the scale 0.0-1.0) from given RGB color.
double get_green(int rgb)
{
    int green((rgb >> 8) & 0xff);
    return green / 255.0;
}

// Returns blue component (in the scale 0.0-1.0) from given RGB color.
double get_blue(int rgb)
{
    return (rgb & 0xff) / 255.0;
}

// Returns the average of red, green and blue components from given RGB color. (Scale: 0.0-1.0)
double get_gray(int rgb)
{
    return (get_blue(rgb) + get_red(rgb) + get_green(rgb)) / 3;
}

// Returns the RGB value of the given red, green and blue components.
int get_RGB(double red, double green, double blue)
{
    int rgb;
    rgb = ((int(red * 255)) << 16) + ((int(green * 255)) << 8) + int(blue * 255);
    return rgb;
}

// Returns the RGB components from given grayscale value (between 0.0 and 1.0).
int get_RGB(double gray)
{
    int grayint(gray * 255);
    grayint = (grayint << 16) + (grayint << 8) + grayint;
    return grayint;
}

// Converts  RGB image to grayscale double image.
GrayImage to_gray(const RGBImage &cimage)
{
    GrayImage gray;
    vector<double> line;
    for (size_t i(0); i < cimage.size(); ++i)
    {
        line.clear();
        for (size_t j(0); j < cimage[i].size(); ++j)
        {
            line.push_back(get_gray(cimage[i][j]));
        }
        gray.push_back(line);
    }

    return {gray}; // TODO MODIFY AND COMPLETE
}

// Converts grayscale double image to an RGB image.
RGBImage to_RGB(const GrayImage &gimage)
{
    RGBImage color;
    vector<int> line;
    for (size_t i(0); i < gimage.size(); ++i)
    {
        line.clear();
        for (size_t j(0); j < gimage[i].size(); ++j)
        {
            line.push_back(get_RGB(gimage[i][j]));
        }
        color.push_back(line);
    }

    return {color}; // TODO MODIFY AND COMPLETE
}

// ***********************************
// TASK 2: FILTER
// ***********************************

// Get a pixel without accessing out of bounds
// return nearest valid pixel color
void clamp(long &val, long max)
{
    if (val < 0)
    {
        val = 0;
    }
    else if (val > max)
    {
        val = max;
    }
}

double convol(const GrayImage &gray, const Kernel &kernel)
{
    double val(0);
    for (size_t i(0); i < kernel.size(); ++i)
    {
        for (size_t j(0); j < kernel[0].size(); ++j)
        {
            val = val + kernel[i][j] * gray[i][j];
        }
    }
    return val;
}

GrayImage subimage(const GrayImage &gray, long x, long y, int a, int b)
{
    GrayImage img;
    vector<double> line;
    long max_i(gray.size() - 1);
    long max_j(gray[0].size() - 1);
    long shiftx((a - 1) / 2);
    long shifty((b - 1) / 2);
    for (long i(y - shifty); i <= y + shifty; ++i)
    {
        line.clear();
        for (long j(x - shiftx); j <= x + shiftx; ++j)
        {
            long clamp_i(i);
            long clamp_j(j);
            clamp(clamp_i, max_i);
            clamp(clamp_j, max_j);
            line.push_back(gray[clamp_i][clamp_j]);
        }
        img.push_back(line);
    }
    return img;
}
// Convolve a single-channel image with the given kernel.
GrayImage filter(const GrayImage &gray, const Kernel &kernel)
{
    GrayImage filtered;
    vector<double> line;
    long max_i(gray.size());
    long max_j(gray[0].size());

    for (long i(0); i < max_i; ++i)
    {
        line.clear();
        for (long j(0); j < max_j; ++j)
        {
            GrayImage smallimage = subimage(gray, j, i, kernel.size(), kernel[0].size());
            double val = convol(smallimage, kernel);

            line.push_back(val);
        }
        filtered.push_back(line);
    }
    return filtered;
}

// Smooth a single-channel image
GrayImage smooth(const GrayImage &gray)
{
    vector<vector<double>> kernelsmooth = {{0.1, 0.1, 0.1},
                                           {0.1, 0.2, 0.1},
                                           {0.1, 0.1, 0.1}};
    return filter(gray, kernelsmooth);
}

// Compute horizontal Sobel filter

GrayImage sobelX(const GrayImage &gray)
{
    vector<vector<double>> kernelsobelx = {{-1.0, 0.0, 1.0},
                                           {-2.0, 0.0, 2.0},
                                           {-1.0, 0.0, 1.0}};
    return filter(gray, kernelsobelx);
}

// Compute vertical Sobel filter

GrayImage sobelY(const GrayImage &gray)
{
    vector<vector<double>> kernelsobely = {{-1.0, -2.0, -1.0},
                                           {0.0, 0.0, 0.0},
                                           {1.0, 2.0, 1.0}};
    return filter(gray, kernelsobely);
}

// Compute the magnitude of combined Sobel filters

GrayImage sobel(const GrayImage &gray)
{
    GrayImage imgx = sobelX(gray);
    GrayImage imgy = sobelY(gray);
    vector<double> empty_line(gray[0].size(), 0);
    GrayImage mag(gray.size(), empty_line);
    for (size_t i(0); i < gray.size(); ++i)
    {
        for (size_t j(0); j < gray[0].size(); ++j)
        {
            mag[i][j] = sqrt((imgx[i][j] * imgx[i][j]) + (imgy[i][j] * imgy[i][j]));
        }
    }
    return mag;
}

// ************************************
// TASK 3: SEAM
// ************************************

Graph create_graph(const GrayImage &gray)
{
    return {}; // TODO MODIFY AND COMPLETE
}

// Return shortest path from Node from to Node to
// The path does NOT include the from and to Node
Path shortest_path(Graph &graph, size_t from, size_t to)
{

    return {}; // TODO MODIFY AND COMPLETE
};

Path find_seam(const GrayImage &gray)
{
    return {}; // TODO MODIFY AND COMPLETE
}

// ***********************************
// TASK 3 provided functions
// Highlight or remove seam from RGB or gray image
// ***********************************

// Draw a seam on a gray image
// return a new gray image with the seam in black
GrayImage highlight_seam(const GrayImage &gray, const Path &seam)
{
    GrayImage result(gray);
    // Paint seam in black
    for (size_t row(0); row < seam.size(); ++row)
    {
        result[row][seam[row]] = 0;
    }
    return result;
}

// Draw a seam on an RGB image
// return a new RGB image with the seam in blue
RGBImage highlight_seam(const RGBImage &image, const Path &seam)
{
    RGBImage result(image);
    // Paint seam in blue
    for (size_t row(0); row < seam.size(); ++row)
    {
        result[row][seam[row]] = 0x000ff;
    }
    return result;
}

// Remove specified seam from a gray-scale image
// return the new gray image (width is decreased by 1)

GrayImage remove_seam(const GrayImage &gray, const Path &seam)
{
    GrayImage result(gray);
    for (size_t row(0); row < seam.size(); ++row)
    {
        result[row].erase(result[row].begin() + seam[row]);
    }
    return result;
}

// Remove specified seam from an RGB image
// return the new RGB image (width is decreased by 1)
RGBImage remove_seam(const RGBImage &image, const Path &seam)
{
    RGBImage result(image);
    for (size_t row(0); row < seam.size(); ++row)
    {
        result[row].erase(result[row].begin() + seam[row]);
    }
    return result;
}
