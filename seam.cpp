#include <cassert>
#include <iostream>
#include <limits>
#include <tgmath.h>
#include <vector>

#include "seam.h"

using namespace std;
constexpr double INF(numeric_limits<double>::max());

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
    return ((int(red * 255)) << 16) + ((int(green * 255)) << 8) + int(blue * 255);
}

// Returns the RGB components from given grayscale value (between 0.0 and 1.0).
int get_RGB(double gray)
{
    int grayint(gray * 255);
    return (grayint << 16) + (grayint << 8) + grayint;
}

GrayImage to_gray(const RGBImage &cimage)
{
    GrayImage gray;
    vector<double> line;
    size_t max_i(cimage.size());
    size_t max_j(cimage[0].size());
    for (size_t i(0); i < max_i; ++i)
    {
        line.clear();
        for (size_t j(0); j < max_j; ++j)
        {
            line.push_back(get_gray(cimage[i][j]));
        }
        gray.push_back(line);
    }
    return gray;
}

RGBImage to_RGB(const GrayImage &gimage)
{
    RGBImage color;
    vector<int> line;
    size_t max_i(gimage.size());
    size_t max_j(gimage[0].size());
    for (size_t i(0); i < max_i; ++i)
    {
        line.clear();
        for (size_t j(0); j < max_j; ++j)
        {
            line.push_back(get_RGB(gimage[i][j]));
        }
        color.push_back(line);
    }
    return color; 
}

// ***********************************
// TASK 2: FILTER
// ***********************************


void clamp(long &val, long max)
{
    val = (val < 0 ? 0 : (val > max ? max : val));
}

//Calcule la valeur de convolution entre un kernel et une portion d'image.
// Le pixel traité de coordonnees row, col prendra cette valeur.
double convol(const GrayImage &gray, const Kernel &kernel, const int row, const int col)
{
    double val(0);
    int nb_rows(gray.size()-1);
    int nb_cols(gray[0].size()-1);
    int kernel_w(kernel.size());
    int kernel_h(kernel[0].size());
    for (size_t i(0); i < kernel_w; ++i)
    {
        for (size_t j(0); j < kernel_h; ++j)
        {
            long clamp_i(i-kernel_w/2+row);
            long clamp_j(j-kernel_h/2+col);
            clamp(clamp_i, nb_rows);
            clamp(clamp_j, nb_cols);
            val = val + kernel[i][i] * gray[clamp_i][clamp_j];
        }
    }
    return val;
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
            // GrayImage smallimage = subimage(gray, j, i, kernel.size(), kernel[0].size());
            double val = convol(gray, kernel, i, j);
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

//Calcule le numero de la node.
inline size_t get_id(size_t row, size_t col, size_t width)
{
    return (row * width + col);
}

//Créé un vector contenant les successors possibles de chaque pixel.
vector<size_t> get_successors(size_t i, size_t j, size_t width, size_t height)
{
    vector<size_t> successors;

    if (i == (height - 1))
    {
        successors = {width * height + 1};
    }
    else if (j == 0)
    {
        successors = {
            get_id(i + 1, j, width),
            get_id(i + 1, j + 1, width)};
    }
    else if (j == (width - 1))
    {
        successors = {
            get_id(i + 1, j - 1, width),
            get_id(i + 1, j, width)};
    }
    else
    {
        successors = {
            get_id(i + 1, j - 1, width),
            get_id(i + 1, j, width),
            get_id(i + 1, j + 1, width)};
    }
    return successors;
}

Graph create_graph(const GrayImage &gray)
{
    Graph graph_im;
    size_t width(gray[0].size());
    for (size_t i(0); i < gray.size(); ++i)
    {
        for (size_t j(0); j < width; ++j)
        {
            vector<size_t> successors;
            successors = get_successors(i, j, width, gray.size());
            double cost = gray[i][j];
            Node nodeij = {successors, cost, INF, 0};
            graph_im.push_back(nodeij);
        }
    }

    vector<size_t> successors;
    for (size_t j(0); j < gray[0].size(); ++j)
    {
        successors.push_back(j);
    }
    Node firstnode = {successors, 0, INF, 0};
    Node lastnode = {{}, 0, INF, 0};
    graph_im.push_back(firstnode);
    graph_im.push_back(lastnode);
    return graph_im;
}

Path shortest_path(Graph &graph, size_t from, size_t to)
{
    graph[from].distance_to_target = graph[from].costs;

    bool modified(true);
    while (modified)
    {
        modified = false;
        for (size_t i(0); i < graph.size(); ++i)
        {
            for (size_t k(0); k < graph[i].successors.size(); ++k)
            {
                size_t n(graph[i].successors[k]);
                double new_distance(graph[i].distance_to_target + graph[n].costs);
                if (graph[n].distance_to_target > new_distance)
                {
                    graph[n].distance_to_target = new_distance;
                    graph[n].predecessor_to_target = i;
                    modified = true;
                }
            }
        }
    }

    Path inverse_path;
    size_t node_id = graph[to].predecessor_to_target;
    while (node_id != from)
    {
        inverse_path.push_back(node_id);
        node_id = graph[node_id].predecessor_to_target;
    }
    Path path;
    for (size_t i(inverse_path.size() - 1); i >= 0 and i < inverse_path.size(); --i)
    {
        path.push_back(inverse_path[i]);
    }
    return path;
};

//Calcule la colonne d'une node étant donné son numero.
size_t get_col(size_t number, size_t width)
{
    size_t col = number % width;
    return col;
}

Path find_seam(const GrayImage &gray)
{
    Graph graph_image = create_graph(gray);
    size_t from(gray.size() * gray[0].size());
    size_t to(gray.size() * gray[0].size() + 1);
    Path s_path = shortest_path(graph_image, from, to);
    Path x_coordinates(s_path.size());
    for (size_t i(0); i < s_path.size(); ++i)
    {
        x_coordinates[i] = get_col(s_path[i], gray[0].size());
    }
    return x_coordinates;
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
