#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "matrix.h"
#include "io.h"
#include <initializer_list>
#include <limits>

using std::tuple;
using std::get;
using std::tie;
using std::make_tuple;
using std::string;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;
using std::numeric_limits;

#include "io.h"
#include "matrix.h"

void print_help(const char *argv0)
{
    const char *usage =
R"(where PARAMS are from list:

--align [--gray-world | --unsharp | --autocontrast [<fraction>]]
    align images with one of postprocessing functions

--gaussian <sigma> [<radius>=1]
    gaussian blur of image, 0.1 < sigma < 100, radius = 1, 2, ...

--gaussian-separable <sigma> [<radius>=1]
    same, but gaussian is separable

--sobel-x
    Sobel x derivative of image

--sobel-y
    Sobel y derivative of image

--unsharp
    sharpen image

--gray-world
    gray world color balancing

--autocontrast [<fraction>=0.0]
    autocontrast image. <fraction> of pixels must be croped for robustness

--resize <scale>
    resize image with factor scale. scale is real number > 0

--canny <threshold1> <threshold2>
    apply Canny filter to grayscale image. threshold1 < threshold2,
    both are in 0..360

--custom <kernel_string>
    convolve image with custom kernel, which is given by kernel_string, example:
    kernel_string = '1,2,3;4,5,6;7,8,9' defines kernel of size 3

[<param>=default_val] means that parameter is optional.
)";
    cout << "Usage: " << argv0 << " <input_image_path> <output_image_path> "
         << "PARAMS" << endl;
    cout << usage;
}

/*approximately equal function*/
bool aeq(double x, double y, double eps = 1e-5)
{
    if (abs(x - y) < eps)  return true;
    return false;
}

/*convert radian to neighbourhoud place*/
int radtoplace(double rad)
{
    if (aeq(rad,0)) return 1;
    if (aeq(rad,M_PI/4)) return 2;
    if (aeq(rad,M_PI/2)) return 3;
    if (aeq(rad,M_PI*3/4)) return 4;
    if (aeq(rad,M_PI) || (aeq(rad,-M_PI))) return 5;
    if (aeq(rad,-M_PI*3/4)) return 6;
    if (aeq(rad,-M_PI/2)) return 7;
    if (aeq(rad,-M_PI/4)) return 8;
    return 0;
}

template<typename ValueType>
ValueType read_value(string s)
{
    stringstream ss(s);
    ValueType res;
    ss >> res;
    if (ss.fail() or not ss.eof())
        throw string("bad argument: ") + s;
    return res;
}

template<typename ValueT>
void check_number(string val_name, ValueT val, ValueT from,
                  ValueT to=numeric_limits<ValueT>::max())
{
    if (val < from)
        throw val_name + string(" is too small");
    if (val > to)
        throw val_name + string(" is too big");
}

void check_argc(int argc, int from, int to=numeric_limits<int>::max())
{
    if (argc < from)
        throw string("too few arguments for operation");

    if (argc > to)
        throw string("too many arguments for operation");
}

Matrix<int> parse_kernel(string kernel)
{
    // Kernel parsing implementation here
    return Matrix<int>(0, 0);
}

/*class Gauss
{
public:
    Gauss(double s = 1.4, int r = 1) : sigma(s), radius(r) {}
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        double gauss[size][size];
        uint r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        double div = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                gauss[i][j] = (1/(M_PI*2)*sigma)*exp((i*i+j*j)/(2*sigma*sigma)*(-1.0));
                div += gauss[i][j];
            }
        }
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                sum_r += r * gauss[i][j]/div;
                sum_g += g * gauss[i][j]/div;
                sum_b += b * gauss[i][j]/div;
            }
        }
        return make_tuple(sum_r, sum_g, sum_b);
    }
    // Radius of neighbourhoud, which is passed to that operator
    double sigma;
    int radius;
};*/

template<typename T>
class Custom
{
public:
    Custom(Matrix<T> &kernel, bool saveimg = true) : ker(kernel), radius((ker.n_cols - 1) / 2), save(saveimg) {}
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        /*Matrix<int> ker = {{-1, 0, 1},
                          {-2, 0, 2},
                          {-1, 0, 1}};
                          {{1, 1, 1},
                          {1, 1, 1},
                          {1, 1, 1}};*/
        uint size = 2 * radius + 1;
        int r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        T elem = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                elem = ker(i,j);
                sum_r += r * elem;
                sum_g += g * elem;
                sum_b += b * elem;
            }
        }
        if (save) {
            if (sum_r < 0) sum_r = 0;
            if (sum_r > 255) sum_r = 255;
            if (sum_g < 0) sum_g = 0;
            if (sum_g > 255) sum_g = 255;
            if (sum_b < 0) sum_b = 0;
            if (sum_b > 255) sum_b = 255;
        }
        return make_tuple(sum_r, sum_g, sum_b);
    }
    // Radius of neighbourhoud, which is passed to that operator
    Matrix<T> ker;
    int radius;
    bool save;
};

class Separ_x
{
public:
    Separ_x(Matrix<double> &kernel) : ker(kernel), radius((ker.n_cols - 1) / 2) {}
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        int r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        double elem = 0;
        for (uint i = 0; i < size; ++i) {
            // Tie is useful for taking elements from tuple
            tie(r, g, b) = m(radius, i);
            elem = ker(radius,i);
            sum_r += r * elem;
            sum_g += g * elem;
            sum_b += b * elem;
        }
        return make_tuple(sum_r, sum_g, sum_b);
    }
    // Radius of neighbourhoud, which is passed to that operator
    Matrix<double> ker;
    int radius;
};

class Separ_y
{
public:
    Separ_y(Matrix<double> &kernel) : ker(kernel), radius((ker.n_cols - 1) / 2) {}
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        int r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        double elem = 0;
        for (uint i = 0; i < size; ++i) {
            // Tie is useful for taking elements from tuple
            tie(r, g, b) = m(i, radius);
            elem = ker(i,radius);
            sum_r += r * elem;
            sum_g += g * elem;
            sum_b += b * elem;
        }
        return make_tuple(sum_r, sum_g, sum_b);
    }
    // Radius of neighbourhoud, which is passed to that operator
    Matrix<double> ker;
    int radius;
};

Image gaussian(const Image &src_image, double sigma = 1.4, uint radius = 1)
{
    uint size = 2 * radius + 1;
    Matrix<double> kernel(size, size);
    double div = 0;
    for (uint i = 0; i < size; ++i) {
        for (uint j = 0; j < size; ++j) {
            kernel(i,j) = 1/((M_PI*2)*sigma*sigma)*exp((i*i+j*j)/(2*sigma*sigma)*(-1.0));
            div += kernel(i,j);
        }
    }
    for (uint i = 0; i < size; ++i) {
        for (uint j = 0; j < size; ++j) {
            kernel(i,j) /= div;
        }
    }
    return src_image.unary_map(Custom<double>(kernel));
}

Image gaussian_separable(const Image &src_image, double sigma = 1.4, uint radius = 1)
{
    uint size = 2 * radius + 1;
    Matrix<double> kernel_x(size, size), kernel_y(size, size);
    double div = 0;
    for (uint j = 0; j < size; ++j) {
        kernel_y(j,radius) = 1/(sqrt(M_PI*2)*sigma)*exp((j*j)/(2*sigma*sigma)*(-1.0));
        kernel_x(radius,j) = 1/(sqrt(M_PI*2)*sigma)*exp((j*j)/(2*sigma*sigma)*(-1.0));
        div += kernel_y(j,radius);
    }
    for (uint i = 0; i < size; ++i) {
        kernel_y(i,radius) /= div;
        kernel_x(radius,i) /= div;
    }
    //cout<<gauss.n_cols<<endl;
    return src_image.unary_map(Separ_x(kernel_x)).unary_map(Separ_y(kernel_y));
    //return src_image.unary_map(Separ_x(kernel_x));
}

Image sobel_x(const Image &src_image, bool save = true)
{
    Matrix<int> kernel = {{-1, 0, 1},
                            {-2, 0, 2},
                            {-1, 0, 1}};
    return src_image.unary_map(Custom<int>(kernel, save));
}

Image sobel_y(const Image &src_image, bool save = true)
{
    Matrix<int> kernel = {{1, 2, 1},
                            {0, 0, 0},
                            {-1, -2, -1}};
    return src_image.unary_map(Custom<int>(kernel, save));
}

Image gray_world(const Image &src_image, bool save = true)
{
    int r, g, b;
    double s = 0, sr = 0, sg = 0, sb = 0;
    uint sizei = src_image.n_rows, sizej = src_image.n_cols;
    Image img(src_image.n_rows, src_image.n_cols);
    for (uint i = 0; i < sizei; ++i) {
        for (uint j = 0; j < sizej; ++j) {
            // Tie is useful for taking elements from tuple
            tie(r, g, b) = src_image(i, j);
            sr += r;
            sg += g;
            sb += b;
        }
    }
    sr /= sizei * sizej;
    sg /= sizei * sizej;
    sb /= sizei * sizej;
    s = (sr + sg + sb) / 3;
    for (uint i = 0; i < sizei; ++i) {
        for (uint j = 0; j < sizej; ++j) {
            tie(r, g, b) = src_image(i, j);
            r *= s / sr;
            g *= s / sg;
            b *= s / sb;
            if (save) {
                //if (sum_r < 0) sum_r = 0;
                if (r > 255) r = 255;
                //if (sum_g < 0) sum_g = 0;
                if (g > 255) g = 255;
                //if (sum_b < 0) sum_b = 0;
                if (b > 255) b = 255;
            }
            img(i,j) = make_tuple(r, g, b);
        }
    }
    return img;
}

Image unsharp(const Image &src_image, bool save = true)
{
    Matrix<double> kernel = {{-1/6, -2/3, -1/6},
                            {-2/3, 13/3, -2/3},
                            {-1/6, -2/3, -1/6}};
    return src_image.unary_map(Custom<double>(kernel, save));
}

/*Image aligng(const Image &src_image, string postprocessing, double fraction = 0, bool save = true)
{
    switch (postprocessing) {
        case "--gray-world":
            return gray_world(src_image);
        case "--unsharp"
            return unsharp(src_image);
        case "--autocontrast"
            return ;
        default: return src_image;
    }
}*/

Image canny(const Image &m, int th1, int th2)
{
    Matrix<double> sobelx = {{-1, 0, 1},
                            {-2, 0, 2},
                            {-1, 0, 1}};
    Matrix<double> sobely = {{1, 2, 1},
                            {0, 0, 0},
                            {-1, -2, -1}};
    Image imgx = m.unary_map(Custom<double>(sobelx)),
        imgy = m.unary_map(Custom<double>(sobely));
    Matrix<int> grad(m.n_rows, m.n_cols), bordmap(m.n_rows, m.n_cols);
    Matrix<double> theta(m.n_rows, m.n_cols);
    unsigned int rx, ry;

    for (uint i = 0; i < m.n_rows; ++i) {
        for (uint j = 0; j < m.n_cols; ++j) {
            tie(rx,rx,rx) = imgx(i,j); tie(ry,ry,ry) = imgy(i,j);
            grad(i,j) = sqrt(rx*rx+ry*ry);
            theta(i,j) = atan2(rx, ry);
        }
    }
    for (uint i = 0; i < m.n_rows; ++i) {
        for (uint j = 0; j < m.n_cols; ++j) {
            switch (radtoplace(theta(i,j))) {
                case 1:
                    if ((grad(i,j) <=  grad(i,j + 1)) || (grad(i,j) <= grad(i,j - 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0;
                    if (grad(i,j) > th2)
                        bordmap(i,j) = 1;
                case 2:
                    if ((grad(i,j) <=  grad(i + 1,j + 1)) || (grad(i,j) <= grad(i - 1,j - 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0;
                    if (grad(i,j) > th2)
                        bordmap(i,j) = 1;
                case 3:
                    if ((grad(i,j) <=  grad(i + 1,j)) || (grad(i,j) <= grad(i - 1,j)) || (grad(i,j) < th1))
                        grad(i,j) = 0;
                    if (grad(i,j) > th2)
                        bordmap(i,j) = 1;
                case 4:
                    if ((grad(i,j) <=  grad(i + 1,j + 1)) || (grad(i,j) <= grad(i - 1,j - 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0;
                    if (grad(i,j) > th2)
                        bordmap(i,j) = 1;
                case 5:
                    if ((grad(i,j) <=  grad(i,j - 1)) || (grad(i,j) <= grad(i,j + 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0;
                    if (grad(i,j) > th2)
                        bordmap(i,j) = 1;
                case 6:
                    if ((grad(i,j) <=  grad(i - 1,j - 1)) || (grad(i,j) <= grad(i + 1,j + 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0;
                    if (grad(i,j) > th2)
                        bordmap(i,j) = 1;
                case 7:
                    if ((grad(i,j) <=  grad(i - 1,j)) || (grad(i,j) <= grad(i + 1,j)) || (grad(i,j) < th1))
                        grad(i,j) = 0;
                    if (grad(i,j) > th2)
                        bordmap(i,j) = 1;
                case 8:
                    if ((grad(i,j) <=  grad(i - 1,j + 1)) || (grad(i,j) <= grad(i + 1,j - 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0;
                    if (grad(i,j) > th2)
                        bordmap(i,j) = 1;
                default: ;
            }
        }
    }


    return m;
}

int main(int argc, char **argv)
{
    try {
        check_argc(argc, 2);
        if (string(argv[1]) == "--help") {
            print_help(argv[0]);
            return 0;
        }

        check_argc(argc, 4);
        Image src_image = load_image(argv[1]), dst_image;

        string action(argv[3]);

        if (action == "--sobel-x") {
            check_argc(argc, 4, 4);
            dst_image = sobel_x(src_image);
        } else if (action == "--sobel-y") {
            check_argc(argc, 4, 4);
            dst_image = sobel_y(src_image);
        } else if (action == "--unsharp") {
            check_argc(argc, 4, 4);
            dst_image = unsharp(src_image);
        } else if (action == "--gray-world") {
            check_argc(argc, 4, 4);
            dst_image = gray_world(src_image);
        } else if (action == "--resize") {
            check_argc(argc, 5, 5);
            //double scale = read_value<double>(argv[4]);
            // dst_image = resize(src_image, scale);
        }  else if (action == "--custom") {
            check_argc(argc, 5, 5);
            // Matrix<double> kernel = parse_kernel(argv[4]);
            // Function custom is useful for making concrete linear filtrations
            // like gaussian or sobel. So, we assume that you implement custom
            // and then implement concrete filtrations using this function.
            // For example, sobel_x may look like this:
            // sobel_x (...) {
            //    Matrix<double> kernel = {{-1, 0, 1},
            //                             {-2, 0, 2},
            //                             {-1, 0, 1}};
            //    return custom(src_image, kernel);
            // }
            // dst_image = custom(src_image, kernel);
        } else if (action == "--autocontrast") {
            check_argc(argc, 4, 5);
            double fraction = 0.0;
            if (argc == 5) {
                fraction = read_value<double>(argv[4]);
                check_number("fraction", fraction, 0.0, 0.4);
            }
            // dst_image = autocontrast(src_image, fraction);
        } else if (action == "--gaussian" || action == "--gaussian-separable") {
            check_argc(argc, 5, 6);
            double sigma = read_value<double>(argv[4]);
            check_number("sigma", sigma, 0.1, 100.0);
            int radius = 3 * sigma;
            if (argc == 6) {
                radius = read_value<int>(argv[5]);
                check_number("radius", radius, 1);
            }
            if (action == "--gaussian") {
                dst_image = gaussian(src_image, sigma, radius);
            } else {
                dst_image = gaussian_separable(src_image, sigma, radius);
            }
        } else if (action == "--canny") {
            check_argc(6, 6);
            int threshold1 = read_value<int>(argv[4]);
            check_number("threshold1", threshold1, 0, 360);
            int threshold2 = read_value<int>(argv[5]);
            check_number("threshold2", threshold2, 0, 360);
            if (threshold1 >= threshold2)
                throw string("threshold1 must be less than threshold2");
             dst_image = canny(src_image, threshold1, threshold2);
        } else if (action == "--align") {
            check_argc(argc, 4, 6);
            if (argc == 5) {
                string postprocessing(argv[4]);
                if (postprocessing == "--gray-world" ||
                    postprocessing == "--unsharp") {
                    check_argc(argc, 5, 5);
                    //dst_image = align(src_image, postprocessing);
                } else if (postprocessing == "--autocontrast") {
                    double fraction = 0.0;
                    if (argc == 6) {
                        fraction = read_value<double>(argv[5]);
                        check_number("fraction", fraction, 0.0, 0.4);
                    }
                    // dst_image = align(src_image, postprocessing, fraction);
                } else {
                    throw string("unknown align option ") + postprocessing;
                }
            } else {
                // dst_image = align(src_image);
            }
        } else {
            throw string("unknown action ") + action;
        }
        save_image(dst_image, argv[2]);
    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        cerr << "For help type: " << endl << argv[0] << " --help" << endl;
        return 1;
    }
}
