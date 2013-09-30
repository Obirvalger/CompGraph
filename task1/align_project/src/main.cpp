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

class Gauss
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
};

Image gaussian(const Image &src_image,double sigma = 1.4, uint radius = 1)
{
    return src_image.unary_map(Gauss(sigma,radius));
}

template<typename T>
class Custom
{
public:
    Custom(Matrix<T> &kernel) : ker(kernel), radius(ker.n_rows) {}
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
        int elem = 0, div = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                elem = ker(i,j);
                div += elem;
                //cout<<gauss<<div;
                sum_r += r * elem;
                sum_g += g * elem;
                sum_b += b * elem;
            }
        }
        //auto norm = size * size;
        if (div==0)
        div = 1;
        sum_r /= div;
        sum_g /= div;
        sum_b /= div;
        return make_tuple(sum_r, sum_g, sum_b);
    }
    // Radius of neighbourhoud, which is passed to that operator
    Matrix<T> ker;
    int radius;
};


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
            // dst_image = sobel_x(src_image);
        } else if (action == "--sobel-y") {
            check_argc(argc, 4, 4);
            // dst_image = sobel_y(src_image);
        } else if (action == "--unsharp") {
            check_argc(argc, 4, 4);
            // dst_image = unsharp(src_image);
        } else if (action == "--gray-world") {
            check_argc(argc, 4, 4);
            // dst_image = gray_world(src_image);
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
                // dst_image = gaussian_separable(src_image, sigma, radius);
            }
        } else if (action == "--canny") {
            check_argc(6, 6);
            int threshold1 = read_value<int>(argv[4]);
            check_number("threshold1", threshold1, 0, 360);
            int threshold2 = read_value<int>(argv[5]);
            check_number("threshold2", threshold2, 0, 360);
            if (threshold1 >= threshold2)
                throw string("threshold1 must be less than threshold2");
            // dst_image = canny(src_image, threshold1, threshold2);
        } else if (action == "--align") {
            check_argc(argc, 4, 6);
            if (argc == 5) {
                string postprocessing(argv[4]);
                if (postprocessing == "--gray-world" ||
                    postprocessing == "--unsharp") {
                    check_argc(argc, 5, 5);
                    // dst_image = align(src_image, postprocessing);
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
