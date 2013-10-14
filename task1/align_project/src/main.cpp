#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "matrix.h"
#include "io.h"
#include <initializer_list>
#include <limits>

using std::tuple;
using std::vector;
using std::min;
using std::max;
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


bool pw(int x, int y)
{
    for(uint i = y; i > 0; --i) {
        x *= x;
    }
    return x;
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

Matrix<double> parse_kernel(string str)
{
    // Kernel parsing implementation here
    const char *cs = str.c_str();
    char c = cs[0];
    uint nr = 0, nc = 0, cur = 0, i = 0, j = 0;
    if (str.find_first_not_of("+-0123467895.,Ie") != string::npos)
        { cout << "lol"; throw "bad kernel";}

    while(c != '\0')
         {
        while ((c != ',') && (c != 'I')) {
            c = cs[++i];
        }

        if (c == ',') cur++;
            if (c == 'I') {
                nr++;
                if (nr == 1)
                    nc = cur;
                if (nc != cur)
        { cout << "lol too\n"; throw "bad kernel";}
                cur = 0;
            }
        c = cs[++i];
    }
    nc++;
    cout<<nr<<" "<<nc<<" "/*<<int(cs[str.length()])*/<<"gsr"<<endl;
    Matrix<double> kernel(nr,nc);
    bool f = false;
    c = cs[i = cur = j = 0];
    while(c != '\0') {
        f = false;
        while ((c != ',') && (c != 'I')) {
            cout << c << " ";
            if (!f) {
                kernel(j,cur) = atof(&cs[i]);
                f = true;
            }
            c = cs[++i];
            //kernel(i,j) = atof(&cs[cur]);
        }
        if (c == ',') {
           // kernel(cur,j) = atof(&cs[i + 1])
            f = false;
            cur++;
        }
        if (c == 'I') {
            j++;
            cur = 0;
        }
        c = cs[++i];
    }
    cout<<endl<<endl<<kernel<<endl;
    //string s = "2,3,5.4;6,7,2";
    //scout<<atof(cs)<<" "<<cs<<endl;*/

    return kernel;
}

Image resample(const Image &a, int k)
{
    int oldh = a.n_rows, oldw = a.n_cols, neww, newh;
    if (k > 0) {
        neww = oldw / k;
        newh = oldh / k;
    } else {
        neww = oldw * abs(k);
        newh = oldh * abs(k);
    }
    Image img(newh,neww);
    int i, j;
    int h, w;
    float t;
    float u;
    float tmp;
    float d1, d2, d3, d4;
    uint r1, r2, r3, r4, g1, g2, g3, g4, b1, b2, b3, b4, r, g, b;

    for (j = 0; j < newh; j++) {
        tmp = static_cast<float>(j) / static_cast<float>(newh - 1) * (oldh - 1);
        h = floor(tmp);
        if (h < 0) {
            h = 0;
        } else {
            if (h >= oldh - 1) {
                h = oldh - 2;
            }
        }
        u = tmp - h;

        for (i = 0; i < neww; i++) {

            tmp = static_cast<float>(i) / static_cast<float>(neww - 1) * (oldw - 1);
            w = floor(tmp);
            if (w < 0) {
                w = 0;
            } else {
                if (w >= oldw - 1) {
                    w = oldw - 2;
                }
            }
            t = tmp - w;
            d1 = (1 - t) * (1 - u);
            d2 = t * (1 - u);
            d3 = t * u;
            d4 = (1 - t) * u;
            tie(r1,g1,b1) = a(h,w);
            tie(r2,g2,b2) = a(h,w + 1);
            tie(r3,g3,b3) = a(h + 1,w + 1);
            tie(r4,g4,b4) = a(h + 1,w);
            b = b1 * d1 + b2 * d2 + b3 * d3 + b4 * d4;
            g = g1 * d1 + g2 * d2 + g3 * d3 + g4 * d4;
            r = r1 * d1 + r2 * d2 + r3 * d3 + r4 *d4;
            r = min(int(255),int(r));
            g = min(int(255),int(g));
            b = min(int(255),int(b));
            img(j,i) = make_tuple(r,g,b);
        }
    }
    return img;
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

//Makes image bigger for better filtering
Image bigimg(const Image &m, uint r)
{
    //uint r = 2 * radius + 1;
    uint nr = m.n_rows, nc = m.n_cols, d = 2 * r;
    Image img(nr + d, nc + d);
    for (uint i = 0; i < nr; ++i) {
        for (uint j = 0; j < nc; ++j) {
            img(i + r,j + r) = m(i,j);
        }
    }
    for (uint i = 0; i < r; ++i) {
        for (uint j = r; j < nc + r; ++j) {
            img(r - 1 - i,j) = m(i,j - r);
            img(nr + r + i,j) = m(nr - 1 - i,j - r);
            //cout<<i<<" "<<j<<" "<<nr<<" "<<nc<<endl;
        }
        //cout<<i<<endl;
    }
    //cout<<"sefse"<<endl;
    for (uint i = 0; i < nr + d; ++i) {
        for (uint j = 0; j < r; ++j) {
            img(i,r - 1 - j) = img(i,j + r);
            img(i,nc + r + j) = img(i,nc + r - 1 - j);
        }
    }
    //cout<<"sefse"<<endl;
    return img;
}

Image smallimg(const Image &m, uint r)
{
    //uint r = 2 * radius + 1;
    uint d = 2 * r, nr = m.n_rows - d, nc = m.n_cols - d;
    Image img(nr, nc);
    for (uint i = 0; i < nr; ++i) {
        for (uint j = 0; j < nc; ++j) {
            img(i,j) = m(i + r + 1,j + r + 1);
        }
    }
    //cout<<"sefse"<<endl;
    return img;
}

template<typename T>
class Custom
{
public:
    Custom(Matrix<T> &kernel, bool saveimg = true) : ker(kernel), radius((ker.n_rows - 1) / 2), save(saveimg) {}
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

Image custom(const Image &src_image, Matrix<double> kernel)
{
    return src_image.unary_map(Custom<double>(kernel));
}

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
    //cout<<radius<<" "<<(kernel.n_cols - 1) / 2<<endl;
    return src_image.unary_map(Custom<double>(kernel));
    //return smallimg(bigimg(src_image,radius).unary_map(Gauss(1.4, radius)), radius);
    //return smallimg(custom(bigimg(src_image,radius),kernel,radius), radius);
    //return smallimg(bigimg(src_image,radius), radius);1
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

Image autocontrast(const Image &src_image, double fraction)
{
    uint nr = src_image.n_rows, nc = src_image.n_cols, r, g, b;
    double y, ymin = 0, ymax = 0;
    Image img(nr,nc);
    vector<int> gist(256);
    for (uint i = 0; i < nr; ++i) {
        for (uint j = 0; j < nc; ++j) {
            tie(r,g,b) = src_image(i,j);
            y = 0.2125 * r + 0.7154 * g + 0.0721 * b;
            gist[static_cast<uint>(y)]++;
        }
    }
    for (uint i = 255 * fraction; i < 255 - 255 * fraction; ++i) {
        y = gist[i];
        if (i == 1) {
            ymin = gist[i];
            ymax = gist[i];
        }
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }
    for (uint i = 0; i < nr; ++i) {
        for (uint j = 0; j < nc; ++j) {
            tie(r,g,b) = src_image(i,j);
            y = gist[i];
            r *= (y - ymin) * 255 / (ymax - ymin) / y;
            g *= (y - ymin) * 255 / (ymax - ymin) / y;
            b *= (y - ymin) * 255 / (ymax - ymin) / y;
            img(i,j) = make_tuple(r,g,b);
        }
    }
    return img;
}

Image canny(const Image &src_image, int th1, int th2)
{
    Image m = gaussian(src_image,1.4,2);
    Image imgx = sobel_x(m,false),//m.unary_map(Custom<double>(sobelx)),
        imgy = sobel_y(m,false);//m.unary_map(Custom<double>(sobely));
    vector<vector<int>> bord;
    vector<int> elem(2);
    Matrix<int> bordmap(m.n_rows, m.n_cols);
    Matrix<double> grad(m.n_rows, m.n_cols), theta(m.n_rows, m.n_cols);
    unsigned int rx, ry, nr = m.n_rows, nc = m.n_cols, top, bottom, left, right;

    //find gradient
    for (uint i = 0; i < nr; ++i) {
        for (uint j = 0; j < nc; ++j) {
            tie(rx,rx,rx) = imgx(i,j); tie(ry,ry,ry) = imgy(i,j);
            grad(i,j) = sqrt(rx*rx+ry*ry);
            theta(i,j) = atan2(rx, ry);
        }
    }

    //border map making with great pixels
    for (uint i = 1; i < nr - 1; ++i) {
        for (uint j = 1; j < nc - 1; ++j) {
            switch (radtoplace(theta(i,j))) {
                case 1:
                    if ((grad(i,j) <=  grad(i,j + 1)) || (grad(i,j) <= grad(i,j - 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0; else
                    if (grad(i,j) > th2) {
                        bordmap(i,j) = 1;
                        elem[0] = i;
                        elem[1] = j;
                        bord.push_back(elem);
                    } else bordmap(i,j) = -1;

                case 2:
                    if ((grad(i,j) <=  grad(i + 1,j + 1)) || (grad(i,j) <= grad(i - 1,j - 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0; else
                    if (grad(i,j) > th2) {
                        bordmap(i,j) = 1;
                        elem[0] = i;
                        elem[1] = j;
                        bord.push_back(elem);
                    } else bordmap(i,j) = -1;
                case 3:
                    if ((grad(i,j) <=  grad(i + 1,j)) || (grad(i,j) <= grad(i - 1,j)) || (grad(i,j) < th1))
                        grad(i,j) = 0; else
                    if (grad(i,j) > th2) {
                        bordmap(i,j) = 1;
                        elem[0] = i;
                        elem[1] = j;
                        bord.push_back(elem);
                    } else bordmap(i,j) = -1;
                case 4:
                    if ((grad(i,j) <=  grad(i + 1,j + 1)) || (grad(i,j) <= grad(i - 1,j - 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0; else
                    if (grad(i,j) > th2) {
                        bordmap(i,j) = 1;
                        elem[0] = i;
                        elem[1] = j;
                        bord.push_back(elem);
                    } else bordmap(i,j) = -1;
                case 5:
                    if ((grad(i,j) <=  grad(i,j - 1)) || (grad(i,j) <= grad(i,j + 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0; else
                    if (grad(i,j) > th2) {
                        bordmap(i,j) = 1;
                        elem[0] = i;
                        elem[1] = j;
                        bord.push_back(elem);
                    } else bordmap(i,j) = -1;
                case 6:
                    if ((grad(i,j) <=  grad(i - 1,j - 1)) || (grad(i,j) <= grad(i + 1,j + 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0; else
                    if (grad(i,j) > th2) {
                        bordmap(i,j) = 1;
                        elem[0] = i;
                        elem[1] = j;
                        bord.push_back(elem);
                    } else bordmap(i,j) = -1;
                case 7:
                    if ((grad(i,j) <=  grad(i - 1,j)) || (grad(i,j) <= grad(i + 1,j)) || (grad(i,j) < th1))
                        grad(i,j) = 0; else
                    if (grad(i,j) > th2) {
                        bordmap(i,j) = 1;
                        elem[0] = i;
                        elem[1] = j;
                        bord.push_back(elem);
                    } else bordmap(i,j) = -1;
                case 8:
                    if ((grad(i,j) <=  grad(i - 1,j + 1)) || (grad(i,j) <= grad(i + 1,j - 1)) || (grad(i,j) < th1))
                        grad(i,j) = 0; else
                    if (grad(i,j) > th2) {
                        bordmap(i,j) = 1;
                        elem[0] = i;
                        elem[1] = j;
                        bord.push_back(elem);
                    } else bordmap(i,j) = -1;
                default: ;
            }
        }
    }

    //border map making with middle pixels
    for(uint k = 0; k < nr * nc; ++k) {
        uint i, j;
        if (k == bord.size()) break;
        i = bord[k][0]; j = bord[k][1];
        if ((j != nc) && (bordmap(i,j + 1) == -1)) {
            bordmap(i,j + 1) = 1;
            elem[0] = i;
            elem[1] = j + 1;
            bord.push_back(elem);
        }
        if ((j != nc) && (i != nr) && (bordmap(i + 1,j + 1) == -1)) {
            bordmap(i + 1,j + 1) = 1;
            elem[0] = i + 1;
            elem[1] = j + 1;
            bord.push_back(elem);
        }
        if ((i != 0) && (bordmap(i - 1,j) == -1)) {
            bordmap(i - 1,j) = 1;
            elem[0] = i - 1;
            elem[1] = j;
            bord.push_back(elem);
        }
        if ((i != 0)  && (j != 0) && (bordmap(i - 1,j - 1) == -1)) {
            bordmap(i - 1,j - 1) = 1;
            elem[0] = i - 1;
            elem[1] = j - 1;
            bord.push_back(elem);
        }
        if ((j != 0) && (bordmap(i,j - 1) == -1)) {
            bordmap(i,j - 1) = 1;
            elem[0] = i;
            elem[1] = j - 1;
            bord.push_back(elem);
        }
        if ((j != 0) && (i != nr) && (bordmap(i + 1,j - 1) == -1)) {
            bordmap(i + 1,j - 1) = 1;
            elem[0] = i + 1;
            elem[1] = j - 1;
            bord.push_back(elem);
        }
        if ((i != nr) && (bordmap(i + 1,j) == -1)) {
            bordmap(i + 1,j) = 1;
            elem[0] = i + 1;
            elem[1] = j;
            bord.push_back(elem);
        }
        if ((j != nc) && (i != nr) && (bordmap(i + 1,j + 1) == -1)) {
            bordmap(i + 1,j + 1) = 1;
            elem[0] = i + 1;
            elem[1] = j + 1;
            bord.push_back(elem);
        }
    }
    for (uint i = 0; i < nr; ++i) {
        for (uint j = 0; j < nc; ++j) {
            if(bordmap(i,j) != 1)
                bordmap(i,j) = 0;
        }
    }

    for (uint i = 0; i < nr; ++i) {
        for (uint j = 0; j < nc; ++j) {
            imgx(i,j) = make_tuple(bordmap(i,j) * 255, bordmap(i,j) * 255, bordmap(i,j) * 255);
        }
    }

    //cutting borders
    uint sum = 0, mx = 0, pos1 = 0, pos2 = 0, k = 2;
    for (uint i = 0; i < nr / 100 * 5; ++i) {
        for (uint j = 0; j < nc; ++j) {
            sum += bordmap(i,j);
        }
        if (sum > mx) {
            mx = sum;
            pos1 = i;
            //cout<<"top1"<<endl;
        }
        sum = 0;
    }
    sum = mx = 0;
    for (uint i = 0; i < nr / 100 * 5; ++i) {
        for (uint j = 0; j < nc; ++j) {
            sum += bordmap(i,j);
        }
        if ((i > pos1 + k) && (sum > mx)) {
            mx = sum;
            pos2 = i;
            //cout<<"top2"<<endl;
        }
        sum = 0;
    }
    top = max(pos1,pos2);

    sum = mx = 0;
    pos1 = pos2 = 2;
    for (uint j = 2; j < nc / 100 * 5; ++j) {
        for (uint i = 0; i < nr; ++i) {
            sum += bordmap(i,j);
        }
        if (sum > mx) {
            mx = sum;
            pos1 = j;
        }
        sum = 0;
    }
    sum = mx = 0;
    for (uint j = 2; j < nc / 100 * 5; ++j) {
        for (uint i = 0; i < nr; ++i) {
            sum += bordmap(i,j);
        }
        if ((j > pos1 + k) && (sum > mx)) {
            mx = sum;
            pos2 = j;
        }
        sum = 0;
    }
    left = max(pos1,pos2);
    //cout<<pos1<<" "<<pos2<<endl;

    sum = mx = 0;
    pos1 = pos2 = nr - 1;
    for (uint i = nr - 1; i >= nr / 100 * 95; --i) {
        for (uint j = 0; j < nc; ++j) {
            sum += bordmap(i,j);
        }
        if (sum > mx) {
            mx = sum;
            pos1 = i;
        }
        sum = 0;
    }
    sum = mx = 0;
    for (uint i = nr - 1; i >= nr / 100 * 95; --i) {
        for (uint j = 0; j < nc; ++j) {
            sum += bordmap(i,j);
        }
        if ((i < pos1 - k) && (sum > mx)) {
            mx = sum;
            pos2 = i;
        }
        sum = 0;
    }
    bottom = min(pos1,pos2);

    sum = mx = 0;
    pos1 = pos2 = nc - 3;
    for (uint j = nc - 3; j >= nc / 100 * 95; --j) {
        for (uint i = 0; i < nr; ++i) {
            sum += bordmap(i,j);
        }
        if (sum > mx) {
            mx = sum;
            pos1 = j;
        }
        sum = 0;
    }
    sum = mx = 0;
    for (uint j = nc - 3; j >= nc / 100 * 95; --j) {
        for (uint i = 0; i < nr; ++i) {
            sum += bordmap(i,j);
        }
        if ((j < pos1 - k) && (sum > mx)) {
            mx = sum;
            pos2 = j;
        }
        sum = 0;
    }
    right = min(pos1,pos2);

    //cout<<top<<" "<<left<<" "<<bottom<<" "<<right<<" "<<endl;

    //making cutted image
    Image img(bottom - top, right - left);
    for (uint i = top + 1; i < bottom; ++i) {
        for (uint j = left + 1; j < right; ++j) {
            img(i - top - 1,j - left - 1) = src_image(i,j);
        }
    }
    return img;
}

int mse(Image img1, Image img2, int lim = 15)
{
    int ans, a1 = 0, a2 = 0;
    uint nr = img1.n_rows, nc = img1.n_cols, val = 0, mv1 = 0, mv2 = 0, r1,r2;
    for (int k = 0; k <= lim; ++k) {
        for (uint i = 0; i < nr - k; ++i) {
            for (uint j = 0; j < nc; ++j) {
            tie(r1,r1,r1) = img1(i,j);
            tie(r2,r2,r2) = img2(i + k,j);
            val += (r1 - r2) * (r1 - r2);
            }
        }
        val /= nr * nc;
        if (k == 0) mv1 = val;
        if (val < mv1) {
            mv1 = val;
            a1 = -k;
        }
        val = 0;
    }
    for (int k = 0; k <= lim; ++k) {
        for (uint i = 0; i < nr - k; ++i) {
            for (uint j = 0; j < nc; ++j) {
            tie(r1,r1,r1) = img1(i + k,j);
            tie(r2,r2,r2) = img2(i,j);
            val += (r1 - r2) * (r1 - r2);
            }
        }
        val /= nr * nc;
        if (k == 0) mv2 = val;
        if (val < mv2) {
            mv2 = val;
            a2 = k;
        }
        val = 0;
    }
    if (mv1 < mv2) {
        ans = a1;
    } else {
        ans = a2;
    }
    return ans;
}

//Cross-correlation
int cross(Image img1, Image img2, int lim = 15)
{
    int ans, a1 = 0, a2 = 0;
    uint nr = img1.n_rows, nc = img1.n_cols, val = 0, mv1 = 0, mv2 = 0, r1,r2;
    for (int k = 0; k <= lim; ++k) {
        for (uint i = 0; i < nr - k; ++i) {
            for (uint j = 0; j < nc; ++j) {
            tie(r1,r1,r1) = img1(i,j);
            tie(r2,r2,r2) = img2(i + k,j);
            val += r1 * r2;
            }
        }
        if (k == 0) mv1 = val;
        if (val > mv1) {
            mv1 = val;
            a1 = -k;
        }
        val = 0;
    }
    for (int k = 0; k <= lim; ++k) {
        for (uint i = 0; i < nr - k; ++i) {
            for (uint j = 0; j < nc; ++j) {
            tie(r1,r1,r1) = img1(i + k,j);
            tie(r2,r2,r2) = img2(i,j);
            val += r1 * r2;
            }
        }
        if (k == 0) mv2 = val;
        if (val > mv2) {
            mv2 = val;
            a2 = k;
        }
        val = 0;
    }
    if (mv1 > mv2) {
        ans = a1;
    } else {
        ans = a2;
    }
    return ans;
}

Image align(const Image &src_image, string postprocessing = "", double fraction = 0, bool save = true)
{
    Image ubimg = canny(src_image,100,300);
    uint nr = ubimg.n_rows - ubimg.n_rows % 3, nc = ubimg.n_cols;
    cout<<"separating source image to RGB chanels"<<endl;
    Image imgr(nr / 3,nc), imgg(nr / 3,nc), imgb(nr / 3,nc);

    //separating source image to RGB chanels
    for (uint i = 0; i < nr / 3; ++i){
        for (uint j = 0; j < nc; ++j) {
            imgb(i,j) = ubimg(i,j);
        }
    }
    for (uint i = nr / 3; i < nr / 3 * 2; ++i){
        for (uint j = 0; j < nc; ++j) {
            imgg(i - nr / 3,j) = ubimg(i,j);
        }
        //cout<<i<<endl;
    }
    cout<<"Red begins"<<endl;
    for (uint i = nr / 3 * 2; i < nr; ++i){
        for (uint j = 0; j < nc; ++j) {
            imgr(i - nr / 3 * 2,j) = ubimg(i,j);
        }
        //cout<<i<<endl;
    }

    /*for (uint i = nr - 1; i >= nr / 3 * 2; --i){
        cout<<nr<<" "<<nr / 3 * 2<<" "<<i<<endl;
        for (uint j = 0; j < nc; ++j) {
            imgr(i - nr / 3 * 2,j) = ubimg(i,j);
        }
    }*/
    //Image resimg = resample(src_image, 2);
    //int gb = mse(imgg,imgb, 30), gr = mse(imgg,imgr,30);
    int gb = 15, gr = 15;
    uint r, g, b, k = 2, num = 0;
    Image img(nr / 3, nc);
    uint min_side = min(nr,nc);
    uint l = min_side;
    while (l > 400 * k) {
        l /= k;
        cout<<l<<endl;
        num++;
    }
    //cout<<num<<endl;
    for (uint i = num; i > 0; --i) {
        //cout<<gr<<" "<<gb<<endl;
        gb = k * mse(resample(imgg, pw(k,i)),resample(imgb, pw(k,i)),abs(gb));
        gr = k * mse(resample(imgg, pw(k,i)),resample(imgr, pw(k,i)),abs(gr));
        cout<<gr<<" "<<gb<<endl;
    }
    //cout<<gr<<" "<<gb<<endl;
    gr = mse(imgg,imgr,abs(gr));
    gb = mse(imgg,imgb,abs(gb));
    cout<<gr<<" "<<gb<<endl;
    for (uint i = 0; i < nr / 3; ++i) {
        for (uint j = 0; j < nc; ++j) {
            tie(g,g,g) = imgg(i,j);
            img(i,j) = make_tuple(0,g,0);
        }
    }

    if (gb >= 0) {
        for (uint i = gb; i < nr / 3; ++i) {
            for (uint j = 0; j < nc; ++j) {
                tie(r,g,b) = img(i,j);
                tie(b,b,b) = imgb(i - gb,j);
                img(i,j) = make_tuple(0,g,b);
            }
        }
    } else {
        for (uint i = 0; i < nr / 3 - abs(gb); ++i) {
            for (uint j = 0; j < nc; ++j) {
                tie(r,g,b) = img(i,j);
                tie(b,b,b) = imgb(i - gb,j);
                img(i,j) = make_tuple(0,g,b);
            }
        }
    }

    if (gr >= 0) {
        for (uint i = gr; i < nr / 3; ++i) {
            for (uint j = 0; j < nc; ++j) {
                tie(r,g,b) = img(i,j);
                tie(r,r,r) = imgr(i - gr,j);
                img(i,j) = make_tuple(r,g,b);
            }
        }
    } else {
        for (uint i = 0; i < nr / 3 - abs(gr); ++i) {
            for (uint j = 0; j < nc; ++j) {
                tie(r,g,b) = img(i,j);
                tie(r,r,r) = imgr(i - gr,j);
                img(i,j) = make_tuple(r,g,b);
            }
        }
    }

    /*for (uint i = 0; i < nr / 3; ++i) {
        for (uint j = 0; j < nc; ++j) {
            if (gb >= 0) {
            tie(g,g,g) = imgg(i,j);
            tie(b,b,b) = imgb(i + gb,j);
            } else {
            tie(g,g,g) = imgg(i,j);
            tie(b,b,b) = imgb(i - gb,j);
            }
            if (gr >= 0) {
            tie(g,g,g) = imgg(i,j);
            tie(r,r,r) = imgr(i + gr,j);
            } else {
            tie(g,g,g) = imgg(i + -gr,j);
            tie(r,r,r) = imgr(i,j);
            }
            img(i,j) = make_tuple(r,g,b);
        }
    }*/
    if (postprocessing == "--gray-world")
        return gray_world(img);
    if (postprocessing == "--unsharp")
        return unsharp(img);
    if (postprocessing == "--autocontrast")
        return autocontrast(img,fraction);
    //img = resample(src_image, -8);
    return img;
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
            Matrix<double> kernel = parse_kernel(argv[4]);
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
            dst_image = custom(src_image, kernel);
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
                    dst_image = align(src_image, postprocessing);
                } else if (postprocessing == "--autocontrast") {
                    double fraction = 0.0;
                    if (argc == 6) {
                        fraction = read_value<double>(argv[5]);
                        check_number("fraction", fraction, 0.0, 0.4);
                    }
                    dst_image = align(src_image, postprocessing, fraction);
                } else {
                    throw string("unknown align option ") + postprocessing;
                }
            } else {
                dst_image = align(src_image);
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
