#include <iostream>
#include "matrix.h"
#include "io.h"

//using namespace std;
using std::cout;
using std::endl;

using std::tuple;
using std::get;
using std::tie;
using std::make_tuple;

// Matrix usage example
// Also see: matrix.h, matrix.hpp for comments on how filtering works

class BoxFilterOp
{
public:
    BoxFilterOp(double s = 1.4, int r = 1) : sigma(s), radius(r) {}
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        double gauss[size][size];
        uint r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        double div = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                gauss[i][j] = (1/(sqrt(M_PI*2))*sigma)*exp((i*i+j*j)/(2*sigma*sigma)*(-1.0));
                div += gauss[i][j];
            }
        }
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                gauss[i][j] /= div;
            }
        }
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                sum_r += r * gauss[i][j];
                sum_g += g * gauss[i][j];
                sum_b += b * gauss[i][j];
            }
        }
        return make_tuple(sum_r - 0, sum_g - 0, sum_b - 0);
    }
    // Radius of neighbourhoud, which is passed to that operator
    double sigma;
    int radius;
};

int main(int argc, char **argv)
{
    // Image = Matrix<tuple<uint, uint, uint>>
    // tuple is from c++ 11 standard
    Image img = load_image(argv[1]);
    Image img3;
    img3 = img;
    Image img2 = img3.unary_map(BoxFilterOp(1.4, 7));
    save_image(img2, argv[2]);
}
