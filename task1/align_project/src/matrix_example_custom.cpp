#include <iostream>
#include "matrix.h"
#include "io.h"

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
    BoxFilterOp(Matrix<int> &kernel) : ker(kernel) {}
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
    Matrix<int> ker;
    static const int radius = 1;
};

int main(int argc, char **argv)
{
    // Image = Matrix<tuple<uint, uint, uint>>
    // tuple is from c++ 11 standard
    Matrix<int> ker = {{-1, 0, 1},
                      {-2, 0, 2},
                      {-1, 0, 1}};
                      /*{{1, 2, 1},
                      {0, 0, 0},
                      {-1, -2, -1}};*/
    Image img = load_image(argv[1]);
    Image img3;
    img3 = img;
    Matrix<int> mat(3,3);
    Image img2 = img3.unary_map(BoxFilterOp(ker));
    save_image(img2, argv[2]);
}
