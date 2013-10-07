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
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        uint r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                sum_r += r;
                sum_g += g;
                sum_b += b;
            }
        }
        auto norm = size * size;
        sum_r /= norm;
        sum_g /= norm;
        sum_b /= norm;
        return make_tuple(sum_r, sum_g, sum_b);
    }
    // Radius of neighbourhoud, which is passed to that operator
    static const int radius = 1;
};

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

int main(int argc, char **argv)
{
    // Image = Matrix<tuple<uint, uint, uint>>
    // tuple is from c++ 11 standard
    Image img = load_image(argv[1]);
    Image img3;
    img3 = img;
    Image img2 = img3.unary_map(BoxFilterOp());
    save_image(img2, argv[2]);
}
