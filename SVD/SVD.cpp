#include <iostream>
#include "function.h"
#include <vector>
#include <complex>

using namespace std;
int main()
{
    try {
        Matrix mat = { 2,2,{
            -2,1,
            1,1
        } };

        Matrix U{ mat.getH(),mat.getH() };
        Matrix S{ mat.getH(),mat.getW() };
        Matrix V{ mat.getW(),mat.getW() };

        mat.SVD(U,S,V);

        cout << "U:" << endl << U << endl << endl << "S: " << endl << S << endl << endl << "V: " << endl << V << endl << endl;

        Matrix CHECK_SVD = U * S * V;
        cout << "S*V*D: " << endl << CHECK_SVD << endl << endl;

        Matrix A = { 3,3,{
            3,2,2,
            2,3,-2,
            1,1,-4
        } };
        Eigen_Values(A);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}