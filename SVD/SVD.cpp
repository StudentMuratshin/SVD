#include <iostream>
#include "function.h"
#include <vector>

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
        cout << "S*V*D: " << endl << CHECK_SVD;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}