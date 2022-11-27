#include "function.h"
#include <vector>
#include <iostream>
using namespace std;

void Matrix::SVD()
{
    vector<Matrix> svd;
    Matrix& A = *this;
    if (A.transpose() == A)
    {
        pair<Matrix, Matrix> SU = Eigen_Values_Vectors_symmetrical(A);
        Matrix Sigma = SU.first;
        Matrix U = SU.second;
        Matrix V = U.transpose();

        cout << "U:" << endl << U << endl << endl << "S: " << endl << Sigma << endl << endl << "V: " << endl << V << endl << endl;

        Matrix SVD = U * Sigma * V;
        cout << "S*V*D: " << endl << SVD;
    }
}