#include "function.h"
#include <vector>
#include <iostream>
using namespace std;


void Matrix::SVD(Matrix& U, Matrix& Sigma, Matrix& V)
{
    Matrix& A = *this;
    if (A.transpose() == A)
    {
        pair<Matrix, Matrix> SU = Eigen_Values_Vectors_symmetrical(A);
        Sigma = SU.first;
        U = SU.second;
        V = U.transpose();
    }
    else throw std::runtime_error("Asymmetrical!");
}