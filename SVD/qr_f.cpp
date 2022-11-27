#include "function.h"
#include <vector>
#include <iostream>
using namespace std;

pair<Matrix, Matrix> Matrix::QR()
{
    Matrix& a = *this;
    if (a.getH() != a.getW())
    {
        Matrix a_sq = Make_square(a);
        Matrix Q = Gram_Schmidt(a_sq);
        Matrix R = Q.transpose() * a;
        return pair<Matrix, Matrix>(Q, R);
    }
    else
    {
        Matrix Q = Gram_Schmidt(a);
        Matrix R = Q.transpose() * a;
        return pair<Matrix, Matrix>(Q, R);
    }
}

pair<Matrix, Matrix> Eigen_Values_Vectors_symmetrical(const Matrix B)
{
    if (B.transpose() == B)
    {
        Matrix A = B;
        Matrix V = V.I(B.getW());
        int n = 0;
        double Max_prev = 0, Max = 0;
        do {
            n++;
            Max_prev = A.getDiagonal().getMaxVal();
            pair<Matrix, Matrix> QR = A.QR();
            Matrix Q = QR.first;
            Matrix R = QR.second;
            V = V * Q;
            A = R * Q;
            Max = A.getDiagonal().getMaxVal();
        } while (abs(Max - Max_prev) >= 1e-13);
        for (int i = 0; i < A.getH(); i++)
            for (int j = 0; j < A.getW(); j++)
            {
                A.ref(i, j) = (i != j) ? 0 : A.get(i, j);
            }
        return pair<Matrix, Matrix>(A, V);
    }
    else throw std::runtime_error("Asymmetrical!");
}


Matrix Eigen_Values(const Matrix B)
{
    Matrix A = B;
    Matrix V = V.I(B.getW());
    int n = 0;
    double Max_prev = 0, Max = 0;
    do {
        n++;
        Max_prev = A.getDiagonal().getMaxVal();
        pair<Matrix, Matrix> QR = A.QR();
        Matrix Q = QR.first;
        Matrix R = QR.second;
        V = V * Q;
        A = R * Q;
        Max = A.getDiagonal().getMaxVal();
    } while (abs(Max - Max_prev) >= 1e-5);
    for (int i = 0; i < A.getH(); i++)
        for (int j = 0; j < A.getW(); j++)
        {
            A.ref(i, j) = (i != j) ? 0 : A.get(i, j);
        }
    return A;
}