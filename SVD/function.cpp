#include "function.h"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

using namespace std;



Matrix Gram_Schmidt(Matrix& arr)
{
    Matrix b{ arr.getH(), arr.getW() };
    Matrix e{ arr.getH(), arr.getW() };
    for (int j = 0; j < arr.getW(); j++)
    {
        b.ref(0, j) = arr.get(j, 0);
        e.ref(0, j) = arr.get(j, 0) / pow(Dot_product(arr.getColumn(0), arr.getColumn(0)), 1. / 2.);
    }
    Matrix term = Matrix{ 1, arr.getH() };
    for (int i = 1; i < arr.getH(); i++)
    {
        term = Matrix{ 1, arr.getH() };

        for (int j = 0; j < i; j++)
        {
            term = term - b.getRow(j) * (Dot_product(b.getRow(j), arr.getColumn(i)) / Dot_product(b.getRow(j), b.getRow(j)));
        }
        for (int j = 0; j < b.getW(); j++)
        {
            double x = arr.get(j, i), y = term.get(0, j);
            b.ref(i, j) = x + y;
        }
        for (int j = 0; j < b.getW(); j++)
        {
            double n = Dot_product(b.getRow(i), b.getRow(i));
            double digit = b.get(i, j) / sqrt(n);
            e.ref(i, j) = ((abs(digit) < 1e-10) ? 0 : digit);
            //e.ref(i, j) = digit;
        }
    }
    return e.transpose();
}

double Dot_product(const Matrix a, const Matrix b)
{
    Matrix t = a * b.transpose();
    return t.get(0, 0);
}


Matrix::Matrix(int h, int w, const double* coefficients) {
    size.w = w;
    size.h = h;
    vec.assign(coefficients, coefficients + w * h);
}
Matrix::Matrix(int h, int w) {
    size.w = w;
    size.h = h;
    vec.resize(w * h, 0);
}
Matrix::Matrix(int h, int w, std::initializer_list<double> coefficients) {
    size.w = w;
    size.h = h;
    vec.assign(coefficients.begin(), coefficients.end());
}

Matrix Matrix::I(int s, double a) {
    double* c = new double[s * s];
    for (int y = 0; y < s; y++) {
        for (int x = 0; x < s; x++) {
            c[y * s + x] = ((x == y) ? a : 0.0);
        }
    }
    Matrix i(s, s, c);
    delete[] c;
    return i;
}

double Matrix::get(int i, int j) const {
    return vec[size.w * i + j];
}
double& Matrix::ref(int i, int j) {
    return vec[size.w * i + j];
}
Matrix Matrix::operator*(const Matrix& b) const {
    const Matrix& a = *this;
    if (a.size.w != b.size.h)
        throw std::runtime_error("Can't multiply matricies with incompatible sizes!");
    int h = a.size.h;
    int w = b.size.w;
    std::vector<double> v;
    for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++) {
            double sum = 0;
            for (int k = 0; k < a.size.w; k++)
                sum += a.get(i, k) * b.get(k, j);
            v.push_back(((abs(sum) < 1e-10) ? 0 : sum));
        }
    return Matrix(h, w, v.data());
}

bool Matrix::operator==(const Matrix& b) const {
    const Matrix& a = *this;
    if (a.size.w != b.size.h)
        throw std::runtime_error("Can't multiply matricies with incompatible sizes!");
    int h = a.size.h;
    int w = b.size.w;
    std::vector<double> v;
    for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++) 
        {
            if (a.get(i, j) != b.get(i, j)) return false;
        }
    return true;
}

Matrix Matrix::operator*(const double b) const {
    const Matrix& a = *this;
    int h = a.size.h;
    int w = a.size.w;
    std::vector<double> v;
    for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++) {
            double digit = b * a.get(i, j);
            v.push_back(((abs(digit) < 1e-10) ? 0 : digit));
        }
    return Matrix(h, w, v.data());
}



Matrix Matrix::operator/(const double b) const {
    const Matrix& a = *this;
    int h = a.size.h;
    int w = a.size.w;
    std::vector<double> v;
    for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++) {
            v.push_back(a.get(i, j) / b);
        }
    return Matrix(h, w, v.data());
}

Matrix Matrix::operator+(const Matrix& other)const {
    if (other.size.w != size.w || other.size.h != size.h) {
        throw std::runtime_error("Can't add matricies with incompatible sizes!");
    }
    std::vector<double> c(size.w * size.h);
    for (int i = 0; i < size.h; i++)
        for (int j = 0; j < size.w; j++)
            c[i * size.w + j] = get(i, j) + other.get(i, j);
    return Matrix(size.h, size.w, c.data());
}

Matrix Matrix::operator-(const Matrix& other)const {
    if (other.size.w != size.w || other.size.h != size.h) {
        throw std::runtime_error("Can't add matricies with incompatible sizes!");
    }
    std::vector<double> c(size.w * size.h);
    for (int i = 0; i < size.h; i++)
        for (int j = 0; j < size.w; j++)
            c[i * size.w + j] = get(i, j) - other.get(i, j);
    return Matrix(size.h, size.w, c.data());
}

Matrix Matrix::getColumn(int col) {
    vector<double>col_data;
    for (int i = 0; i < size.w; i++) {
        col_data.push_back(get(i, col));
    }
    return Matrix{ 1, size.w, col_data.data() };
}
Matrix Matrix::getRow(int col)
{
    vector<double>col_data;
    for (int i = 0; i < size.w; i++) {
        col_data.push_back(get(col, i));
    }
    return Matrix{ 1, size.h, col_data.data() };
}

Matrix Matrix::transpose() const {
    vector<double>col_data;
    if (size.h == 1 || size.w == 1)
    {
        for (int i = 0; i < size.h; i++) {
            for (int j = 0; j < size.w; j++) {
                col_data.push_back(get(i, j));
            }
        }
        return Matrix{ size.w,size.h,col_data.data() };
    }
    else if (size.h == size.w)
    {
        for (int i = 0; i < size.h; i++) {
            for (int j = 0; j < size.w; j++) {
                col_data.push_back(get(j, i));
            }
        }
        return Matrix{ size.w,size.h,col_data.data() };
    }
    else
    {
        for (int i = 0; i < size.w; i++) {
            for (int j = 0; j < size.h; j++) {
                col_data.push_back(get(j, i));
            }
        }
        return Matrix{ size.w,size.h,col_data.data() };
    }
}

int Matrix::getW() const {
    return size.w;
}

int Matrix::getH() const {
    return size.h;
}


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

ostream& operator<<(ostream& out, Matrix& x)
{
    for (int i = 0; i < x.getH(); i++) {
        for (int j = 0; j < x.getW(); j++)
            out << x.get(i, j) << " ";
        out << endl;
    }
    return out;
}

Matrix Make_square(const Matrix b)
{
    int n = min(b.getH(), b.getW());
    Matrix res{ n , n };
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res.ref(i, j) = b.get(i, j);
        }
    }
    return res;
}

pair<Matrix,Matrix> Eigen_Values(const Matrix B)
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
        } while (abs(Max - Max_prev) >= 1e-5);
        return pair<Matrix, Matrix>(A.getDiagonal(), V);
    }
    else
    {
        Matrix A = B;
        double Max_prev = 0, Max = 0;
        do {
            Max_prev = A.getDiagonal().getMaxVal();
            pair<Matrix, Matrix> QR = A.QR();
            Matrix Q = QR.first;
            Matrix R = QR.second;
            A = R * Q;
            Max = A.getDiagonal().getMaxVal();
        } while (abs(Max - Max_prev) >= 1e-5);
        //return A.getDiagonal();
    }
}

double Matrix::getMaxVal()
{
    Matrix& A = *this;
    double Max = 0;
    for (int i = 0; i < A.getH(); i++)
    {
        for (int j = 0; j < A.getW(); j++)
        {
            if (A.get(i, j) > Max)
            {
                Max = A.get(i, j);
            }
        }
    }
    return Max;
}

Matrix Matrix::getDiagonal()
{
    vector <double> col_data;
    for (int i = 0; i < min(size.w, size.h); i++)
    {
        col_data.push_back(get(i, i));
    }
    return Matrix{ 1, min(size.w, size.h), col_data.data() };
}

Matrix solve_gauss(const Matrix ai, const Matrix bi) {
    int n = bi.getW();
    double* A = new double[ai.getH() * ai.getW()];
    double* b = new double[bi.getW()];
    double* x = new double[bi.getW()];

    for (int i = 0; i < ai.getH(); i++)
    {
        for (int j = 0; j < ai.getW(); j++)
        {
            A[i * n + j] = ai.get(i, j);
        }
    }
    for (int i = 0; i < ai.getH(); i++)
    {
        b[i] = bi.get(0, i);
    }

    for (int i = 0; i < n; ++i) {
        double pivot = A[i + i * n];
        if (abs(pivot) < DBL_EPSILON) {
            double max = 0.;
            int max_index = i;
            for (int j = i + 1; j < n; ++j) {
                if (abs(A[i + j * n]) > max) {
                    max_index = j;
                    max = A[j + j * n];
                }
            }

            if (max < DBL_EPSILON) {
                throw std::invalid_argument("det(A) = 0");
            }
            for (int k = i; k < n; ++k) {
                std::swap(A[k + max_index * n], A[k + i * n]);
            }
            std::swap(b[i], b[max_index]);
            i--;
            continue;
        }
        for (int j = i; j < n; ++j) {
            A[j + i * n] /= pivot;
        }
        b[i] /= pivot;

        for (int j = i + 1; j < n; ++j) {
            pivot = A[i + j * n];
            for (int k = 0; k < n; ++k) {
                A[k + j * n] -= pivot * A[k + n * i];
            }
            b[j] -= pivot * b[i];
        }
    }

    //backward
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[j + i * n] * x[j];
        }
    }
    return Matrix{ 1 , bi.getW(), x };
}