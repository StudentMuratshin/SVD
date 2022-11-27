#pragma once
#include <vector>
#include <iostream>
using namespace std;

class Matrix {
	struct MatrixSize {
		int w;
		int h;
	} size;
	std::vector<double> vec;
public:

	Matrix(int h, int w, const double* coefficients);
	Matrix(int h, int w);
	Matrix(int h, int w, std::initializer_list<double> coefficients);
	static Matrix I(int s, double a = 1.0);

	double get(int i, int j)const;
	double& ref(int i, int j);

	Matrix operator*(const Matrix& b) const;

	bool operator==(const Matrix& b) const;

	Matrix operator*(const double b) const;

	friend ostream& operator<<(ostream& out, Matrix& x);

	Matrix operator/(const double b) const;

	Matrix operator+(const Matrix& other)const;

	Matrix operator-(const Matrix& other)const;

	Matrix getColumn(int col);

	Matrix getDiagonal();

	Matrix getRow(int col);

	void SVD();

	double getMaxVal();

	pair<Matrix, Matrix> QR();

	Matrix transpose()const;
	int getW() const;
	int getH() const;

};

Matrix Make_square(const Matrix b);

Matrix Gram_Schmidt(Matrix& arr);

double Dot_product(const Matrix a, const Matrix b);

pair<Matrix, Matrix> Eigen_Values_Vectors_symmetrical(const Matrix B);

Matrix Eigen_Values(const Matrix B);
