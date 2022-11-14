#include <iostream>
#include <Eigen/SVD>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

int main()
{
    MatrixXf m = MatrixXf::Random(3, 2);
    cout << "Here is the matrix m:" << endl << m << endl;
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
    Vector3f rhs(1, 0, 0);
    cout << "Now consider this rhs vector:" << endl << rhs << endl;
    cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;
    cout << svd.matrixU() * svd.singularValues() * svd.matrixV();
}