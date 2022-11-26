#include <iostream>
#include "function.h"
#include <vector>

using namespace std;
int main()
{
    Matrix mat = { 3,3,{
        -6,5.5,-1,
        5.5,1,-2,
        -1,-2,-3
    } };
    pair <Matrix,Matrix> ei = Eigen_Values(mat);
    cout << ei.first << endl << endl << ei.second;
}