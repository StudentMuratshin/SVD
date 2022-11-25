#include <iostream>
#include "function.h"
#include <vector>

using namespace std;
int main()
{
    Matrix mat = { 2,3,{
        3,2,2,
        2,3,-2
    } };
    
    Matrix sq = Eigen_Values(mat * mat.transpose());
    cout << sq;

}