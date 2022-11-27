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
        
        mat.SVD();


    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}