#include "MatrixC.h"
#include <iostream>
using namespace std;
MatrixC::MatrixC()
{
    int size = 4;
    croNN = new double**[size]();
    matrixC = new double*[size]();
    for(int i =0; i < size; i++)
    {
        croNN[i] = new double*[size]();
        matrixC[i] = new double[size]();
        for(int j = 0; j < size; j++)
            croNN[i][j] = new double[size]();
    }
}
MatrixC::~MatrixC()
{
    //dtor
}
void MatrixC::calculateMatrixC(N nObject, double c, double ro, double *detJakobian )
{
    int size = 4;
    double **matrixN = nObject.getMatrixN();
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            for(int l = 0; l < size; l++)
            {
                croNN[i][j][l] = matrixN[i][l] * matrixN[i][j]*c*ro*detJakobian[i];
            }
        }
    }
    double suma = 0;
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            for(int l = 0; l < size; l++)
            {
                suma +=croNN[l][j][i];
            }
            matrixC[i][j] = suma;
            suma = 0;
        }
    }
}
double** MatrixC::getMatrixC()
{
    return matrixC;
}
