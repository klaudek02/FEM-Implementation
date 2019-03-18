#ifndef MATRIX_H
#define MATRIX_H

#include "Jakobian.h"
class MatrixH
{
    double **dNdX;
    double **dNdY;
    double ***dNdXdNdXT;
    double ***dNdYdNdYT;
    double ***matrixSumWithDetK;
    double **matrixH;
    double **matrixH_BC;
    public:
        MatrixH();
        virtual ~MatrixH();
        void calculateMatrixH(Jakobian, N, double);
        void calculateMatrixH_BC(Node**, double alfa);
        double **getMatrixH();
    protected:
    private:
};

#endif // MATRIX_H
