#ifndef MATRIXC_H
#define MATRIXC_H

#include "N.h"
class MatrixC
{
    double ***croNN;
    double **matrixC;
    public:
        MatrixC();
        virtual ~MatrixC();
        void calculateMatrixC(N, double, double, double*);
        double **getMatrixC();
};

#endif // MATRIXC_H
