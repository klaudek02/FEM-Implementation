#ifndef JAKOBIAN_H
#define JAKOBIAN_H

#include "N.h"
#include <vector>
#include "Node.h"
class Jakobian
{
    public:
    double **matrixJakobian;
    double *detJakobian;
    double **invertibleMatrixJakobian;
        Jakobian();
        virtual ~Jakobian();
        double ** getInvertibleMatrixJakobian();
        double *getDetJakobian();
        void calculateJakobian(Node**, N);
        double calculateJakobianMatrixX(double *, Node**);
        double calculateJakobianMatrixY(double *, Node**);
        double calculateDetJ(double *);
        void printMatrixJakobian();
        void printDetJakobian();
        void printInvertibleMatrixJakobian();
    private:
};

#endif // JAKOBIAN_H
