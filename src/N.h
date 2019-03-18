#ifndef N_H
#define N_H

#include "NodeLocal.h"

class N
{
    double **matrixN;
    double **dNdKsi;
    double **dNdEta;
    public:
        N();
        virtual ~N();
        void calculateN();
        double **getMatrixN();
        double **getDnDKsi();
        double **getDnDEta();
        void printMatrixN();
        void printDnDKsi();
        void printDnDEta();

        double countN1(double, double);
        double countN2(double, double);
        double countN3(double, double);
        double countN4(double, double);
 private:
        double derivateOfShapeFunctionKsi1(double);
        double derivateOfShapeFunctionKsi2(double);
        double derivateOfShapeFunctionKsi3(double);
        double derivateOfShapeFunctionKsi4(double);

        double derivateOfShapeFunctionEta1(double);
        double derivateOfShapeFunctionEta2(double);
        double derivateOfShapeFunctionEta3(double);
        double derivateOfShapeFunctionEta4(double);
};

#endif // N_H
