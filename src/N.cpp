#include "N.h"
#include <iostream>
#include <cmath>
N::N()
{
    int size = 4;
    matrixN = new double*[size]();
    dNdEta =new double*[size]();
    dNdKsi = new double*[size]();
    for(int i = 0; i < size; i++){
        matrixN[i]= new double[size]();
        dNdEta[i]= new double[size]();
        dNdKsi[i]= new double[size]();
    }
}
double **N::getDnDEta()
{
    return dNdEta;
}
double **N::getDnDKsi(){
    return dNdKsi;
}
double **N::getMatrixN()
{
    return matrixN;
}
N::~N()
{
    //dtor
}
void N::calculateN()
{
    double *eta = new double[4];
    double *ksi = new double[4];

    ksi[0] =-1/sqrt(3); eta[0]=-1/sqrt(3);

    ksi[1]=1/sqrt(3); eta[1]=-1/sqrt(3);

    ksi[2]=1/sqrt(3); eta[2]=1/sqrt(3);

    ksi[3]=-1/sqrt(3); eta[3]=1/sqrt(3);
    int size = 4;
    for(int i = 0; i < size; i++)
    {
        matrixN[i][0] = countN1(ksi[i],eta[i]);
        matrixN[i][1] = countN2(ksi[i],eta[i]);
        matrixN[i][2] = countN3(ksi[i],eta[i]);
        matrixN[i][3] = countN4(ksi[i],eta[i]);
    }
    for(int i = 0; i < size; i++)
    {
        dNdEta[i][0] = derivateOfShapeFunctionEta1(ksi[i]);
        dNdEta[i][1] = derivateOfShapeFunctionEta2(ksi[i]);
        dNdEta[i][2] = derivateOfShapeFunctionEta3(ksi[i]);
        dNdEta[i][3] = derivateOfShapeFunctionEta4(ksi[i]);

        dNdKsi[i][0] = derivateOfShapeFunctionKsi1(eta[i]);
        dNdKsi[i][1] = derivateOfShapeFunctionKsi2(eta[i]);
        dNdKsi[i][2] = derivateOfShapeFunctionKsi3(eta[i]);
        dNdKsi[i][3] = derivateOfShapeFunctionKsi4(eta[i]);
    }
}
double N::countN1(double ksi, double eta)
{
    return 0.25*(1-ksi)*(1-eta);
}
double N::countN2(double ksi, double eta)
{
    return 0.25*(1+ksi)*(1-eta);
}
double N::countN3(double ksi, double eta)
{
    return 0.25*(1+ksi)*(1+eta);
}
double N::countN4(double ksi, double eta)
{
    return 0.25*(1-ksi)*(1+eta);
}
double N::derivateOfShapeFunctionKsi1(double eta)
{
    return -0.25*(1-eta);
}
double N::derivateOfShapeFunctionKsi2(double eta)
{
    return 0.25*(1-eta);
}
double N::derivateOfShapeFunctionKsi3(double eta)
{
    return 0.25*(1+eta);
}
double N::derivateOfShapeFunctionKsi4(double eta)
{
    return -0.25*(1+eta);
}
double N::derivateOfShapeFunctionEta1(double ksi)
{
    return -0.25*(1-ksi);
}
double N::derivateOfShapeFunctionEta2(double ksi)
{
    return -0.25*(1+ksi);
}
double N::derivateOfShapeFunctionEta3(double ksi)
{
    return 0.25*(1+ksi);
}
double N::derivateOfShapeFunctionEta4(double ksi)
{
    return 0.25*(1-ksi);
}
void N::printMatrixN()
{
    int size = 4;
    for(int i = 0; i < size;i++){
        for(int j = 0; j < size; j++){
                std::cout<<matrixN[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
}
void N::printDnDEta()
{
    int size = 4;
    for(int i = 0; i < size;i++){
        for(int j = 0; j < size; j++){
                std::cout<<dNdEta[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
}
void N::printDnDKsi()
{
    int size = 4;

    for(int i = 0; i < size;i++){
        for(int j = 0; j < size; j++){
                std::cout<<dNdKsi[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
}
