#include "Jakobian.h"
#include <iostream>
Jakobian::Jakobian()
{
    int size = 4;
    matrixJakobian = new double*[size];
    invertibleMatrixJakobian =new double*[size];
    detJakobian = new double[size];
    for(int i = 0; i < size; i++){
        matrixJakobian[i]= new double[size];
        invertibleMatrixJakobian[i]= new double[size];
    }
}

Jakobian::~Jakobian()
{
    //dtor
}

void Jakobian::calculateJakobian(Node**nodes,N nObject)
{
    int size = 4;
    double **dNdKsi = nObject.getDnDKsi();
    double **dNdEta = nObject.getDnDEta();
    for(int i = 0; i < size; i++)
    {
        matrixJakobian[i][0] = calculateJakobianMatrixX(dNdKsi[0],nodes);
        matrixJakobian[i][1] = calculateJakobianMatrixY(dNdKsi[1],nodes);
        matrixJakobian[i][2] = calculateJakobianMatrixX(dNdEta[2],nodes);
        matrixJakobian[i][3] = calculateJakobianMatrixY(dNdEta[3],nodes);
    }
    for(int i =0; i < size; i++)
        detJakobian[i] = calculateDetJ(matrixJakobian[i]);
    for(int i = 0; i < size; i++)
    {
        invertibleMatrixJakobian[i][0] = matrixJakobian[i][3]/detJakobian[i];
        invertibleMatrixJakobian[i][1] = -matrixJakobian[i][1]/detJakobian[i];
        invertibleMatrixJakobian[i][2] = -matrixJakobian[i][2]/detJakobian[i];
        invertibleMatrixJakobian[i][3] = matrixJakobian[i][0]/detJakobian[i];
    }
}

double Jakobian::calculateJakobianMatrixX(double *a, Node**node)
{
    double suma = 0;
    for(int i = 0 ; i < 4; i++)
    {
        suma+=a[i]*node[i]->x;
    }
    return suma;
}
double Jakobian::calculateJakobianMatrixY(double *a, Node**node)
{
    double suma = 0;
    for(int i = 0 ; i < 4; i++)
    {
        suma+=a[i]*node[i]->y;
    }
    return suma;
}
double Jakobian::calculateDetJ(double *a)
{
    return a[0]*a[3]-a[1]*a[2];
}
void Jakobian::printDetJakobian()
{
    int size = 4;
    for(int i = 0; i < size; i++)
        std::cout<<detJakobian[i]<<" ";
    std::cout<<endl<<endl;
}
void Jakobian::printMatrixJakobian()
{
    int size = 4;
    for(int i = 0; i < size; i++)
    {
        for(int j = 0 ;j < size; j++)
            std::cout<<matrixJakobian[i][j]<<" ";
        std::cout<<endl;
    }
    std::cout<<endl<<endl;
}
void Jakobian::printInvertibleMatrixJakobian()
{
    int size = 4;
    for(int i = 0; i < size; i++)
    {
        for(int j = 0 ; j < size; j++)
            std::cout<<invertibleMatrixJakobian[i][j]<<" ";
        std::cout<<endl;
    }
    std::cout<<endl<<endl;
}
double ** Jakobian::getInvertibleMatrixJakobian()
{
    return invertibleMatrixJakobian;
}
double *Jakobian::getDetJakobian()
{
    return detJakobian;
}
