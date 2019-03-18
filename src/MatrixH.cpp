#include "MatrixH.h"
#include "PrintUtil.cpp"
using namespace std;
#include <iostream>
#include <math.h>
MatrixH::MatrixH()
{
    int size = 4;
    dNdX = new double*[size]();
    dNdY = new double*[size]();
    dNdXdNdXT = new double**[size]();
    dNdYdNdYT = new double**[size]();
    matrixSumWithDetK = new double**[size]();
    matrixH = new double*[size]();
    matrixH_BC = new double*[size];
    for(int i = 0; i < size; i ++)
    {
        dNdX[i] = new double[size];
        dNdY[i] = new double[size];
        matrixH[i] = new double[size];
        matrixH_BC[i] = new double[size];
        dNdYdNdYT[i] = new double*[size];
        dNdXdNdXT[i] = new double*[size];
        matrixSumWithDetK[i] = new double*[size];
        for (int j=0; j<size; j++)
        {
            dNdYdNdYT[i][j]=new double[size];
            dNdXdNdXT[i][j]=new double[size];
            matrixSumWithDetK[i][j]=new double[size];
            matrixH_BC[i][j] = 0;
        }
    }
}
MatrixH::~MatrixH()
{
    //dtor
}
void MatrixH::calculateMatrixH(Jakobian jakobian, N nObject, double k)
{
    int size = 4;
    double **invertibleMatrixJakobian = jakobian.getInvertibleMatrixJakobian();
    double *detJakobian = jakobian.getDetJakobian();
    double **dNdKsi = nObject.getDnDKsi();
    double **dNdEta = nObject.getDnDEta();
    for(int i =0; i < size; i++)
    {
        dNdY[i][0]=invertibleMatrixJakobian[i][2]*dNdKsi[i][0]+invertibleMatrixJakobian[i][3]*dNdEta[i][0];
        dNdY[i][1]=invertibleMatrixJakobian[i][2]*dNdKsi[i][1]+invertibleMatrixJakobian[i][3]*dNdEta[i][1];
        dNdY[i][2]=invertibleMatrixJakobian[i][2]*dNdKsi[i][2]+invertibleMatrixJakobian[i][3]*dNdEta[i][2];
        dNdY[i][3]=invertibleMatrixJakobian[i][2]*dNdKsi[i][3]+invertibleMatrixJakobian[i][3]*dNdEta[i][3];

        dNdX[i][0]=invertibleMatrixJakobian[i][0]*dNdKsi[i][0]+invertibleMatrixJakobian[i][1]*dNdEta[i][0];
        dNdX[i][1]=invertibleMatrixJakobian[i][0]*dNdKsi[i][1]+invertibleMatrixJakobian[i][1]*dNdEta[i][1];
        dNdX[i][2]=invertibleMatrixJakobian[i][0]*dNdKsi[i][2]+invertibleMatrixJakobian[i][1]*dNdEta[i][2];
        dNdX[i][3]=invertibleMatrixJakobian[i][0]*dNdKsi[i][3]+invertibleMatrixJakobian[i][1]*dNdEta[i][3];

    }
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            for (int b=0; b<size; b++)
            {
                dNdXdNdXT[i][j][b]=dNdX[i][b]*dNdX[i][j];
                dNdYdNdYT[i][j][b]=dNdY[i][b]*dNdY[i][j];
            }
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            for (int b=0; b<size; b++)
                matrixSumWithDetK[i][j][b] = (dNdXdNdXT[i][j][b]+dNdYdNdYT[i][j][b])*k*detJakobian[i];
    for(int i =0; i < size; i++)
        for(int j = 0; j< size; j++)
            for(int l = 0; l < size; l++)
                matrixH[i][j] += matrixSumWithDetK[l][i][j];
}
void MatrixH::calculateMatrixH_BC(Node** nodes, double alfa)
{
    double sideLength;
    double N[2][4], pc0[4][4], pc1[4][4];
    double pkt[8][2] = {{-1/sqrt(3),-1},{1/sqrt(3),-1},{1,-1/sqrt(3)},{1,1/sqrt(3)},
        {1/sqrt(3),1},{-1/sqrt(3),1},{-1,1/sqrt(3)},{-1,-1/sqrt(3)}};
    for(int l = 0; l < 4; l++)
    {
        int firstNode = l;
        int secondNode = l + 1;
        if(secondNode>=4) secondNode = 0;
        sideLength = sqrt(pow(nodes[firstNode]->x-nodes[secondNode]->x,2) +
                         pow(nodes[firstNode]->y-nodes[secondNode]->y,2));
        if ((nodes[firstNode]->bc && nodes[secondNode]->bc))
        {
            N[0][0]=0.25*(1-pkt[l*2][0])*(1-pkt[l*2][1]);
            N[0][1]=0.25*(1+pkt[l*2][0])*(1-pkt[l*2][1]);
            N[0][2]=0.25*(1+pkt[l*2][0])*(1+pkt[l*2][1]);
            N[0][3]=0.25*(1-pkt[l*2][0])*(1+pkt[l*2][1]);

            N[1][0]=0.25*(1-pkt[l*2+1][0])*(1-pkt[l*2+1][1]);
            N[1][1]=0.25*(1+pkt[l*2+1][0])*(1-pkt[l*2+1][1]);
            N[1][2]=0.25*(1+pkt[l*2+1][0])*(1+pkt[l*2+1][1]);
            N[1][3]=0.25*(1-pkt[l*2+1][0])*(1+pkt[l*2+1][1]);

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    pc0[i][j] = N[0][j] * N[0][i] * alfa;

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    pc1[i][j] = N[1][j] * N[1][i] * alfa;

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    matrixH_BC[i][j] += (pc0[i][j] + pc1[i][j]) * (sideLength / 2);
        }
    }
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            matrixH[i][j] += matrixH_BC[i][j];
}
double **MatrixH::getMatrixH()
{
    return matrixH;
}

