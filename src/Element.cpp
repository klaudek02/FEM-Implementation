#include "Element.h"
#include <iostream>
#include <cmath>
Element::Element(int id, int * nodesId)
{
    this->id = id;
    this->nodesId = nodesId;
    this->thermalConductivity = InputHelper::getThermalConductivity();
    this->alfa = InputHelper::getAlfa();
    this->c = InputHelper::getSpecificHeat();
    this->ro = InputHelper::getDensity();
    this->ambientTemperature = InputHelper::getAmbientTemperature();
    vectorP = new double[4];
    for(int i = 0; i < 4; i++)
        vectorP[i] = 0;
}
Element::~Element()
{

}

int Element::getId()
{
    return id;
}

int *Element::getNodesId()
{
    return nodesId;
}
void Element::calculate(Node** node)
{
    n.calculateN();
    calculateXpYp(node);
    jakobian.calculateJakobian(node,n);
    matrixH.calculateMatrixH(jakobian,n,thermalConductivity);
    matrixH.calculateMatrixH_BC(node, alfa);
    matrixC.calculateMatrixC(n,c,ro,jakobian.getDetJakobian());
    calculateP(node);
}
void Element::calculateXpYp(Node** node)
{
    xP = new double[4];
    yP = new double[4];
    double **N = n.getMatrixN();
    for(int i = 0; i < 4; i++)
    {
        xP[i] = node[i]->x*N[0][i] + node[i]->x*N[1][i] + node[i]->x*N[2][i] + node[i]->x*N[3][i];
        yP[i] = node[i]->y*N[0][i] + node[i]->y*N[1][i] + node[i]->y*N[2][i] + node[i]->y*N[3][i];
    }
}
Jakobian Element::getJakobian()
{
    return jakobian;
}
double **Element::getMatrixH()
{
    return matrixH.getMatrixH();
}
double **Element::getMatrixC()
{
    return matrixC.getMatrixC();
}
void Element::calculateP(Node** nodes)
{
    double sideLength;
    double N[2][4];
    double pkt[8][2] = {{-1/sqrt(3),-1},{1/sqrt(3),-1},{1,-1/sqrt(3)},{1,1/sqrt(3)},
        {1/sqrt(3),1},{-1/sqrt(3),1},{-1,1/sqrt(3)},{-1,-1/sqrt(3)}};
    for (int l = 0; l < 4; l++)
    {
        int firstNode = l;
        int secondNode = l + 1;
        if (secondNode >= 4) secondNode = 0;

        sideLength=sqrt(pow(nodes[firstNode]->x-nodes[secondNode]->x,2) +
                         pow(nodes[firstNode]->y-nodes[secondNode]->y,2));
        if(nodes[firstNode]->bc && nodes[secondNode]->bc)
        {
            N[0][0] = 0.25*(1-pkt[l*2][0])*(1-pkt[l*2][1]);
            N[0][1] = 0.25*(1+pkt[l*2][0])*(1-pkt[l*2][1]);
            N[0][2] = 0.25*(1+pkt[l*2][0])*(1+pkt[l*2][1]);
            N[0][3] = 0.25*(1-pkt[l*2][0])*(1+pkt[l*2][1]);

            N[1][0] = 0.25*(1-pkt[l*2+1][0])*(1-pkt[l*2+1][1]);
            N[1][1] = 0.25*(1+pkt[l*2+1][0])*(1-pkt[l*2+1][1]);
            N[1][2] = 0.25*(1+pkt[l*2+1][0])*(1+pkt[l*2+1][1]);
            N[1][3] = 0.25*(1-pkt[l*2+1][0])*(1+pkt[l*2+1][1]);

            for (int i = 0; i< 2; i++)
                for (int j = 0; j < 4; j++)
                    vectorP[j]+=(N[i][j]*alfa*ambientTemperature*sideLength/2);
        }
    }
}
double* Element::getVectorP()
{
    return vectorP;
}


