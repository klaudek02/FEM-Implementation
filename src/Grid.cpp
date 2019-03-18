#include "Grid.h"
#include <iostream>
#include <cmath>

Grid::Grid()
{
    int nH = InputHelper::getVerticalNodeNumber();
    int nL = InputHelper::getHorizontalNodeNumber();
    double t0 = InputHelper::getTemperature();
    double gridHeight = InputHelper::getGridHeight();
    double gridWidth = InputHelper::getGridWidth();
    double nodeHeight = gridHeight/ (nH-1);
    double nodeWidth = gridWidth/ (nL-1);

    for (int i = 0; i < nL; i++)
        for (int j = 0; j < nH; j++)
        {
            Node *node = new Node(i*nodeWidth, j*nodeHeight, t0);
            if ((i == 0 || i == nL - 1) || (j == 0 || j == nH - 1))
                node->bc = true;
            nodes.push_back(node);
        }
    for (int i = 0, id = 0; i < nL - 1; i++)
        for (int j = 1; j < nH; j++)
        {
            int * nodesId = new int[4];
            nodesId[0] = j +  nH * i;
            nodesId[1] = j +  nH  * (i+1);
            nodesId[2] = j + 1 + nH * (i+1);
            nodesId[3] = j + 1 + nH * i;
            elements.push_back(new Element(++id, nodesId));
        }
    int size = nodes.size();
    globalH = new double*[size];
    globalC = new double*[size];
    globalP = new double[size];
    for(int i = 0; i < size; i++)
    {
        globalP[i] = 0;
        globalH[i] = new double[size];
        globalC[i] = new double[size];
        for(int j = 0; j < size; j++)
        {
            globalH[i][j] = 0;
            globalC[i][j] = 0;
        }
    }
    time = InputHelper::getSimulationTime();
    dTau = InputHelper::getSimulationStepTime();
}

void Grid::calculate()
{
    for(int i =0; i < elements.size(); i++)
    {
        Node** nodesArray=getByIndex(elements[i]->getNodesId());
        elements[i]->calculate(nodesArray);
    }
    makeGlobalFromLocal();
    generateRest();

}
Grid::~Grid()
{
}
void Grid::showGrid()
{
    for (int i = 0; i < elements.size(); i++)
    {
        int *nodesId = elements[i]->getNodesId();
        cout<<"["<<elements[i]->getId()<<"]: "<<"\t [";
        for(int i = 0; i < 4; i++)
            cout<<nodesId[i]<<" ";
        cout<<"] "<<endl;

    }
}
Node** Grid::getByIndex(int *index)
{
    int size = 4;
    Node **tempNode = new Node*[size];
    Node **returnArray = new Node*[size];
    for(int i = 0; i < size; i++)
        returnArray[i] = nodes[index[i]-1];
    return returnArray;
}
void Grid::makeGlobalFromLocal()
{
    for(int i = 0; i < elements.size(); i++)
    {
        double **localH = elements[i]->getMatrixH();
        double **localC = elements[i]->getMatrixC();
        double *localP = elements[i]->getVectorP();
        int *nodesId = elements[i]->getNodesId();
        for(int j = 0; j < 4; j++)
        {
            globalP[nodesId[j]-1] += localP[j];
            for(int k = 0; k < 4; k++)
            {
                globalH[nodesId[j]-1][nodesId[k]-1] +=  localH[j][k];
                globalC[nodesId[j]-1][nodesId[k]-1] +=  localC[j][k];
            }
        }
    }
}
void Grid::generateRest()
{
    int size = nodes.size();

    double **CdT = new double*[size];
    double *CdTT0P = new double[size];
    double **HCdT = new double*[size];
    double **Equations = new double*[size];

    for(int i = 0; i < size; i++)
    {
        CdT[i] = new double[size];
        HCdT[i] = new double[size];
        Equations[i] = new double[size+1];
    }
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            CdT[i][j] = globalC[i][j]/dTau;
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            HCdT[i][j] = globalH[i][j] + CdT[i][j];

    double *t1;
    t1 = new double[size];
    int iterations = time/dTau;
    for(int x = 0; x < iterations; x++)
    {
        for(int i = 0; i<size; i++ )
        {
            CdTT0P[i] = 0.0;
            for(int j = 0; j < size; j++)
                CdTT0P[i] += (CdT[i][j]*nodes[j]->t);
        }
        for(int i = 0; i < size; i++)
            CdTT0P[i] += globalP[i];
        for(int i = 0; i < size; i++)
        {
            for(int j = 0; j < size+1; j++)
            {
                if(j == size)
                    Equations[i][j] = CdTT0P[i];
                else
                    Equations[i][j] = HCdT[i][j];
            }
        }
        for(int i = 0; i < size; i++)
            t1[i] = 0.0;
        gauss(size, Equations, t1);
        double minT = returnMin(size, t1);
        double maxT = returnMax(size, t1);
        for(int i = 0; i < size; i++)
            nodes[i]->t = t1[i];
        cout<<"Min: "<<minT<<"Max: "<<maxT<<endl;
        if(x == 19)
            break;
    }
}
double Grid::returnMin(int n, double *tablica)
{
    double min = tablica[0];
    for(int i = 1; i <n; i++)
    {
        if(min>tablica[i])
            min = tablica[i];
    }
    return min;
}
double Grid::returnMax(int n, double *tablica)
{
    double max = tablica[0];
    for(int i = 0; i <n; i++)
    {
        if(max<tablica[i])
            max=tablica[i];
    }
    return max;
}
bool Grid::gauss(int n, double ** AB, double *X)
{

    int i,j,k;
    double m,s;
    double eps = 1e-12;
    // eliminacja wspólczynników
    for(i = 0; i < n - 1; i++)
    {
        for(j = i + 1; j < n; j++)
        {
            if(fabs(AB[i][i]) < eps)
                return false;
            m = -AB[j][i] / AB[i][i];
            for(k = i + 1; k <= n; k++)
                AB[j][k] += m * AB[i][k];
        }
    }

    // wyliczanie niewiadomych
    for(i = n - 1; i >= 0; i--)
    {
        s = AB[i][n];
        for(j = n - 1; j >= i + 1; j--)
            s -= AB[i][j] * X[j];
        if(fabs(AB[i][i]) < eps)
            return false;
        X[i] = s / AB[i][i];
    }
    return true;
}
