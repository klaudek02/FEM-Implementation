#pragma once
#include "Node.h"
#include "InputHelper.h"
#include "Jakobian.h"
#include "MatrixH.h"
#include "MatrixC.h"
#include <vector>
class Element
{
    int id;
	int *nodesId;
	double thermalConductivity;
	double alfa;
	double c;
	double ro;
	double ambientTemperature;
	N n;
	Jakobian jakobian;
	MatrixH matrixH;
	MatrixC matrixC;
	double *vectorP;
	double *xP;
	double *yP;
public:
	Element(int id, int *nodesId);
	void calculate(Node**);
	~Element();
    int getId();
	int *getNodesId();
	Jakobian getJakobian();
	double **getMatrixH();
	double **getMatrixC();
	void calculateP(Node**);
	double *getVectorP();
	void calculateXpYp(Node**);
};
