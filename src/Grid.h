#pragma once
#include "Node.h"
#include "Element.h"
#include <vector>
#include "InputHelper.h"
#include "N.h"
using namespace std;
class Grid
{
	vector<Node*> nodes;
	vector<Element*> elements;
	double **globalH;
	double **globalC;
	double *globalP;
	double time;
	double dTau;
public:
	Grid();
	~Grid();
	void showGrid();
	void calculate();
	Node** getByIndex(int *index);
	void makeGlobalFromLocal();
	void generateRest();
	bool gauss(int, double **,double*);
	double returnMin(int, double*);
	double returnMax(int, double*);
};
