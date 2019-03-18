#pragma once
#include <fstream>
#include <string>
using namespace std;
class InputHelper
{
	static double gridHeight;
	static double gridWidth;
	static int verticalNodeNumber;
	static int horizontalNodeNumber;
	static double thermalConductivity;
	static double temperature;
	static double simulationTime;
	static double simulationStepTime;
	static double ambientTemperature;
	static double alfa;
	static double specificHeat;
	static double density;
public:
	InputHelper();
	~InputHelper();
	static double getGridHeight();
	static double getGridWidth();
	static int getVerticalNodeNumber();
	static int getHorizontalNodeNumber();
	static double getTemperature();
	static double getThermalConductivity();
	static double getSimulationTime();
	static double getSimulationStepTime();
	static double getAmbientTemperature();
	static double getAlfa();
	static double getSpecificHeat();
	static double getDensity();
};
