#include "InputHelper.h"
#include <iostream>
InputHelper::InputHelper()
{
	fstream file;
	file.open("input.txt", ios::in);
	if (file.good())
	{
		string line;
        //wysokosc siatki
		getline(file, line, '=');
		getline(file, line, '\n');
		gridHeight = stod(line);

        //szerokosc siatki
		getline(file, line, '=');
		getline(file, line, '\n');
		gridWidth = stod(line);

        //liczba nodow w pionie
		getline(file, line, '=');
		getline(file, line, '\n');
		verticalNodeNumber = stoi(line);

        //liczba nodow w poziomie
		getline(file, line, '=');
		getline(file, line, '\n');
		horizontalNodeNumber = stoi(line);

        //wspolczynnik przewodzenia
		getline(file, line, '=');
		getline(file, line, '\n');
		thermalConductivity = stod(line);

        //temperatura poczatkowa
		getline(file, line, '=');
		getline(file, line, '\n');
		temperature = stod(line);

		//czas symulacji
        getline(file, line, '=');
		getline(file, line, '\n');
		simulationTime = stod(line);

		//przeskok
		getline(file, line, '=');
		getline(file, line, '\n');
		simulationStepTime = stod(line);

		//temperatura otoczenia
		getline(file, line, '=');
		getline(file, line, '\n');
		ambientTemperature = stod(line);


		getline(file, line, '=');
		getline(file, line, '\n');
		alfa = stod(line);

		//pojemnosc cieplna
		getline(file, line, '=');
		getline(file, line, '\n');
		specificHeat = stod(line);

		//gestosc
		getline(file, line, '=');
		getline(file, line, '\n');
		density = stod(line);

		file.close();
	}
	else{cout<<"file bad"<<endl;}
}

InputHelper::~InputHelper(){}

double InputHelper::getGridHeight()
{
	return gridHeight;
}
double InputHelper::getGridWidth()
{
    return gridWidth;
}
int InputHelper::getVerticalNodeNumber()
{
	return verticalNodeNumber;
}

int InputHelper::getHorizontalNodeNumber()
{
	return horizontalNodeNumber;
}

double InputHelper::getTemperature()
{
	return temperature;
}

double InputHelper::getThermalConductivity()
{
	return thermalConductivity;
}
double InputHelper::getSimulationTime()
{
    return simulationTime;
}
double InputHelper::getSimulationStepTime()
{
    return simulationStepTime;
}
double InputHelper::getAmbientTemperature()
{
    return ambientTemperature;
}
double InputHelper::getAlfa()
{
    return alfa;
}
double InputHelper::getSpecificHeat()
{
    return specificHeat;
}
double InputHelper::getDensity()
{
    return density;
}

double InputHelper::gridHeight;
double InputHelper::gridWidth;
int InputHelper::verticalNodeNumber;
int InputHelper::horizontalNodeNumber;
double InputHelper::thermalConductivity;
double InputHelper::temperature;
double InputHelper::simulationTime;
double InputHelper::simulationStepTime;
double InputHelper::ambientTemperature;
double InputHelper::alfa;
double InputHelper::specificHeat;
double InputHelper::density;
