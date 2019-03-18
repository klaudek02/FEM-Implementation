#pragma once
#include "InputHelper.h"
class Node
{

public:
    double x;
    double y;
    double t;
    bool bc = false;
    Node(double x, double y, double t);
    ~Node();
};
