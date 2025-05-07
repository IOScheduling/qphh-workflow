#pragma once


#ifndef FRAME_HGA_H
#include "common.h"
#define FRAME_HGA_H


chromosome runHGA(string Model, string ECfileName, string XmlFile, string RscAlcFile, double& SchTime, int& iteration);


#endif //FRAME_HGA_H
