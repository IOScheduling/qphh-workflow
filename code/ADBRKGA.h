#pragma once


#ifndef FRAME_ADBRKGA_H
#include "common.h"
#define FRAME_ADBRKGA_H

chromosome runADBRKGA(string Model, string ECfileName, string XmlFile, string RscAlcFile, double& SchTime, int& iteration);

#endif //FRAME_ADBRKGA_H
