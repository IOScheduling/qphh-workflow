#pragma once

#include <fstream>
#ifndef FRAME_CONFIG_H
#include "common.h"
#define FRAME_CONFIG_H
#define ITERATIONUM 100

int readID(string id);
void ReadFile(string filename, string RscAlcFile);
void TaskCombine();
void DeleteFirstLineInFile(string fileName);
void ConfigParameter_HGA();
void ConfigParameter_NGA();
void ConfigParameter_LWSGA();
void ConfigParameter_CGA();
void ConfigParameter_HPSO();
void ConfigParameter_TSEDA();
void ConfigParameter_ADBRKGA();
void ConfigParameter_QPHH();
void ClearALL();

#endif //FRAME_CONFIG_H
