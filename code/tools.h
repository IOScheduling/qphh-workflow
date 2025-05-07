#pragma once
//
// Created by qqq on 2021/9/15.
//

#ifndef FRAME_TOOLS_H
#include "common.h"
#define FRAME_TOOLS_H
void CalculateLevelList();
void CalculateDescendants();
void CalculateAncestors();
bool SortValueByDescend(pair<int, double>& a, pair<int, double>& b);
double CalculatePowerByLoad(double wl, int HTid);
bool SortByEnergyConsumption(chromosome& a, chromosome& b);
bool SortByEnergyConsumptionByDescend(chromosome& a, chromosome& b);

double RandomDouble(int start, int end);
double RandomDouble2(int start, int end);
void IndexSortByValueOnAscend(vector<int>& ind, vector<double>& value);
void IndexSortByValueOnAscend(vector<int>& ind, vector<int>& value);
void IndexSortByValueOnDescend(vector<int>& ind, vector<double>& value);


#endif //FRAME_TOOLS_H
