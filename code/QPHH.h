#pragma once


#ifndef FRAME_QPHH_H
#include "common.h"
#define FRAME_QPHH_H
#define NUMSACTIONS 3
#define NUMSTATES 13
#define NUMRANDCHROM 50
#define NUMELITECHROM 20
#define MAXDOUBLE numeric_limits<double>::max()


chromosome runQPHH(string Model, string ECfileName,string XmlFile, string RscAlcFile, double& SchTime, int& iteration, double PopSizeRate, int selectCount,double epslion, double gamma);
void printVector(vector<int> vec);
void swapTwoTasksInSchedule(chromosome& Chrom, set<pair<int, int>>& swappedPairs);
chromosome insertTaskInSchedule(chromosome& Chrom);
chromosome insertTaskInScheduleFrontHal(chromosome& Chrom);
chromosome insertTaskInScheduleBackHal(chromosome& Chrom);
void insertTaskToBstPosition(chromosome& Chrom);
chromosome insertTaskAfterNearestParent(chromosome& Chrom);
void swapTwoTasksInSameLvl(chromosome& Chrom, set<pair<int, int>>& swappedPairs);
chromosome insertTaskPairsToBstPosition(chromosome& NewChrom);
#endif //FRAME_QPHH_H
#pragma once
