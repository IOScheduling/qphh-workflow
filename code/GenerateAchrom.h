#pragma once


#ifndef FRAME_GENERATEACHROM_H
#include "common.h"
#define FRAME_GENERATEACHROM_H

//void UpdateITL(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& EndTime);
void UpdateITL(set<double>& ITLofRscId, double& StartTime, double& EndTime);
void W_Cal_Average_S(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b_S(vector<double>& RankList, vector<double>& ExeTime);
void Calculate_Rank_t_S(vector<double>& RankList, vector<double>& ExeTime);
void IntChr(chromosome& chrom);
void ModifyRscAlcLstByCode_RK(chromosome& chrom);
void GnrCode_RK(chromosome& chrom);
void GnrRscAlcTskSchLstFromCode_RK(chromosome& chrom);
vector<int> GnrSS_TS();
vector<int> GnrSS_Lvl();
chromosome GnrChr_HEFT_t_S(vector<double> rank_t);
chromosome GnrChr_HEFT_b_t_S(vector<double> rank_b_t);
chromosome GnrChr_HEFT_S(vector<double> Rank_b);
chromosome GnrChr_HEFT_b_ADBRKGA_S(vector<double> Rank_b);
chromosome GnrChr_HEFT_t_ADBRKGA_S(vector<double> Rank_t);
chromosome GnrChr_DIHEFT_S(vector<double>& Rank_b);
chromosome GnrChr_IHEFT3_b_S(vector<double> Rank_b);
chromosome GnrChr_IHEFT3_t_S(vector<double> Rank_t);
double IHEFT3_S(chromosome& ch);
double ClcAvrReadyTime_S(int TskId, chromosome& chrom);
vector<double> GnrDecimalsByAscend();
chromosome GnrChr_Lvl_Ran();
chromosome GnrChr_HMEC_S(vector<double> Rank_b);
chromosome GnrChr_HEFT_Baseline_S();
chromosome GnrPrtByRank_Rnd_S(vector<double>& Rank);
chromosome GnrPrtByRank_EFT_S(vector<double>& Rank);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR);
double GnrML_Evl_EFT_S(chromosome& ch);
double GnrML_Evl_MEC_S_QPHH(chromosome& ch);
double HrsDcd_EFT_S(chromosome& ch);
double NrmDcd_S(chromosome& ch);
double HrsDcd_CTP_S(chromosome& ch);
double GnrML_Evl_MEC_S(chromosome& ch);
double DcdEvl_S(chromosome& ch);
void AdpDcd_S(chromosome& chrom, double& CurTime, double& TotalTime);
void SelectRsc_EFT_S(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime);
double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
double CalculateEnergy2(chromosome& Chrom);
double CalculateEnergy1(chromosome& Chrom);
double CalculateECByDelta2(chromosome& Chrom, int& NumBer);
double CalculateECByDelta1(chromosome& Chrom, int& NumBer);
void RepairPriorityAndGnrSchOrd(chromosome& chrom);
void RepairMapAndGnrRscAlcLst(chromosome& ch);
int FindNearestRscId(int TaskId, double value);
void UpdateParticle(chromosome& ch, chromosome& Pbest, chromosome& Gbest, double& runtime, double& SchTime);
double IFBD_S(chromosome& ch);
double IFBS_S(chromosome& ch);
void LBCA_S(chromosome& chrom);
void LBCA_IFBS(chromosome& ch);
void InitProModelOfResAlc(vector<vector<double>>& PMR);
void InitProModelOfTskSch(vector<vector<double>>& PMS, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR);
chromosome GnrTskLstOfChr(vector<vector<double> >& PMS, vector<double>& eta_TSO);
void UpdatePMR(vector<vector<double>>& PMR, chromosome& bstChrom);
void UpdatePMS(vector<vector<double>>& PMS, chromosome& bstChrom);



chromosome GetTskSchLst(vector<double> Rank_b);
void Calculate_Rank_Rand_S(vector<double>& RankList, vector<double>& RandPri);
vector<int> GetOriginSequence(vector<double> seq);
#endif //FRAME_GENERATEACHROM_H
