#pragma once


#ifndef FRAME_GENOPERATOR_H
#include "common.h"
#define FRAME_GENOPERATOR_H

chromosome Crossover_NGA(chromosome& pop1, chromosome& pop2, bool& flag);
void CrsSS_TS_L(chromosome& ch1, chromosome& ch2, int CrossPoint);
void CrsSS_TS_R(chromosome& ch1, chromosome& ch2, int CrossPoint);
void Mutation_NGA(chromosome& chrom);
void GnrTskSchLst_HGA_S(chromosome& chrom);
void Crossover_HGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_HGA(chromosome& chrom);
void RscLoadAdjust_HGA(vector<chromosome>& chromosomes);
void SelectionTournament(int& parent_1, int& parent_2, int& NumOfChromPerPop);
void CrsMS_SP(chromosome& ch1, chromosome& ch2);
void CrsMS_DP(chromosome& ch1, chromosome& ch2);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A, int& NumOfChromPerPop);
int SelectChrom(vector<double>& A);
void MtnSS_TS(chromosome& a);
void MtnMS_MP(chromosome& ch);
void CrsMS_MP(chromosome& chrom1, chromosome& chrom2);
void Crossover_LWSGA(chromosome& ch1, chromosome& ch2);
void Crs_Lvl(chromosome& chrom1, chromosome& chrom2);
void CrsSS_ExcTskInLvl(chromosome& chrom1, chromosome& chrom2);
void Mutation_LWSGA(chromosome& ch);
void MtnSS_ExcTskInLvl(chromosome& chrom);
void Mtn_rebuild_level(chromosome& ch);
void Crossover_CGA(chromosome& pop1, chromosome& pop2);
void Mutation_CGA(chromosome& a);
void MtnMS_SP(chromosome& ch);
chromosome Crs_BPUC(chromosome& chrom1, chromosome& chrom2);
#endif //FRAME_GENOPERATOR_H
