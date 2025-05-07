

#include "common.h"
#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"
#include "GenOperator.h"
chromosome runHEFT(string ECfileName, string XmlFile, string RscAlcFile, double& SchTime) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    CalculateLevelList();
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    //    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average_S(ww);
    //    C_Cal_Average(cc);
    Calculate_Rank_b_S(Rank_b, ww);
    chromosome Chrom_HEFT_b = GnrChr_HEFT_S(Rank_b);
    SchTime = (double)(clock() - start) / CLOCKS_PER_SEC;
    return Chrom_HEFT_b;
}