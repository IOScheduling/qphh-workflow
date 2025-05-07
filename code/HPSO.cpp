

#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"
#include "GenOperator.h"

chromosome runHPSO(string Model, string ECfileName, string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);              //read model information
    ConfigParameter_HPSO();
    CalculateLevelList();
    //{calcualte the rank_b of tasks}
    vector<double> ww(comConst.NumOfTsk, 0);
    //    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    W_Cal_Average_S(ww);                           //calculate the average execution time of tasks
    //    C_Cal_Average(cc);                           //calculate the average transfer time among tasks
    vector<double> Rank_t(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    Calculate_Rank_t_S(Rank_t, ww);            //calcualte the rank_t
    Calculate_Rank_b_S(Rank_b, ww);            //calcualte the rank_b
    vector<double> Rank_b1(comConst.NumOfTsk, 0);
    double MaxRank_b = 0;
    for (int i = 0; i < Rank_b.size(); ++i) {
        if (MaxRank_b + PrecisionValue < Rank_b[i])
            MaxRank_b = Rank_b[i];
    }
    for (int i = 0; i < Rank_b.size(); ++i) {
        Rank_b1[i] = MaxRank_b - Rank_b[i];
    }
    vector<chromosome> population;
    population.push_back(GnrPrtByRank_EFT_S(Rank_b1));
    population.push_back(GnrPrtByRank_EFT_S(Rank_t));

    for (int i = 1; i < Parameter_HPSO.NumOfChromPerPop; ++i) {
        population.push_back(GnrPrtByRank_Rnd_S(Rank_b1));
        population.push_back(GnrPrtByRank_Rnd_S(Rank_t));
    }
    sort(population.begin(), population.end(), SortByEnergyConsumption);
    vector<chromosome>::iterator iter = population.begin();
    advance(iter, Parameter_HPSO.NumOfChromPerPop);
    population.assign(population.begin(), iter);
    vector<chromosome> Pbest = population;
    chromosome Gbest = population[0];
    while (1) {
        ++iteration;

        double originBstEC = Gbest.EnergyConsumption;


        double RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (RunTime >= SchTime) {
            SchTime = RunTime;
            break;
        }
#pragma omp parallel for
        for (int i = 0; i < Parameter_HPSO.NumOfChromPerPop; ++i) {
            UpdateParticle(population[i], Pbest[i], Gbest, RunTime, SchTime);
            DcdEvl_S(population[i]);
            CalculateEnergy1(population[i]);
        }
        for (int i = 0; i < Parameter_HPSO.NumOfChromPerPop; ++i) {
            if (population[i].EnergyConsumption + PrecisionValue < Pbest[i].EnergyConsumption) {
                Pbest[i] = population[i];
                if (Pbest[i].EnergyConsumption + PrecisionValue < Gbest.EnergyConsumption)
                {
                    Gbest = Pbest[i];
                }
            }
        }

        std::ofstream outFile(ECfileName, std::ios::app);
        if (!outFile) {
            std::cerr << "无法打开文件: " << ECfileName << std::endl;
            break;
        }

        string algName = "HPSO";
        chromosome BstChrom = Gbest;
        if (BstChrom.EnergyConsumption + PrecisionValue <= originBstEC)
        {
            outFile << algName << "\t" << Model << "\t" << comConst.NumOfTsk << "\t" << iteration << "\t" << RunTime << "\t"
                << BstChrom.EnergyConsumption << endl;
        }

        
    }
    return Gbest;
}