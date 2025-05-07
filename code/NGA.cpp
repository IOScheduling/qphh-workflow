

#include "NGA.h"
#include "common.h"
#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"
#include "GenOperator.h"

chromosome runNGA(string Model, string ECfileName, string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);       //read model information
    ConfigParameter_NGA();                       //set the parameter values
    CalculateLevelList();                        //calculate the levels of tasks
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> Rank_t(comConst.NumOfTsk, 0);
    vector<double> Rank_b_t(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    //    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average_S(ww);
    //    C_Cal_Average(cc);
    Calculate_Rank_b_S(Rank_b, ww);
    Calculate_Rank_t_S(Rank_t, ww);

    for (int i = 0; i < comConst.NumOfTsk; ++i)
        Rank_b_t[i] = Rank_b[i] + Rank_t[i];
    //{initialize the population}
    vector<chromosome> population(Parameter_NGA.NumOfChromPerPop);
    //seed  HEFT_b_t, HEFT_t, HEFT_b into population;
    population[0] = GnrChr_HEFT_b_t_S(Rank_b_t);
    population[1] = GnrChr_HEFT_t_S(Rank_t);
    population[2] = GnrChr_HEFT_S(Rank_b);
#pragma omp parallel for
    for (int n = 3; n < Parameter_NGA.NumOfChromPerPop; ++n) {
        chromosome TemChrom;
        IntChr(TemChrom);
        TemChrom.TskSchLst = GnrSS_TS();
        GnrML_Evl_EFT_S(TemChrom);
        CalculateEnergy1(TemChrom);
        //        GnrML_Evl_MEC_S(TemChrom);
        population[n] = TemChrom;
    }
    sort(population.begin(), population.end(), SortByEnergyConsumption);

    while (1) {
        ++iteration;
        double originBstEC = population[0].EnergyConsumption;


        double RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (RunTime >= SchTime) {
            SchTime = RunTime;
            break;
        }
        int NumOfElites = int(Parameter_NGA.NumOfChromPerPop * Parameter_NGA.EliteRate);
        vector<chromosome> NextPopulation(Parameter_NGA.NumOfChromPerPop);
        //{Copy elite to new population}
#pragma omp parallel for
        for (int n = 0; n < NumOfElites; ++n) {
            NextPopulation[n] = population[n];
        }
        //{crossover}
        bool flag = true;
#pragma omp parallel for
        for (int n = NumOfElites; n < Parameter_NGA.NumOfChromPerPop; ++n) {
            int Pop1 = rand() % Parameter_NGA.NumOfChromPerPop;
            int Pop2 = NumOfElites + rand() % (Parameter_NGA.NumOfChromPerPop - NumOfElites);
            while (Pop2 == Pop1) {
                Pop1 = rand() % Parameter_NGA.NumOfChromPerPop;
            }
            NextPopulation[n] = Crossover_NGA(population[Pop1], population[Pop2], flag);
            flag = !flag;
        }
        //{mutation}
#pragma omp parallel for
        for (int n = 0; n < Parameter_NGA.NumOfChromPerPop; ++n) {
            if (RandomDouble(0, 1) < Parameter_NGA.MutationRate) {
                Mutation_NGA(NextPopulation[n]);
            }
        }
#pragma omp parallel for
        for (int n = 0; n < Parameter_NGA.NumOfChromPerPop; ++n) {
            GnrML_Evl_EFT_S(NextPopulation[n]);
            CalculateEnergy1(NextPopulation[n]);
            //            GnrML_Evl_MEC_S(NextPopulation[n]); //if using it, the ch.HtUseTime must be initialized as {9999999999,-1}
        }
        sort(NextPopulation.begin(), NextPopulation.end(), SortByEnergyConsumption);
        population = NextPopulation;

        std::ofstream outFile(ECfileName, std::ios::app);
        if (!outFile) {
            std::cerr << "无法打开文件: " << ECfileName << std::endl;
            break;
        }

        string algName = "NGA";
        chromosome BstChrom = population[0];
        if (BstChrom.EnergyConsumption + PrecisionValue <= originBstEC)
        {
            outFile << algName << "\t" << Model << "\t" << comConst.NumOfTsk << "\t" << iteration << "\t" << RunTime << "\t"
                << BstChrom.EnergyConsumption << endl;
        }
    }
    return population[0];
}

