#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"
#include "GenOperator.h"
#include <fstream>

chromosome runHGA(string Model, string ECfileName, string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);                  //read model information
    ConfigParameter_HGA();                          //set the parameter values
    CalculateLevelList();                           //calculate the levels of tasks
    //{calcualte the rank_b of tasks}
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    //    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average_S(ww);                             //calculate the average execution time of tasks
    //    C_Cal_Average(cc);                             //calculate the average transfer time among tasks
    Calculate_Rank_b_S(Rank_b, ww);    //calcualte the rank_b
    vector<chromosome> population(Parameter_HGA.NumOfChromPerPop);
    //    #pragma omp parallel for
    for (int n = 0; n < Parameter_HGA.NumOfChromPerPop - 1; ++n) {
        chromosome chrom;
        IntChr(chrom);
        for (int j = 0; j < comConst.NumOfTsk; ++j) {
            chrom.RscAlcLst[j] = Tasks[j].ElgRsc[rand() % Tasks[j].ElgRsc.size()];
        }
        GnrTskSchLst_HGA_S(chrom);
        DcdEvl_S(chrom);
        CalculateEnergy1(chrom);
        population[n] = chrom;
    }
    population[Parameter_HGA.NumOfChromPerPop - 1] = GnrChr_HEFT_S(Rank_b);
    sort(population.begin(), population.end(), SortByEnergyConsumption); //sorting
    //{Ensure that the elite are even numbers}
    int NumOfElite = int(Parameter_HGA.NumOfChromPerPop * Parameter_HGA.EliteRate);
    if (NumOfElite % 2 == 1) {
        ++NumOfElite;
    }
    while (1) {
        double originBstEC = population[0].EnergyConsumption;

        ++iteration;
        double RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (RunTime >= SchTime) {
            SchTime = RunTime;
            break;
        }
        vector<chromosome> NewPopulation(Parameter_HGA.NumOfChromPerPop);
        //{Elitism: Copy elite to new population}
#pragma omp parallel for
        for (int n = 0; n < NumOfElite; ++n) {
            NewPopulation[n] = population[n];
        }
        //{selection, crossover, and mutation}
#pragma omp parallel for
        for (int n = NumOfElite; n < Parameter_HGA.NumOfChromPerPop; n += 2) {
            int parent1 = -1;
            int parent2 = -1;
            SelectionTournament(parent1, parent2, Parameter_HGA.NumOfChromPerPop);
            chromosome TemChromosome1 = population[parent1];
            chromosome TemChromosome2 = population[parent2];
            Crossover_HGA(TemChromosome1, TemChromosome2);
            if (RandomDouble(0, 1) < Parameter_HGA.MutationRate) {
                Mutation_HGA(TemChromosome1);
            }
            if (RandomDouble(0, 1) < Parameter_HGA.MutationRate) {
                Mutation_HGA(TemChromosome2);
            }
            NewPopulation[n] = TemChromosome1;
            NewPopulation[n + 1] = TemChromosome2;
        }
        RscLoadAdjust_HGA(NewPopulation);
        sort(NewPopulation.begin(), NewPopulation.end(), SortByEnergyConsumption);
        population = NewPopulation;

        std::ofstream outFile(ECfileName, std::ios::app);
        if (!outFile) {
            std::cerr << "无法打开文件: " << ECfileName << std::endl;
            break;
        }

        string algName = "HGA";
        chromosome BstChrom = population[0];
        if (BstChrom.EnergyConsumption + PrecisionValue <= originBstEC)
        {
            outFile << algName << "\t" << Model << "\t" << comConst.NumOfTsk << "\t" << iteration << "\t" << RunTime << "\t"
                << BstChrom.EnergyConsumption << endl;
        }

    }
    return population[0];
}
