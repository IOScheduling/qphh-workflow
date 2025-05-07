#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"
#include "GenOperator.h"

chromosome runLWSGA(string Model, string ECfileName, string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);                                      //read model information
    ConfigParameter_LWSGA();                                            //set the parameter values
    CalculateLevelList();                                               //calculate the levels of tasks
    vector<chromosome> population(Parameter_LWSGA.NumOfChromPerPop);
    //#pragma omp parallel for
    for (int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; ++n) {
        chromosome chrom;
        IntChr(chrom);
        chrom.TskSchLst = GnrSS_Lvl();
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            chrom.RscAlcLst[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()];
        }
        population[n] = chrom;
    }
    //#pragma omp parallel for
    for (int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; ++n) {
        DcdEvl_S(population[n]);                        //decoding
        CalculateEnergy1(population[n]);
    }
    sort(population.begin(), population.end(), SortByEnergyConsumption);//sorting
    while (1) {
        double originBstEC = population[0].EnergyConsumption;

        ++iteration;

        double RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (RunTime >= SchTime) {
            SchTime = RunTime;
            break;
        }
        vector<chromosome> NewPopulation(Parameter_LWSGA.NumOfChromPerPop);
#pragma omp parallel for
        for (int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; n += 2) {
            int parent1 = -1;
            int parent2 = -1;
            SelectionTournament(parent1, parent2, Parameter_LWSGA.NumOfChromPerPop); //select two chromosomes using the tournament method
            double rand = RandomDouble(0, 1);
            chromosome TemChromosome1 = population[parent1];
            chromosome TemChromosome2 = population[parent2];
            if (rand < Parameter_LWSGA.CrossoverRate) {
                Crossover_LWSGA(TemChromosome1, TemChromosome2); //crossover
            }
            else {
                Mutation_LWSGA(TemChromosome1);                      //mutation
                Mutation_LWSGA(TemChromosome2);
            }
            NewPopulation[n] = TemChromosome1;
            NewPopulation[n + 1] = TemChromosome2;
        }
#pragma omp parallel for
        for (int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; ++n) {
            DcdEvl_S(NewPopulation[n]);
            CalculateEnergy1(NewPopulation[n]);
        }
        //{generate the next population}
        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end());
        sort(population.begin(), population.end(), SortByEnergyConsumption);
        population.resize(Parameter_LWSGA.NumOfChromPerPop);

        std::ofstream outFile(ECfileName, std::ios::app);
        if (!outFile) {
            std::cerr << "无法打开文件: " << ECfileName << std::endl;
            break;
        }

        string algName = "LWSGA";
        chromosome BstChrom = population[0];
        if (BstChrom.EnergyConsumption + PrecisionValue <= originBstEC)
        {
            outFile << algName << "\t" << Model << "\t" << comConst.NumOfTsk << "\t" << iteration << "\t" << RunTime << "\t"
                << BstChrom.EnergyConsumption << endl;
        }
        
    }
    return population[0];
}
