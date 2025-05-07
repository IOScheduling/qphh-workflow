
#include "TSEDA.h"
#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"

chromosome runTSEDA(string Model, string ECfileName, string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    string ECfileName2 = "../TSEDA_ITEC2.txt";
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_TSEDA();//配置变量
    CalculateLevelList();//计算层级列表
    CalculateDescendants();//计算每个任务所有后代任务
    CalculateAncestors();//计算所有任务的父代任务
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    //    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    W_Cal_Average_S(ww);//计算平均执行时间，并将其保存到ww中
    //    C_Cal_Average(cc);
    Calculate_Rank_b_S(Rank_b, ww);// 计算任务优先级并保存到 Rank_b 中

    //找到Rank_b向量中的最大值
    double MaxRank_b = Rank_b[0];
    for (int i = 1; i < comConst.NumOfTsk; ++i) {
        if (MaxRank_b + PrecisionValue < Rank_b[i]) {
            MaxRank_b = Rank_b[i];
        }
    }

    vector<int> NumOfAncestors(comConst.NumOfTsk);
    vector<int> NumOfDescendants(comConst.NumOfTsk);
    vector<int> NumOfNonDescendants(comConst.NumOfTsk);

     //    #pragma omp parallel for
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfAncestors[i] = Ancestors[i].size();//父节点数量
        NumOfDescendants[i] = Descendants[i].size();//子节点数量
        NumOfNonDescendants[i] = comConst.NumOfTsk - NumOfDescendants[i];//非子节点数量
    }

    vector<vector<double>> PMR(comConst.NumOfTsk, vector<double>(comConst.NumOfRsc, 0));
    InitProModelOfResAlc(PMR);    //初始化资源分配概率模型 论文公式16 PMR[i][j] represents the probability that task i is assigned to resource j;
    vector<vector<double>> PMS(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    InitProModelOfTskSch(PMS, NumOfAncestors, NumOfNonDescendants); //初始化任务调度概率模型 论文公式17 PMS[i][k] represents the probability that the k-th scheduled task is task i

    vector<chromosome> Population(Parameter_TSEDA.NumOfChromPerPop);
    chromosome BstChrom = GnrChr_HMEC_S(Rank_b);//根据最小能耗生成初始最优个体
    chromosome Chrom_HEFT = GnrChr_HEFT_S(Rank_b);//根据HEFT最短时间策略生成个体

    // 如果 HEFT 个体的能耗更低，则更新最优个体
    if (Chrom_HEFT.EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption)
        BstChrom = Chrom_HEFT;

    vector<double> eta_TSO(comConst.NumOfTsk);// 初始化任务调度概率修正因子
    double RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;// 计算运行时间
    
    // 第一阶段：优化资源分配   迭代0.75*时间
    while (RunTime + PrecisionValue < Parameter_TSEDA.RunTimeRatioOfStg1 * SchTime) {

        double originBstEC = BstChrom.EnergyConsumption;
    //while (iteration < 0.5 * ITERATIONUM) {

        //多线程并行执行

        // 使用 OpenMP 并行计算任务调度概率修正因子
#pragma omp parallel for
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            //eta_TSO[i] = pow(当前任务优先级 / 最高优先级， （1 - 运行时间 / 总调度时间） * 调节参数）
            //pow(x, y) 表示计算 x 的 y 次幂。
            eta_TSO[i] = pow(Rank_b[i] / MaxRank_b, (1 - RunTime / SchTime) * Parameter_TSEDA.fdhi);
        }

        // 使用 OpenMP 并行生成种群并评估能耗
#pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfChromPerPop; ++n) {
            Population[n] = GnrTskLstOfChr(PMS, eta_TSO);// 生成任务调度列表
            GnrML_Evl_MEC_S(Population[n]);// 生成chromosome并评估个体能耗
        }

        // 更新最优个体
        for (int n = 0; n < Parameter_TSEDA.NumOfChromPerPop; ++n) {
            if (Population[n].EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption)
                BstChrom = Population[n];
        }

        UpdatePMR(PMR, BstChrom);  // 更新资源分配概率模型
        UpdatePMS(PMS, BstChrom);  // 更新任务调度概率模型
        cout << "iteration1: " << iteration << endl;
        ++iteration;// 增加迭代次数

        RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;// 更新运行时间

        std::ofstream outFile(ECfileName, std::ios::app);
        if (!outFile) {
            std::cerr << "无法打开文件: " << ECfileName << std::endl;
            break;
        }

        string algName = "TSEDA";
        if (BstChrom.EnergyConsumption + PrecisionValue < originBstEC)
        {
            outFile << algName << "\t" << Model << "\t" << comConst.NumOfTsk << "\t" << iteration << "\t" << RunTime << "\t"
                << BstChrom.EnergyConsumption << endl;
        }
        outFile.close();

    }

    
    // 第二阶段：优化任务调度 迭代剩余0.25 * 时间
    while (RunTime + PrecisionValue < SchTime) {

        double originBstEC = BstChrom.EnergyConsumption;
    //while (iteration < ITERATIONUM) {

        // 使用 OpenMP 并行计算任务调度概率修正因子
#pragma omp parallel for
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            eta_TSO[i] = pow(Rank_b[i] / MaxRank_b, (1 - RunTime / SchTime) * Parameter_TSEDA.fdhi);
        }

        // 使用 OpenMP 并行生成种群、分配资源并评估个体
#pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfChromPerPop; ++n) {
            Population[n] = GnrTskLstOfChr(PMS, eta_TSO);//生成任务调度列表
            GnrRscLstOfChr(Population[n], PMR);//评估个体能耗
            DcdEvl_S(Population[n]);// 解码并评估个体
            CalculateEnergy1(Population[n]);// 计算能耗
        }

        // 对种群按能耗排序
        sort(Population.begin(), Population.end(), SortByEnergyConsumption);

        // 使用 OpenMP 并行改进种群
#pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfImproveOfPop; ++n) {
            IFBD_S(Population[n]); //IFBSS;// 改进个体 IFBSS(Iterative Forward and Backward Scheduling Strategy)
            LBCA_S(Population[n]); //LBCAS;// 改进个体 LBCAS(Load Balance and Communication Avoidance Strategy)
        }

        // 更新最优个体
        for (int n = 0; n < Parameter_TSEDA.NumOfImproveOfPop; ++n) {
            if (Population[n].EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption)
                BstChrom = Population[n];
        }

        UpdatePMR(PMR, BstChrom); // 更新资源分配概率模型
        UpdatePMS(PMS, BstChrom);// 更新任务调度概率模型
        cout << "iteration2: " << iteration << endl;
        ++iteration;// 增加迭代次数
        RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;// 更新运行时间

        std::ofstream outFile(ECfileName2, std::ios::app);
        if (!outFile) {
            std::cerr << "无法打开文件: " << ECfileName2 << std::endl;
            break;
        }
        cout << BstChrom.EnergyConsumption << endl;
        string algName = "TSEDA";
        if (BstChrom.EnergyConsumption + PrecisionValue <= originBstEC)
        {
            outFile << algName << "\t" << Model << "\t" << comConst.NumOfTsk << "\t" << iteration << "\t" << RunTime << "\t"
                << BstChrom.EnergyConsumption << endl;
        }
        outFile.close();

    }
    SchTime = RunTime;// 更新调度时间
    return BstChrom;// 返回最优个体
}
