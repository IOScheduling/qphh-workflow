#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "NGA.h"
#include "CGA.h"
#include "HGA.h"
#include "TSEDA.h"
#include "LWSGA.h"
#include "HEFT.h"
#include "HMEC.h"
#include "HPSO.h"
#include "ADBRKGA.h"
#include "QPHH.h"

using namespace std;

void saveChromosomeToFile(const chromosome& Chrom, const string& filename) {
    ofstream outFile(filename);

    if (!outFile) {
        cerr << "Unable to open file for writing: " << filename << endl;
        return;
    }

    outFile << "Resource List(RscAlcLst): ";
    for (int vm: Chrom.RscAlcLst) {
        outFile << vm << " " ;
    }
    outFile << endl;

    outFile << "Task Schedule List (TskSchLst): ";

    vector<int> taskPosition(Chrom.TskSchLst.size());
    for (int i = 0; i < Chrom.TskSchLst.size(); ++i) {
        taskPosition[Chrom.TskSchLst[i]] = i;
    }

    for (int task : taskPosition) {
        outFile << task << " ";
    }
    outFile << endl;

    outFile << "Start Times: ";
    for (double startTime : Chrom.StartTime) {
        outFile << startTime << " ";
    }
    outFile << endl;

    outFile << "End Times: ";
    for (double endTime : Chrom.EndTime) {
        outFile << endTime << " ";
    }
    outFile << endl;

    outFile << "MakeSpan: " << Chrom.MakeSpan << endl;
    outFile << "Energy Consumption: " << Chrom.EnergyConsumption << endl;

    outFile.close();
    cout << "Chromosome data has been saved to " << filename << endl;
}

int main() {



    srand((int)time(0));
    map<string, double> SchTime;
    //CDT=20

    SchTime["CyberShake30_1.0"] = 0.407;
    SchTime["CyberShake50_1.0"] = 0.51;
    SchTime["CyberShake100_1.0"] = 2.747;
    SchTime["CyberShake500_1.0"] = 28.292;
    SchTime["CyberShake1000_1.0"] = 270.906;
    SchTime["CyberShake1500_1.0"] = 648.458;
    SchTime["Epigenomics24_1.0"] = 0.297;
    SchTime["Epigenomics47_1.0"] = 0.551;
    SchTime["Epigenomics100_1.0"] = 5.021;
    SchTime["Epigenomics497_1.0"] = 57.266;
    SchTime["Epigenomics997_1.0"] = 495.177;
    SchTime["Epigenomics1499_1.0"] = 1234.469;
    SchTime["Ligo30_1.0"] = 0.181;
    SchTime["Ligo50_1.0"] = 0.222;
    SchTime["Ligo100_1.0"] = 1.262;
    SchTime["Ligo500_1.0"] = 60.009;
    SchTime["Ligo1000_1.0"] = 214.492;
    SchTime["Ligo1500_1.0"] = 298.262;
    SchTime["Montage25_1.0"] = 0.185;
    SchTime["Montage50_1.0"] = 0.535;
    SchTime["Montage100_1.0"] =0.613;
    SchTime["Montage500_1.0"] = 12.758;
    SchTime["Montage1000_1.0"] = 144.047;
    SchTime["Montage1500_1.0"] = 266.813;
    SchTime["Sipht29_1.0"] = 0.303;
    SchTime["Sipht58_1.0"] = 0.742;
    SchTime["Sipht97_1.0"] = 3.046;
    SchTime["Sipht484_1.0"] = 61.818;
    SchTime["Sipht968_1.0"] = 426.015;
    SchTime["Sipht1452_1.0"] = 1484.371;

    ofstream outfile("../result.txt", ios::out);
    outfile.close();
    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
        }
        getline(iFile, StrLine);
        if (StrLine.size() < 1) {
            cout << endl << "Empty input file" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(StrLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + ".xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        cout << endl << Model << " " << NumOfTask << " " << RscAvlRatio << " ";

        string ITECfileName = "ITEC.txt";



        double HEFT_SchTime = 0;
        string HEFT_ECfileName = "../HEFT_" + ITECfileName;
        chromosome HEFT_Result = runHEFT(HEFT_ECfileName, XmlFile, RscAlcFile, HEFT_SchTime);
        ClearALL();

        double HGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HGA_Iteration = 0;
        string HGA_ECfileName = "../HGA_" + ITECfileName;
        chromosome HGA_Result = runHGA(Model, HGA_ECfileName, XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);
        ClearALL();

        double NGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int NGA_Iteration = 0;
        string NGA_ECfileName = "../NGA_" + ITECfileName;
        chromosome NGA_Result = runNGA(Model, NGA_ECfileName, XmlFile, RscAlcFile, NGA_SchTime, NGA_Iteration);
        ClearALL();

        double LWSGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int LWSGA_Iteration = 0;
        string LWSGA_ECfileName = "../LWSGA_" + ITECfileName;
        chromosome LWSGA_Result = runLWSGA(Model,LWSGA_ECfileName, XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);
        ClearALL();

        double ADBRKGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int ADBRKGA_Iteration = 0;
        string ADBRKGA_ECfileName = "../ADBRKGA_" + ITECfileName;
        chromosome ADBRKGA_Result = runADBRKGA(Model,ADBRKGA_ECfileName, XmlFile, RscAlcFile, ADBRKGA_SchTime, ADBRKGA_Iteration);
        ClearALL();

        double HPSO_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HPSO_Iteration = 0;
        string HPSO_ECfileName = "../HPSO_" + ITECfileName;
        chromosome HPSO_Result = runHPSO(Model, HPSO_ECfileName, XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration);
        ClearALL();

        double TSEDA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int TSEDA_Iteration = 0;
        string TSEDA_ECfileName = "../TSEDA_" + ITECfileName;
        chromosome TSEDA_Result = runTSEDA(Model, TSEDA_ECfileName, XmlFile, RscAlcFile, TSEDA_SchTime, TSEDA_Iteration);
        ClearALL();
        
        double QPHH_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int QPHH_Iteration = 0;
        double popSize, epslion, gamma;
        int selectCount;
        popSize = 50;
        epslion = 0.3;
        selectCount = 15;
        gamma = 0.5;
        string QPHH_ECfileName = "../QPHH_" + ITECfileName;
        chromosome QPHH_Result = runQPHH(Model, QPHH_ECfileName,XmlFile, RscAlcFile, QPHH_SchTime, QPHH_Iteration, popSize,selectCount, epslion,gamma);
        ClearALL();


        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }

        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);

        outfile << Model << "\t" << NumOfTask << "\t" << RscAvlRatio << "\t "
            << HEFT_Result.EnergyConsumption << "\t"  << HGA_Result.EnergyConsumption << "\t" << NGA_Result.EnergyConsumption << "\t"
            << LWSGA_Result.EnergyConsumption << "\t" << ADBRKGA_Result.EnergyConsumption << "\t" << HPSO_Result.EnergyConsumption << "\t"
            << TSEDA_Result.EnergyConsumption << "\t"   << QPHH_Result.EnergyConsumption << "\t"
            << HEFT_Result.MakeSpan << "\t"<< HGA_Result.MakeSpan << "\t" << NGA_Result.MakeSpan << "\t"
            << LWSGA_Result.MakeSpan << "\t" << ADBRKGA_Result.MakeSpan << "\t" << HPSO_Result.MakeSpan << "\t "
            << TSEDA_Result.MakeSpan << "\t " << QPHH_Result.MakeSpan << "\t"  
            << HEFT_SchTime << "\t "  << HGA_SchTime << "\t" << NGA_SchTime << "\t "
            << LWSGA_SchTime << "\t " << ADBRKGA_SchTime << "\t" << HPSO_SchTime << "\t "
            << TSEDA_SchTime << "\t " << QPHH_SchTime << "\t"
            << HGA_Iteration << "\t" << NGA_Iteration << "\t"  << LWSGA_Iteration << "\t"
            << ADBRKGA_Iteration << "\t " << HPSO_Iteration << "\t " << TSEDA_Iteration << "\t" << QPHH_Iteration << "\t"
            << endl;

        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}
