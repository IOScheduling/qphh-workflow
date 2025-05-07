
#include "config.h"
#include <fstream>
#include <sstream>
#include "GenerateAchrom.h"
#include "pugixml.hpp"
#include "tools.h"
void DeleteFirstLineInFile(string fileName) {
    vector<string> VecContent;
    string StrLine;
    ifstream iFile(fileName);
    if (!iFile) {
        cout << fileName << "filelist open failed!\n";
        exit(0);
    }
    //Read all content in the document into "VecContent"
    while (iFile) {
        getline(iFile, StrLine);
        VecContent.push_back(StrLine);
    }
    iFile.close();
    VecContent.erase(VecContent.begin()); // delete the first line
    ofstream oFile(fileName);
    if (!oFile) {
        cout << fileName << "filelist open failed!\n";
        exit(0);
    }
    vector<string>::const_iterator iter = VecContent.begin();
    //{Rewrite the contents of "vecContent" into the file}.
    for (; VecContent.end() != iter; ++iter) {
        oFile.write((*iter).c_str(), (*iter).size());
        oFile << '\n';
    }
    oFile.close();
}

int readID(string id) {
    int ret = 0;
    for (int i = 0; i < id.length(); i++) {
        if (id[i] <= '9' && id[i] >= '0')
            ret = ret * 10 + id[i] - '0';
    }
    return ret;
}

//read information
void ReadFile(string XmlFile, string RscAlcFile) {
    string FilePath = "../data";
    string XmlPath = FilePath + "/" + XmlFile;
    pugi::xml_document doc;
    doc.load_file((XmlPath).c_str());
    pugi::xml_node root = doc.child("adag");
    for (pugi::xml_node job = root.child("job"); job; job = job.next_sibling("job")) {
        Task task = Task();
        task.length = XY_MAX(fabs(atof(job.attribute("runtime").value())), 0.001);// read the length of task
        //{read file}
        for (pugi::xml_node uses = job.child("uses"); uses; uses = uses.next_sibling("uses")) {
            vfile file = vfile();
            file.source = -1;
            if (!strcmp(uses.attribute("link").value(), "input")) {    //read input file
                file.FileName = uses.attribute("file").value();
                file.size = fabs(atof(uses.attribute("size").value()));
                task.IFile.push_back(file);
                task.OrignalInputFileSizeSum = task.OrignalInputFileSizeSum + file.size;    //calculate the size of all input files
            }
            else {                                                          //read output file
                file.FileName = uses.attribute("file").value();
                file.size = fabs(atof(uses.attribute("size").value()));
                task.OFile.push_back(file);
                task.OFileSizeSum = task.OFileSizeSum + file.size;    //calculate the size of all output file
            }
        }
        Tasks.push_back(task);
    }
    comConst.NumOfTsk = Tasks.size();

    // read the info of relation between tasks
    for (pugi::xml_node child = root.child("child"); child; child = child.next_sibling("child")) {
        int ChildIndex = readID(child.attribute("ref").value());
        for (pugi::xml_node parent = child.child("parent"); parent; parent = parent.next_sibling("parent")) {
            int ParentIndex = readID(parent.attribute("ref").value());
            Tasks[ChildIndex].parents.push_back(ParentIndex);
            Tasks[ParentIndex].children.push_back(ChildIndex);
        }
    }

    //{calculate the transfer data size among tasks}
    ParChildTranFileSizeSum.resize(comConst.NumOfTsk);
    for (int k = 0; k < comConst.NumOfTsk; ++k) {
        ParChildTranFileSizeSum[k].resize(comConst.NumOfTsk, 0);
    }
    for (int i = 0; i < Tasks.size(); ++i) {
        if (Tasks[i].parents.size() == 0) {
            Tasks[i].ExternalInputFileSizeSum = Tasks[i].OrignalInputFileSizeSum;
            Tasks[i].IFileSizeSum = Tasks[i].OrignalInputFileSizeSum;   //revised-xy-2022.2.17
            continue;
        }
        for (int p = 0; p < Tasks[i].IFile.size(); ++p) {               //two loop (for p and j) can be switched in order
            string IName = Tasks[i].IFile[p].FileName;
            int flag = 0;
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {         //Traverse the parent task
                int Parent = Tasks[i].parents[j];
                for (int q = 0; q < Tasks[Parent].OFile.size(); ++q) {  //Traverse the output files of the parent task
                    string OName = Tasks[Parent].OFile[q].FileName;
                    if (IName.compare(OName) == 0) {                // judge whether two file names are the same; 0: same; -1: not same
                        ParChildTranFileSizeSum[Parent][i] += Tasks[i].IFile[p].size;
                        //If multiple identical files from different parent tasks are transferred to the same child task, the "source" records the last parent task
                        Tasks[i].IFile[p].source = Parent;
                        ++flag;
                        break;
                    }
                }
            }
            if (flag == 0) {
                Tasks[i].ExternalInputFileSizeSum += Tasks[i].IFile[p].size;
            }

        }
        for (int Parent : Tasks[i].parents) {
            Tasks[i].IFileSizeSum += ParChildTranFileSizeSum[Parent][i];
        }
        Tasks[i].IFileSizeSum += Tasks[i].ExternalInputFileSizeSum;
       
    }

    int nmbOfRscCnf = 1;

      /* if (XmlFile == "CyberShake_30.xml"){
        nmbOfRscCnf = 5;
    } else if (XmlFile == "CyberShake_50.xml"){

    } else if (XmlFile == "CyberShake_100.xml"){
        nmbOfRscCnf = 17;
    } else if (XmlFile == "CyberShake_1000.xml") {
           nmbOfRscCnf = 170;
    } else if (XmlFile == "Epigenomics_24.xml"){
        nmbOfRscCnf = 2;
    } else if (XmlFile == "Epigenomics_47.xml"){
        nmbOfRscCnf = 4;
    } else if (XmlFile == "Epigenomics_100.xml"){
        nmbOfRscCnf = 8;
    } else if (XmlFile == "Epigenomics_997.xml") {
           nmbOfRscCnf = 80;
    } else if (XmlFile == "Ligo_30.xml"){
        nmbOfRscCnf = 3;
    } else if (XmlFile == "Ligo_50.xml"){
        nmbOfRscCnf = 4;
    } else if (XmlFile == "Ligo_100.xml"){
        nmbOfRscCnf = 8;
    } else if (XmlFile == "Ligo_1000.xml") {
        nmbOfRscCnf = 80;
    } else if (XmlFile == "Montage_25.xml"){
        nmbOfRscCnf = 3;
    } else if (XmlFile == "Montage_50.xml"){
        nmbOfRscCnf = 10;
    } else if (XmlFile == "Montage_100.xml"){
        nmbOfRscCnf = 21;
    } else if (XmlFile == "Montage_1000.xml") {
       nmbOfRscCnf = 210;
    } else if (XmlFile == "Sipht_29.xml"){
        nmbOfRscCnf = 7;
    } else if (XmlFile == "Sipht_58.xml"){
        nmbOfRscCnf = 14;
    } else if (XmlFile == "Sipht_97.xml"){
        nmbOfRscCnf = 25;
}    else if (XmlFile == "Sipht_1000.xml") {
        nmbOfRscCnf = 250;
    } else {
        cout << endl << "XmlFile does not exit!";
    }*/

    for (int k = 0; k < nmbOfRscCnf; ++k) {
        Resource Rsc_0 = Resource(0 + 3 * k, 2.0, 2000);
        Resource Rsc_1 = Resource(1 + 3 * k, 2.0, 2000);
        Resource Rsc_2 = Resource(2 + 3 * k, 2.0, 2000);
        Resource Rsc_3 = Resource(0 + 3 * k, 4.0, 4000);
        Resource Rsc_4 = Resource(0 + 3 * k, 4.0, 4000);
        Resource Rsc_5 = Resource(1 + 3 * k, 4.0, 4000);
        Resource Rsc_6 = Resource(2 + 3 * k, 4.0, 4000);
        Resource Rsc_7 = Resource(1 + 3 * k, 8.0, 8000);
        Resource Rsc_8 = Resource(2 + 3 * k, 8.0, 8000);
        Resource Rsc_9 = Resource(2 + 3 * k, 8.0, 8000);
        Rscs.push_back(Rsc_0);
        Rscs.push_back(Rsc_1);
        Rscs.push_back(Rsc_2);
        Rscs.push_back(Rsc_3);
        Rscs.push_back(Rsc_4);
        Rscs.push_back(Rsc_5);
        Rscs.push_back(Rsc_6);
        Rscs.push_back(Rsc_7);
        Rscs.push_back(Rsc_8);
        Rscs.push_back(Rsc_9);
        vector<int> VMlist0 = { 0 + 10 * k, 3 + 10 * k, 4 + 10 * k };
        vector<int> VMlist1 = { 1 + 10 * k, 5 + 10 * k, 7 + 10 * k };
        vector<int> VMlist2 = { 2 + 10 * k, 6 + 10 * k, 8 + 10 * k, 9 + 10 * k };
        HT ht0 = HT(VMlist0, 10, 10000);
        HT ht1 = HT(VMlist1, 14, 16000);
        HT ht2 = HT(VMlist2, 22.2, 24000);
        HstSet.push_back(ht0);
        HstSet.push_back(ht1);
        HstSet.push_back(ht2);
    }
    comConst.NumOfRsc = Rscs.size();
    //    cout << endl << "nmbOfRsc = " << comConst.NumOfRsc;

        //read the RscAlc file to task data structure
        //in the RscAlc file, each resource can perform at least one task and each task can be performed by at least one resource
    char line[8192] = { 0 };
    int TskIndex = 0;
    int RscIndex = 0;
    string RscAlcPath = FilePath + "/" + RscAlcFile;  //RscAlcPath xy4
    ifstream fin(RscAlcPath, ios::in);
    if (!fin) {
        cout << "Error at open Rsc file" << endl;
        exit(0);
    }
    else {
        while (fin.getline(line, sizeof(line))) {
            stringstream Word(line);
            while (1) {
                Word >> TskIndex;
                if (Word.fail()) break;
                Tasks[TskIndex].ElgRsc.push_back(RscIndex);
                Rscs[RscIndex].ElgTsk.push_back(TskIndex);
            }
            ++RscIndex;
            if (RscIndex == comConst.NumOfRsc) break;
        }
        if (RscIndex < comConst.NumOfRsc) {
            cout << endl << "the number of lines in the RscAlcFile is less then the number of resources!";
            exit(0);
        }
    }
    ModelScale = 0;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ModelScale += Tasks[i].ElgRsc.size();
    }
    fin.close();
}

void ConfigParameter_LWSGA() {
    Parameter_LWSGA.NumOfChromPerPop = 70;
    Parameter_LWSGA.CrossoverRate = 0.7;
}

void ConfigParameter_NGA() {
    Parameter_NGA.NumOfChromPerPop = 4 * Tasks.size();
    Parameter_NGA.MutationRate = 0.2;
    Parameter_NGA.CrossoverRate = 0.8;
    Parameter_NGA.EliteRate = 0.2;
}

void ConfigParameter_HGA() {
    Parameter_HGA.NumOfChromPerPop = 100;
    Parameter_HGA.MutationRate = 0.02;
    Parameter_HGA.EliteRate = 0.2;
}

void ConfigParameter_HPSO() {
    Parameter_HPSO.NumOfChromPerPop = 20;
    Parameter_HPSO.InertiaWeight = 0.5;
    Parameter_HPSO.c1 = 2;
    Parameter_HPSO.c2 = 2;
}

void ConfigParameter_CGA() {
    Parameter_CGA.NumOfChromPerPop = 2 * Tasks.size();
    Parameter_CGA.MutationRate = 0.2;
    Parameter_CGA.CrossoverRate = 0.8;
}

void ConfigParameter_TSEDA() {
    Parameter_TSEDA.NumOfChromPerPop = Tasks.size() * 1.8;
    if (Parameter_TSEDA.NumOfChromPerPop % 2 == 1) {
        ++Parameter_TSEDA.NumOfChromPerPop;
    }
    Parameter_TSEDA.theta1 = 0.35;
    Parameter_TSEDA.theta2 = 0.25;
    Parameter_TSEDA.fdhi = 0.8;
    Parameter_TSEDA.NumOfImproveOfPop = ceil(Parameter_TSEDA.NumOfChromPerPop * 0.03);
    Parameter_TSEDA.RunTimeRatioOfStg1 = 0.75; //0.80;
}

void ConfigParameter_ADBRKGA() {
    Parameter_ADBRKGA.NumOfChromPerPop = Tasks.size() * 2.5;        //set the population size
    if (Parameter_ADBRKGA.NumOfChromPerPop % 2 == 1) {
        ++Parameter_ADBRKGA.NumOfChromPerPop;
    }
    Parameter_ADBRKGA.alpha = 3;
    Parameter_ADBRKGA.beta = 0.8;
    Parameter_ADBRKGA.BiasesRate = 0.75;
    Parameter_ADBRKGA.ImmigrationRate = 0.07;
    Parameter_ADBRKGA.ImprovementRate = 0.07;
}
void ConfigParameter_QPHH() {
    Parameter_QPHH.NumOfChromPerPop = Tasks.size() * 1.8;
    if (Parameter_TSEDA.NumOfChromPerPop % 2 == 1) {
        ++Parameter_TSEDA.NumOfChromPerPop;
    }
}

void ClearALL() {
    Tasks.clear();
    TskLstInLvl.clear();
    LevelIdOfTask.clear();
    Rscs.clear();
    ParChildTranFileSizeSum.clear();
    Descendants.clear();
    Ancestors.clear();
    HstSet.clear();
}