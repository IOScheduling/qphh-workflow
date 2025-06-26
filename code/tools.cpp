
#include "common.h"
#include "GenOperator.h"
#include "GenerateAchrom.h"
#include "tools.h"
using namespace std;



void CalculateLevelList() {
    LevelIdOfTask.resize(comConst.NumOfTsk); 
    vector<int> InDegree;   
    vector<int> stk;       
    InDegree.assign(comConst.NumOfTsk, 0);

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        InDegree[i] = Tasks[i].parents.size();
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (InDegree[i] == 0) stk.push_back(i);
    }

    int MaxLevel = 0;
    while (!stk.empty()) {
        int v = stk[0];
        LevelIdOfTask[v] = 0;

        for (int i = 0; i < Tasks[v].parents.size(); ++i) {
            if (LevelIdOfTask[Tasks[v].parents[i]] >= LevelIdOfTask[v]) {
                LevelIdOfTask[v] = LevelIdOfTask[Tasks[v].parents[i]] + 1;
            }
        }
        if (LevelIdOfTask[v] + 1 > MaxLevel) {
            MaxLevel = LevelIdOfTask[v] + 1;
            TskLstInLvl.resize(MaxLevel);
        }
        TskLstInLvl[LevelIdOfTask[v]].push_back(v);
        stk.erase(stk.begin());

        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            InDegree[Tasks[v].children[i]]--;
            if (InDegree[Tasks[v].children[i]] == 0) {
                stk.push_back(Tasks[v].children[i]);
            }
        }
    }
}


void CalculateDescendants() {
    Descendants.resize(comConst.NumOfTsk);

    for (int i = TskLstInLvl.size() - 2; i >= 0; --i) {
        for (int taskId : TskLstInLvl[i]) {
            for (int childId : Tasks[taskId].children) {
                Descendants[taskId].insert(childId);
                Descendants[taskId].insert(Descendants[childId].begin(), Descendants[childId].end());
            }
        }
    }
}


void CalculateAncestors() {
    Ancestors.resize(comConst.NumOfTsk);
    for (int i = 1; i < TskLstInLvl.size(); ++i) {
        for (int taskId : TskLstInLvl[i]) {
            for (int parentId : Tasks[taskId].parents) {
                Ancestors[taskId].insert(parentId);
                Ancestors[taskId].insert(Ancestors[parentId].begin(), Ancestors[parentId].end());
            }
        }
    }
}

void IndexSortByValueOnAscend(vector<int>& ind, vector<double>& value) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] < value[v2]; });
}

void IndexSortByValueOnAscend(vector<int>& ind, vector<int>& value) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] < value[v2]; });
}

void IndexSortByValueOnDescend(vector<int>& ind, vector<double>& value) {
    vector<double> result;
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] > value[v2]; });
}


bool SortValueByDescend(pair<int, double>& a, pair<int, double>& b) {
    return a.second > b.second + PrecisionValue;
}

bool SortByEnergyConsumptionByDescend(chromosome& a, chromosome& b) {
    return a.EnergyConsumption > b.EnergyConsumption + PrecisionValue;
}

bool SortByEnergyConsumption(chromosome& a, chromosome& b) {
    return a.EnergyConsumption + PrecisionValue < b.EnergyConsumption;
}


//{generate a random number in [Start,End) }
double RandomDouble(int start, int end) {
    double ret = rand() % ((end - start) * 1000) / 1000.0 + start;
    return ret;
}

//{generate a random number in [Start,End] }
double RandomDouble2(int start, int end) {
    double ret = start + rand() % (end - start) + rand() % 1001 / 1000.0;
    return ret;
}

double CalculatePowerByLoad(double ld, int HTid) {
    if (HTid % 3 == 0) {
        // {NEC Corporation Express5800/GT110f-S}
        if (ld <= 0.1) {
            return ld * 44 + 15.9;
        }
        else if (ld <= 0.2) {
            return ld * 21 + 18.2;
        }
        else if (ld <= 0.3) {
            return ld * 20 + 18.4;
        }
        else if (ld <= 0.4) {
            return ld * 28 + 16;
        }
        else if (ld <= 0.5) {
            return ld * 26 + 16.8;
        }
        else if (ld <= 0.6) {
            return ld * 32 + 13.8;
        }
        else if (ld <= 0.7) {
            return ld * 38 + 10.2;
        }
        else if (ld <= 0.8) {
            return ld * 27 + 17.9;
        }
        else if (ld <= 0.9) {
            return ld * 31 + 14.7;
        }
        else if (ld <= 1.0 + PrecisionValue) {
            return ld * 25 + 20.1;
        }
        else {
            cout << "load can not is larger than 1" << endl; exit(0);
        }
    }
    else if (HTid % 3 == 1) {
        //{FUJITSU Server PRIMERGY RX1330 M3}
        if (ld <= 0.1) {
            return ld * 44 + 13.1;
        }
        else if (ld <= 0.2) {
            return ld * 30 + 14.5;
        }
        else if (ld <= 0.3) {
            return ld * 25 + 15.5;
        }
        else if (ld <= 0.4) {
            return ld * 24 + 15.8;
        }
        else if (ld <= 0.5) {
            return ld * 28 + 14.2;
        }
        else if (ld <= 0.6) {
            return ld * 37 + 9.7;
        }
        else if (ld <= 0.7) {
            return ld * 54 - 0.5;
        }
        else if (ld <= 0.8) {
            return ld * 72 - 13.1;
        }
        else if (ld <= 0.9) {
            return ld * 70 - 11.5;
        }
        else if (ld <= 1.0 + PrecisionValue) {
            return ld * 46 + 10.1;
        }
        else {
            cout << "load can not is larger than 1" << endl; exit(0);
        }
    }
    else if (HTid % 3 == 2) {
        //{L: Lenovo Global Technology ThinkSystem SR150}-1
        if (ld <= 0.1) {
            //return ld * 3.2 + 16.8;
            return ld * 32 + 16.8;//20
        }
        else if (ld <= 0.2) {
            return ld * 28 + 17.2;//22.8
        }
        else if (ld <= 0.3) {
            return ld * 32 + 16.4;//26
        }
        else if (ld <= 0.4) {
            return ld * 39 + 14.3;//29.9
        }
        else if (ld <= 0.5) {
            return ld * 51 + 9.5;//35
        }
        else if (ld <= 0.6) {
            return ld * 54 + 8;//40.4
        }
        else if (ld <= 0.7) {
            return ld * 59 + 5;//46.3
        }
        else if (ld <= 0.8) {
            return ld * 93 - 18.8;//55.6
        }
        else if (ld <= 0.9) {
            return ld * 132 - 50;//68.8
        }
        else if (ld <= 1.0 + PrecisionValue) {
            return ld * 184 - 96.8; // 87.2
        }
        else {
            cout << "load can not is larger than 1" << endl; exit(0);
        }
        //{F: FUJITSU Server PRIMERGY TX1330 M4}-2
//        if (ld <= 0.1) {
//            return ld*33 + 14.9;
//        } else if (ld <= 0.2) {
//            return ld*18 + 16.4;
//        } else if (ld <= 0.3) {
//            return ld*28 + 14.4;
//        } else if (ld <= 0.4) {
//            return ld*35 + 12.3;
//        } else if (ld <= 0.5) {
//            return ld*45 + 8.3;
//        } else if (ld <= 0.6) {
//            return ld*50 + 5.8;
//        } else if (ld <= 0.7) {
//            return ld*68 - 5;
//        } else if (ld <= 0.8) {
//            return ld*84 - 16.2;
//        } else if (ld <= 0.9) {
//            return ld*104 - 32.2;
//        } else if (ld <= 1.0 + PrecisionValue) {
//            return ld*115 - 42.1;
//        } else {
//            cout << "load can not is larger than 1" << endl; exit(0);
//        }
    }
    else {
        cout << "host configurations are wrong!";  exit(0);
    }

}
