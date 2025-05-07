
#include "classAndVarDefine.h"
vector<Task> Tasks;
vector<Task> TaskO;
vector<int> Oid;
vector<int> Nid;
vector<vector<int>> TskLstInLvl;
vector<vector<double> > ParChildTranFileSizeSum;
vector<int> LevelIdOfTask;
ComConst comConst;
vector<set<int> > Descendants;
vector<set<int> > Ancestors;
vector<Resource> Rscs;
vector<HT> HstSet;
double ModelScale;
double SchTime;
Paramet_HGA Parameter_HGA;
Paramet_NGA Parameter_NGA;
Paramet_LWSGA Parameter_LWSGA;
Paramet_CGA Parameter_CGA;
Paramet_HPSO Parameter_HPSO;
Paramet_TSEDA Parameter_TSEDA;
Paramet_ADBRKGA Parameter_ADBRKGA;
Paramet_QPHH Parameter_QPHH;
