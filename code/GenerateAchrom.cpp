

#include "GenerateAchrom.h"
#include <cstdlib>
#include "tools.h"
#include "common.h"
#include "GenOperator.h"
#include "IntervalTree.h"
#include <thread>
#include <chrono>



struct Event {
    double time;   
    bool isStart;  
    int taskId;    
    bool operator<(const Event& other) const {
        return time < other.time;
    }
};




void W_Cal_Average_S(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0;
        double AllTransferData = Tasks[i].IFileSizeSum + Tasks[i].OFileSizeSum;
        for (int RscId : Tasks[i].ElgRsc)
            sum += Tasks[i].length / Rscs[RscId].pc + AllTransferData / VALUE * 8 / Rscs[RscId].bw;
            
        w[i] = sum / Tasks[i].ElgRsc.size();
    }
}
//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (Tasks[i].parents.size() == 0) {
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    }
                    else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double)(Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}



void Calculate_Rank_Rand_S(vector<double>& RankList, vector<double>& RandPri) {


    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size() - 1].size(); ++i) {
        int TaskId = TskLstInLvl[TskLstInLvl.size() - 1][i];
        RankList[TaskId] = RandPri[TaskId];
    }

    for (int i = TskLstInLvl.size() - 2; i >= 0; --i) {

        for (int TaskId : TskLstInLvl[i]) {


            for (int Child : Tasks[TaskId].children) {
                if (RankList[TaskId] + PrecisionValue < RankList[Child]) {
                    RankList[TaskId] = RankList[Child];
                }
            }
            RankList[TaskId] = RankList[TaskId] + RandPri[TaskId];
        }
    }
}




void Calculate_Rank_b_S(vector<double>& RankList, vector<double>& ExeTime) {

    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size() - 1].size(); ++i) {
        int TaskId = TskLstInLvl[TskLstInLvl.size() - 1][i];
        RankList[TaskId] = ExeTime[TaskId];
    }


    for (int i = TskLstInLvl.size() - 2; i >= 0; --i) {
        for (int TaskId : TskLstInLvl[i]) {
            for (int Child : Tasks[TaskId].children) {
                if (RankList[TaskId] + PrecisionValue < RankList[Child]) {
                    RankList[TaskId] = RankList[Child];
                }
            }
            RankList[TaskId] = RankList[TaskId] + ExeTime[TaskId];
        }
    }
}

void Calculate_Rank_t_S(vector<double>& RankList, vector<double>& ExeTime) {
    for (int i = 1; i < TskLstInLvl.size(); ++i) {
        for (int TaskId : TskLstInLvl[i]) {
            for (int Parent : Tasks[TaskId].parents) {
                double tem = RankList[Parent] + ExeTime[Parent];
                if (RankList[TaskId] < tem) {
                    RankList[TaskId] = tem;
                }
            }
        }
    }
}

chromosome GnrChr_HEFT_t_S(vector<double> Rank_t) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_t);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[i];
    }
    GnrML_Evl_EFT_S(TemChrom);
    CalculateEnergy1(TemChrom);
    return TemChrom;
}

//{in the SS, tasks are arranged according to the level form small to large, and the tasks in the same level are arranged in descend of rank_b_t}
chromosome GnrChr_HEFT_b_t_S(vector<double> Rank_b_t) {
    chromosome TemChrom;
    IntChr(TemChrom);
    int cur = 0;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        vector<double> TemList(TskLstInLvl[i].size());
        for (int j = 0; j < TskLstInLvl[i].size(); ++j)
            TemList[j] = Rank_b_t[TskLstInLvl[i][j]];
        vector<int> ind(TskLstInLvl[i].size());
        IndexSortByValueOnAscend(ind, TemList);
        for (int j = TskLstInLvl[i].size() - 1; j > -1; j--)
            TemChrom.TskSchLst[cur++] = TskLstInLvl[i][ind[j]];
    }
    GnrML_Evl_EFT_S(TemChrom);
    CalculateEnergy1(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HEFT_S(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_b);

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrML_Evl_EFT_S(TemChrom);
    CalculateEnergy1(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HEFT_b_ADBRKGA_S(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    vector<double> Decimals = GnrDecimalsByAscend();
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.Code_RK[ind[comConst.NumOfTsk - i - 1]] = Decimals[i];
    }
    HrsDcd_EFT_S(TemChrom);
    CalculateEnergy1(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HEFT_t_ADBRKGA_S(vector<double> Rank_t) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    vector<double> Decimals = GnrDecimalsByAscend();
    IndexSortByValueOnAscend(ind, Rank_t);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.Code_RK[ind[i]] = Decimals[i];
    }
    HrsDcd_EFT_S(TemChrom);
    CalculateEnergy1(TemChrom);
    return TemChrom;
}

chromosome GnrChr_DIHEFT_S(vector<double>& Rank_b) {
    chromosome chrom;
    IntChr(chrom);
    vector<int> upr(comConst.NumOfTsk, 0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }
    vector<set<double>> ITL;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    int IndexCount = 0;
    vector<double > AvrRtSet(comConst.NumOfTsk, 0.0);
    while (!RTI.empty()) {
        double max_RANK = 0;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            double priority = Rank_b[*lit] + AvrRtSet[*lit];
            if (priority - PrecisionValue > max_RANK) {
                pit = lit; max_RANK = priority;
            }
        }
        int CrnTskId = *pit;
        RTI.erase(pit);
        chrom.TskSchLst[IndexCount] = CrnTskId;
        IndexCount++;
        int FinalRscForCrnTask = -1, NeedProcessChildTask = -1;
        double FinalStartTimeOfCrnTask = 0, FinalEndTimeOfCrnTask = InfiniteValue, TemMaxRnk = -1;
        vector<int> ReadyChildTaskSet;
        for (int childId : Tasks[CrnTskId].children) {
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                if (Rank_b[childId] - PrecisionValue > TemMaxRnk) {
                    NeedProcessChildTask = childId;
                    TemMaxRnk = Rank_b[NeedProcessChildTask];
                }
                ReadyChildTaskSet.push_back(childId);
            }
        }
        if (NeedProcessChildTask == -1) {
            SelectRsc_EFT_S(chrom, ITL, CrnTskId, FinalRscForCrnTask, FinalStartTimeOfCrnTask, FinalEndTimeOfCrnTask);
            chrom.StartTime[CrnTskId] = FinalStartTimeOfCrnTask;
            chrom.EndTime[CrnTskId] = FinalEndTimeOfCrnTask;
            chrom.RscAlcLst[CrnTskId] = FinalRscForCrnTask;
            chrom.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].first = XY_MIN(chrom.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].first, chrom.StartTime[CrnTskId]);
            chrom.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].second = XY_MAX(chrom.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].second, chrom.EndTime[CrnTskId]);
            chrom.MakeSpan = XY_MAX(chrom.MakeSpan, chrom.EndTime[CrnTskId]);
            UpdateITL(ITL[FinalRscForCrnTask], FinalStartTimeOfCrnTask, FinalEndTimeOfCrnTask);
        }
        else {  
            int FinalRscForChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = InfiniteValue;
            vector<set<double>> ITL_CrnTaskScheduled;
            for (int RscForCrnTsk : Tasks[CrnTskId].ElgRsc) {
                double ReadyTimeOfCrnTask = 0, TransferDataSize = Tasks[CrnTskId].ExternalInputFileSizeSum;
                chrom.RscAlcLst[CrnTskId] = RscForCrnTsk;
                for (int i2 = 0; i2 < Tasks[CrnTskId].parents.size(); ++i2) {
                    int ParentTask = Tasks[CrnTskId].parents[i2];
                    int RscOfParentTask = chrom.RscAlcLst[ParentTask];
                    if (RscForCrnTsk != RscOfParentTask) {
                        TransferDataSize += ParChildTranFileSizeSum[ParentTask][CrnTskId];
                    }
                    if (ReadyTimeOfCrnTask + PrecisionValue < chrom.EndTime[ParentTask]) {
                        ReadyTimeOfCrnTask = chrom.EndTime[ParentTask];
                    }
                }
                double ExeTimeOfCrnTask = Tasks[CrnTskId].length / Rscs[RscForCrnTsk].pc + (TransferDataSize + Tasks[CrnTskId].OFileSizeSum) / VALUE * 8 / Rscs[RscForCrnTsk].bw;
                double CrnTaskStartTime = FindIdleTimeSlot(ITL[RscForCrnTsk], ExeTimeOfCrnTask, ReadyTimeOfCrnTask);
                double CrnTaskEndTime = CrnTaskStartTime + ExeTimeOfCrnTask;
                vector<set<double> > TemITL = ITL;
                UpdateITL(TemITL[RscForCrnTsk], CrnTaskStartTime, CrnTaskEndTime);
                //                chrom.StartTime[CrnTskId] = CrnTaskStartTime;
                chrom.EndTime[CrnTskId] = CrnTaskEndTime;
                for (int RscOfChildTsk : Tasks[NeedProcessChildTask].ElgRsc) {
                    double ChildReadyTime = 0, TransferDataSize = Tasks[NeedProcessChildTask].ExternalInputFileSizeSum;
                    for (int TemParentTask : Tasks[NeedProcessChildTask].parents) {
                        int RscOfTemParTsk = chrom.RscAlcLst[TemParentTask];
                        if (RscOfChildTsk != RscOfTemParTsk) {
                            TransferDataSize += ParChildTranFileSizeSum[TemParentTask][NeedProcessChildTask]; // / VALUE * 8 / (XY_MIN(Rscs[RscOfChildTsk].bw, Rscs[RscOfTemParTsk].bw));
                        }
                        if (ChildReadyTime + PrecisionValue < chrom.EndTime[TemParentTask]) {
                            ChildReadyTime = chrom.EndTime[TemParentTask];
                        }
                    }
                    double ChildExeTime = Tasks[NeedProcessChildTask].length / Rscs[RscOfChildTsk].pc + (TransferDataSize + Tasks[NeedProcessChildTask].OFileSizeSum) / VALUE * 8 / Rscs[RscOfChildTsk].bw;
                    double ChildStartTime = FindIdleTimeSlot(TemITL[RscOfChildTsk], ChildExeTime, ChildReadyTime);
                    double ChildEndTime = ChildStartTime + ChildExeTime;
                    if (FinalEndTimeOfChildTask > ChildEndTime + PrecisionValue) {
                        FinalRscForChildTask = RscOfChildTsk;
                        FinalStartTimeOfChildTask = ChildStartTime;
                        FinalEndTimeOfChildTask = ChildEndTime;
                        FinalStartTimeOfCrnTask = CrnTaskStartTime;
                        FinalEndTimeOfCrnTask = CrnTaskEndTime;
                        FinalRscForCrnTask = RscForCrnTsk;
                        ITL_CrnTaskScheduled = TemITL;
                    }
                }
            }
            chrom.RscAlcLst[CrnTskId] = FinalRscForCrnTask;
            chrom.StartTime[CrnTskId] = FinalStartTimeOfCrnTask;
            chrom.EndTime[CrnTskId] = FinalEndTimeOfCrnTask;
            chrom.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].first = XY_MIN(chrom.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].first, chrom.StartTime[CrnTskId]);
            chrom.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].second = XY_MAX(chrom.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].second, chrom.EndTime[CrnTskId]);
            chrom.TskSchLst[IndexCount] = NeedProcessChildTask;
            IndexCount++;
            chrom.RscAlcLst[NeedProcessChildTask] = FinalRscForChildTask;
            chrom.StartTime[NeedProcessChildTask] = FinalStartTimeOfChildTask;
            chrom.EndTime[NeedProcessChildTask] = FinalEndTimeOfChildTask;
            chrom.HtUseTime[Rscs[FinalRscForChildTask].Hostid].first = XY_MIN(chrom.HtUseTime[Rscs[FinalRscForChildTask].Hostid].first, chrom.StartTime[NeedProcessChildTask]);
            chrom.HtUseTime[Rscs[FinalRscForChildTask].Hostid].second = XY_MAX(chrom.HtUseTime[Rscs[FinalRscForChildTask].Hostid].second, chrom.EndTime[NeedProcessChildTask]);
            chrom.MakeSpan = XY_MAX(chrom.MakeSpan, chrom.EndTime[NeedProcessChildTask]);
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL[FinalRscForChildTask], FinalStartTimeOfChildTask, chrom.EndTime[NeedProcessChildTask]);
            for (int TskId : ReadyChildTaskSet) { 
                if (TskId != NeedProcessChildTask) {
                    RTI.push_back(TskId);
                    AvrRtSet[TskId] = ClcAvrReadyTime_S(TskId, chrom);
                }
            }
            for (int childId : Tasks[NeedProcessChildTask].children) {
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    RTI.push_back(childId);
                    AvrRtSet[childId] = ClcAvrReadyTime_S(childId, chrom);
                }
            }
        }
    }
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[chrom.TskSchLst[i]] = chrom.RscAlcLst[chrom.TskSchLst[i]] + Decimals[i];
    }
    CalculateEnergy1(chrom);
    return chrom;
}

double ClcAvrReadyTime_S(int TskId, chromosome& chrom) {
    double AvrRt = 0;
    for (int parent : Tasks[TskId].parents) {
        if (AvrRt + PrecisionValue < chrom.EndTime[parent]) {
            AvrRt = chrom.EndTime[parent];
        }
    }
    return AvrRt;
}

vector<double> GnrDecimalsByAscend() {
    vector<double> decimals(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        decimals[i] = (i + RandomDouble(0, 1)) / comConst.NumOfTsk;
    }
    return decimals;
}

chromosome GnrChr_IHEFT3_b_S(vector<double> Rank_b) {
    chromosome Chrom_Rank_b;
    IntChr(Chrom_Rank_b);
    IndexSortByValueOnDescend(Chrom_Rank_b.TskSchLst, Rank_b);
    IHEFT3_S(Chrom_Rank_b);
    CalculateEnergy1(Chrom_Rank_b);
    return Chrom_Rank_b;
}

chromosome GnrChr_IHEFT3_t_S(vector<double> Rank_t) {
    chromosome Chrom_Rank_t;
    IntChr(Chrom_Rank_t);
    IndexSortByValueOnAscend(Chrom_Rank_t.TskSchLst, Rank_t);
    IHEFT3_S(Chrom_Rank_t);
    CalculateEnergy1(Chrom_Rank_t);
    return Chrom_Rank_t;
}

double IHEFT3_S(chromosome& ch) {
    list <int> TemTskSchLst;
    TemTskSchLst.assign(ch.TskSchLst.begin(), ch.TskSchLst.end());
    ch.RscAlcLst.resize(comConst.NumOfTsk, -1);
    ch.TskSchLst.resize(comConst.NumOfTsk, -1);
    vector<set<double>> ITL;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    vector<int> upr(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
    }
    int IndexCount = 0;
    while (!TemTskSchLst.empty()) {
        int CrnTask = TemTskSchLst.front();
        ch.TskSchLst[IndexCount] = CrnTask;
        IndexCount++;
        TemTskSchLst.erase(TemTskSchLst.begin());
        int FinalRscForCrnTask = -1;
        double FinalStartTimeOfCrnTask = 0;
        double FinalEndTimeOfCrnTask = InfiniteValue;
        vector<int> NeedProcessChildTaskSet;
        for (int childId : Tasks[CrnTask].children) {
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                NeedProcessChildTaskSet.push_back(childId);
            }
        }
        if (NeedProcessChildTaskSet.empty()) {
            SelectRsc_EFT_S(ch, ITL, CrnTask, FinalRscForCrnTask, FinalStartTimeOfCrnTask, FinalEndTimeOfCrnTask);
            ch.StartTime[CrnTask] = FinalStartTimeOfCrnTask;
            ch.EndTime[CrnTask] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTask] = FinalRscForCrnTask;
            ch.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].first, ch.StartTime[CrnTask]);
            ch.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].second, ch.EndTime[CrnTask]);
            ch.MakeSpan = XY_MAX(ch.MakeSpan, ch.EndTime[CrnTask]);
            UpdateITL(ITL[FinalRscForCrnTask], FinalStartTimeOfCrnTask, FinalEndTimeOfCrnTask);
        }
        else {
            int FinalChildTask = -1, FinalRscForChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = InfiniteValue;
            vector<set<double> > ITL_CrnTaskScheduled;
            for (int CrnTaskRsc : Tasks[CrnTask].ElgRsc) {
                double ReadyTimeOfCrnTask = 0, TransferDataSize = Tasks[CrnTask].ExternalInputFileSizeSum;
                ch.RscAlcLst[CrnTask] = CrnTaskRsc;
                for (int ParentTask : Tasks[CrnTask].parents) {
                    int RscOfParentTask = ch.RscAlcLst[ParentTask];
                    if (CrnTaskRsc != RscOfParentTask) {
                        TransferDataSize += ParChildTranFileSizeSum[ParentTask][CrnTask];
                    }
                    if (ReadyTimeOfCrnTask + PrecisionValue < ch.EndTime[ParentTask]) {
                        ReadyTimeOfCrnTask = ch.EndTime[ParentTask];
                    }
                }
                double ExeTimeOfCrnTask = Tasks[CrnTask].length / Rscs[CrnTaskRsc].pc + (TransferDataSize + Tasks[CrnTask].OFileSizeSum) / VALUE * 8 / Rscs[CrnTaskRsc].bw;
                double CrnTaskStartTime = FindIdleTimeSlot(ITL[CrnTaskRsc], ExeTimeOfCrnTask, ReadyTimeOfCrnTask);
                double CrnTaskEndTime = CrnTaskStartTime + ExeTimeOfCrnTask;
                vector<set<double> > TemITL = ITL;
                UpdateITL(TemITL[CrnTaskRsc], CrnTaskStartTime, CrnTaskEndTime);
                //                ch.StartTime[CrnTask] = CrnTaskStartTime;
                ch.EndTime[CrnTask] = CrnTaskEndTime;
                for (int TemChildTask : NeedProcessChildTaskSet) {
                    for (int TemChildRsc : Tasks[TemChildTask].ElgRsc) {
                        double TemChildReadyTime = 0, TransferDataSize = Tasks[TemChildTask].ExternalInputFileSizeSum;
                        for (int TemParentTask : Tasks[TemChildTask].parents) {
                            int TemParRsc = ch.RscAlcLst[TemParentTask];
                            if (TemChildRsc != TemParRsc) {
                                TransferDataSize += ParChildTranFileSizeSum[TemParentTask][TemChildTask];
                            }
                            if (TemChildReadyTime + PrecisionValue < ch.EndTime[TemParentTask]) {
                                TemChildReadyTime = ch.EndTime[TemParentTask];
                            }
                        }
                        double TemChildExeTime = Tasks[TemChildTask].length / Rscs[TemChildRsc].pc + (TransferDataSize + Tasks[TemChildTask].OFileSizeSum) / VALUE * 8 / Rscs[TemChildRsc].bw;
                        double TemChildStartTime = FindIdleTimeSlot(TemITL[TemChildRsc], TemChildExeTime, TemChildReadyTime);
                        double TemChildEndTime = TemChildStartTime + TemChildExeTime;
                        if (FinalEndTimeOfChildTask > TemChildEndTime + PrecisionValue) {
                            FinalChildTask = TemChildTask;
                            FinalRscForChildTask = TemChildRsc;
                            FinalStartTimeOfChildTask = TemChildStartTime;
                            FinalEndTimeOfChildTask = TemChildEndTime;
                            FinalStartTimeOfCrnTask = CrnTaskStartTime;
                            FinalEndTimeOfCrnTask = CrnTaskEndTime;
                            FinalRscForCrnTask = CrnTaskRsc;
                            ITL_CrnTaskScheduled = TemITL;
                        }
                    }
                }
            }
            ch.StartTime[CrnTask] = FinalStartTimeOfCrnTask;
            ch.EndTime[CrnTask] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTask] = FinalRscForCrnTask;
            ch.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].first, ch.StartTime[CrnTask]);
            ch.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[FinalRscForCrnTask].Hostid].second, ch.EndTime[CrnTask]);
            ch.TskSchLst[IndexCount] = FinalChildTask;
            IndexCount++;
            ch.RscAlcLst[FinalChildTask] = FinalRscForChildTask;
            ch.StartTime[FinalChildTask] = FinalStartTimeOfChildTask;
            ch.EndTime[FinalChildTask] = FinalEndTimeOfChildTask;
            ch.HtUseTime[Rscs[FinalRscForChildTask].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[FinalRscForChildTask].Hostid].first, ch.StartTime[FinalChildTask]);
            ch.HtUseTime[Rscs[FinalRscForChildTask].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[FinalRscForChildTask].Hostid].second, ch.EndTime[FinalChildTask]);
            ch.MakeSpan = XY_MAX(ch.MakeSpan, ch.EndTime[FinalChildTask]);
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL[FinalRscForChildTask], FinalStartTimeOfChildTask, ch.EndTime[FinalChildTask]);
            TemTskSchLst.erase(find(TemTskSchLst.begin(), TemTskSchLst.end(), FinalChildTask));
            for (int childId : Tasks[FinalChildTask].children) {
                upr[childId] = upr[childId] - 1;
            }
        }
    }
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.Code_RK[ch.TskSchLst[i]] = ch.RscAlcLst[ch.TskSchLst[i]] + Decimals[i];
    }
    return ch.MakeSpan;
}

chromosome GnrChr_Lvl_Ran() {
    chromosome chrom;
    IntChr(chrom);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()] + (LevelIdOfTask[i] + (rand() % 10000) / 10000.0) / TskLstInLvl.size();
    }
    return chrom;
}

chromosome GnrChr_HMEC_S(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrML_Evl_MEC_S(TemChrom);
    return TemChrom;
}



chromosome GetTskSchLst(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    return TemChrom;
}


vector<int> GetOriginSequence(vector<double> seq) {
    vector<int> ind(comConst.NumOfTsk);
    vector<int> TempSeq(comConst.NumOfTsk);

    IndexSortByValueOnAscend(ind,seq);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TempSeq[i] = ind[comConst.NumOfTsk - i - 1];
    }
    return TempSeq;
}



double GnrML_Evl_MEC_S_New(chromosome& ch) {

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.RscAlcLst[i] = -1;
    }

    vector<set<double> > ITL;    
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }

    ch.MakeSpan = 0;
    ch.EnergyConsumption = 0;

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1;
        int TaskIndex = ch.TskSchLst[i];
        double FinalEC = 9999999999, FinalEndTime = 9999999999, FinalStartTime = 0;
        for (int v : Tasks[TaskIndex].ElgRsc) {
            double ReadyTime = 0;
            double TransferDataSize = Tasks[TaskIndex].ExternalInputFileSizeSum;

            for (int ParentIndex : Tasks[TaskIndex].parents) {
                int ParentRscIndex = ch.RscAlcLst[ParentIndex];

                if (v != ParentRscIndex) {
                    TransferDataSize += ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                }

                if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]) {
                    ReadyTime = ch.EndTime[ParentIndex];
                }
            }

            double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc + (TransferDataSize + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[v].bw;
            ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[v], ExeTime, ReadyTime);
            ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExeTime;
            ch.RscAlcLst[TaskIndex] = v;
            double CurEC = CalculateECByDelta2(ch, i);

            if (CurEC + PrecisionValue < FinalEC) {
                FinalStartTime = ch.StartTime[TaskIndex];
                FinalEndTime = ch.EndTime[TaskIndex];
                FinalEC = CurEC;
                RscIndex = v;
            }
        }

        ch.EnergyConsumption = FinalEC;
        ch.StartTime[TaskIndex] = FinalStartTime;
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        ch.HtUseTime[Rscs[RscIndex].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[RscIndex].Hostid].first, ch.StartTime[TaskIndex]);
        ch.HtUseTime[Rscs[RscIndex].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[RscIndex].Hostid].second, ch.EndTime[TaskIndex]);
        ch.MakeSpan = XY_MAX(ch.MakeSpan, FinalEndTime);
        UpdateITL(ITL[RscIndex], FinalStartTime, FinalEndTime);

    }

    return ch.EnergyConsumption;
}




double GnrML_Evl_MEC_S(chromosome& ch) {

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.RscAlcLst[i] = -1;
    }   

    vector<set<double> > ITL;    
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }

    ch.MakeSpan = 0;
    ch.EnergyConsumption = 0;


    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1;
        int TaskIndex = ch.TskSchLst[i];
        double FinalEC = 9999999999, FinalEndTime = 9999999999, FinalStartTime = 0;
        for (int v : Tasks[TaskIndex].ElgRsc) {
            double ReadyTime = 0;
            double TransferDataSize = Tasks[TaskIndex].ExternalInputFileSizeSum;

            for (int ParentIndex : Tasks[TaskIndex].parents) {
                int ParentRscIndex = ch.RscAlcLst[ParentIndex];

                if (v != ParentRscIndex) {
                    TransferDataSize += ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                }

                if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]) {
                    ReadyTime = ch.EndTime[ParentIndex];
                }
            }

            double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc + (TransferDataSize + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[v].bw;
            ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[v], ExeTime, ReadyTime);
            ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExeTime;
            ch.RscAlcLst[TaskIndex] = v;
            double CurEC = CalculateECByDelta2(ch, i);

            if (CurEC + PrecisionValue < FinalEC) {
                FinalStartTime = ch.StartTime[TaskIndex];
                FinalEndTime = ch.EndTime[TaskIndex];
                FinalEC = CurEC;
                RscIndex = v;
            }
        }

        ch.EnergyConsumption = FinalEC;
        ch.StartTime[TaskIndex] = FinalStartTime;
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        ch.HtUseTime[Rscs[RscIndex].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[RscIndex].Hostid].first, ch.StartTime[TaskIndex]);
        ch.HtUseTime[Rscs[RscIndex].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[RscIndex].Hostid].second, ch.EndTime[TaskIndex]);
        ch.MakeSpan = XY_MAX(ch.MakeSpan, FinalEndTime);
        UpdateITL(ITL[RscIndex], FinalStartTime, FinalEndTime);

    }

    return ch.EnergyConsumption;
}





double GnrML_Evl_MEC_S_QPHH(chromosome& ch) {

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.RscAlcLst[i] = -1;
    }

    vector<set<double> > ITL;    
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    vector<int> RscUseNum(comConst.NumOfRsc, 0);
    ch.MakeSpan = 0;
    ch.EnergyConsumption = 0;
    double theta = 1.0;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1;
        int TaskIndex = ch.TskSchLst[i];
        double FinalEC = 9999999999, FinalEndTime = 9999999999, FinalStartTime = 0;

        double SelectTheta = ((double)rand() / (RAND_MAX)) / 1.0;
        if (SelectTheta <= theta){
            for (int v : Tasks[TaskIndex].ElgRsc) {
                double ReadyTime = 0;
                double TransferDataSize = Tasks[TaskIndex].ExternalInputFileSizeSum;

                for (int ParentIndex : Tasks[TaskIndex].parents) {
                    int ParentRscIndex = ch.RscAlcLst[ParentIndex];

                    if (v != ParentRscIndex) {
                        TransferDataSize += ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                    }

                    if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]) {
                        ReadyTime = ch.EndTime[ParentIndex];
                    }
                }

                double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc + (TransferDataSize + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[v].bw;
                ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[v], ExeTime, ReadyTime);
                ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExeTime;
                ch.RscAlcLst[TaskIndex] = v;
                double CurEC = CalculateECByDelta2(ch, i);

                if (CurEC + PrecisionValue < FinalEC) {
                    FinalStartTime = ch.StartTime[TaskIndex];
                    FinalEndTime = ch.EndTime[TaskIndex];
                    FinalEC = CurEC;
                    RscIndex = v;
                }
            }
            RscUseNum[RscIndex]++;
        }
        else {
            for (int v : Tasks[TaskIndex].ElgRsc) { 
                if (RscUseNum[v] != 0) {
                    double ReadyTime = 0;
                    double TransferDataSize = Tasks[TaskIndex].ExternalInputFileSizeSum;

                    for (int ParentIndex : Tasks[TaskIndex].parents) {
                        int ParentRscIndex = ch.RscAlcLst[ParentIndex];

                        if (v != ParentRscIndex) {
                            TransferDataSize += ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                        }
                        if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]) {
                            ReadyTime = ch.EndTime[ParentIndex];
                        }
                    }

                    double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc + (TransferDataSize + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[v].bw;
                    ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[v], ExeTime, ReadyTime);
                    ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExeTime;
                    ch.RscAlcLst[TaskIndex] = v;
                    double CurEC = CalculateECByDelta2(ch, i);

                    if (CurEC + PrecisionValue < FinalEC) {
                        FinalStartTime = ch.StartTime[TaskIndex];
                        FinalEndTime = ch.EndTime[TaskIndex];
                        FinalEC = CurEC;
                        RscIndex = v; 
                    }
                }
                
            }
            RscUseNum[RscIndex]++;
        }

        ch.EnergyConsumption = FinalEC;
        ch.StartTime[TaskIndex] = FinalStartTime;
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        ch.HtUseTime[Rscs[RscIndex].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[RscIndex].Hostid].first, ch.StartTime[TaskIndex]);
        ch.HtUseTime[Rscs[RscIndex].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[RscIndex].Hostid].second, ch.EndTime[TaskIndex]);
        ch.MakeSpan = XY_MAX(ch.MakeSpan, FinalEndTime);
        UpdateITL(ITL[RscIndex], FinalStartTime, FinalEndTime);

        theta = theta * 0.95;
    }
    
    return ch.EnergyConsumption;
}


double DcdEvl_S(chromosome& ch) {
    for (int k = 0; k < HstSet.size(); ++k) {
        ch.HtUseTime[k] = { 9999999999,-1 };
    }
    ch.MakeSpan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);  a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i], RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the task and resource (Rsc) allocated to this task
        double ReadyTime = 0;
        double TransferSize = Tasks[TaskIndex].ExternalInputFileSizeSum; //the size of input files that need to obtain from share database;
        for (int ParentTask : Tasks[TaskIndex].parents) {
            if (RscIndex != ch.RscAlcLst[ParentTask]) {
                TransferSize = TransferSize + ParChildTranFileSizeSum[ParentTask][TaskIndex];
            }
        }
        if (ch.IsFrw) {   //forward-loading
            for (int ParentTask : Tasks[TaskIndex].parents) {
                if (ReadyTime + PrecisionValue < ch.EndTime[ParentTask]) {
                    ReadyTime = ch.EndTime[ParentTask];
                }
            }
        }
        else {         //backward-loading
            for (int ChildTask : Tasks[TaskIndex].children) {
                if (ReadyTime + PrecisionValue < ch.EndTime[ChildTask]) {
                    ReadyTime = ch.EndTime[ChildTask];
                }
            }
        }
        double ExecutionTime = (TransferSize + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[RscIndex].bw + Tasks[TaskIndex].length / Rscs[RscIndex].pc;
        ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[RscIndex], ExecutionTime, ReadyTime);
        ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExecutionTime;
        ch.HtUseTime[Rscs[RscIndex].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[RscIndex].Hostid].first, ch.StartTime[TaskIndex]);
        ch.HtUseTime[Rscs[RscIndex].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[RscIndex].Hostid].second, ch.EndTime[TaskIndex]);
        if (ch.MakeSpan + PrecisionValue < ch.EndTime[TaskIndex]) {
            ch.MakeSpan = ch.EndTime[TaskIndex];
        }
        //{update ITL}
        UpdateITL(ITL[RscIndex], ch.StartTime[TaskIndex], ch.EndTime[TaskIndex]);
    }
    return ch.MakeSpan;
}

void AdpDcd_S(chromosome& chrom, double& CurTime, double& TotalTime) {
    vector<double> SP(3);
    SP[0] = pow(CurTime / TotalTime, Parameter_ADBRKGA.alpha);
    SP[1] = Parameter_ADBRKGA.beta * (1 - SP[0]);
    SP[2] = (1 - Parameter_ADBRKGA.beta) * (1 - SP[0]);
    double RandNum = double(rand() % 1000) / 1000;
    if (RandNum < SP[0]) {
        NrmDcd_S(chrom);
    }
    if (RandNum >= SP[0] && RandNum < (SP[0] + SP[1])) {
        HrsDcd_EFT_S(chrom);
        //        GnrRscAlcTskSchLstFromCode_RK(chrom);
        //        GnrML_Evl_MEC_S(chrom);
        //        ModifyRscAlcLstByCode_RK(chrom);
    }
    if (RandNum >= (SP[0] + SP[1])) {
        HrsDcd_CTP_S(chrom);
    }
}

double NrmDcd_S(chromosome& ch) {
    vector<int > upr(comConst.NumOfTsk, -1);
    list<int> RTI;
    if (ch.IsFrw)
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i] = Tasks[i].parents.size();
            if (upr[i] == 0)  RTI.push_back(i);
        }
    else
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i] = Tasks[i].children.size();
            if (upr[i] == 0)  RTI.push_back(i);
        }
    //generate resource allocation list and task scheduling order list
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.RscAlcLst[i] = floor(ch.Code_RK[i]);
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            int decimal = ch.Code_RK[*lit] - floor(ch.Code_RK[*lit]);
            if (decimal < tmp) {
                tmp = decimal; pit = lit; 
            }
        }
        ch.TskSchLst[i] = *pit;
        RTI.erase(pit);

        if (ch.IsFrw)
            for (int childId : Tasks[ch.TskSchLst[i]].children) {
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0)   RTI.push_back(childId);
            }
        else
            for (int parentId : Tasks[ch.TskSchLst[i]].parents) {
                upr[parentId] = upr[parentId] - 1;
                if (upr[parentId] == 0)  RTI.push_back(parentId);
            }
    }
    DcdEvl_S(ch);
    return ch.MakeSpan;
}

double HrsDcd_CTP_S(chromosome& ch) {
    vector<double> ExeTime(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<int> ind(comConst.NumOfTsk);
    vector<vector<double>> TransferTime(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    //{calculate the transfer time between tasks when resource(Rsc) allocation has been determined}
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = floor(ch.Code_RK[i]);
        double TransferDataSize = Tasks[i].ExternalInputFileSizeSum;
        for (int parent : Tasks[i].parents) {
            int ParRsc = floor(ch.Code_RK[parent]);
            if (ParRsc != RscIndex) {
                TransferDataSize = TransferDataSize + ParChildTranFileSizeSum[parent][i];
            }
        }
        //+ Tasks[i].ExternalInputFileSizeSum -revised-xy in 2022.02.18
        ExeTime[i] = Tasks[i].length / Rscs[RscIndex].pc + (TransferDataSize + Tasks[i].OFileSizeSum + Tasks[i].ExternalInputFileSizeSum) / VALUE * 8 / Rscs[RscIndex].bw;
    }
    Calculate_Rank_b_S(Rank_b, ExeTime);
    IndexSortByValueOnAscend(ind, Rank_b);
    vector<double> Decimals = GnrDecimalsByAscend(); 
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.TskSchLst[i] = TaskIndex;
        ch.Code_RK[TaskIndex] = floor(ch.Code_RK[TaskIndex]) + Decimals[i];
    }
    NrmDcd_S(ch);
    return ch.MakeSpan;
}

void GnrCode_RK(chromosome& chrom) {
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[chrom.TskSchLst[i]] = chrom.RscAlcLst[chrom.TskSchLst[i]] + Decimals[i];
    }
}


void ModifyRscAlcLstByCode_RK(chromosome& chrom) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[i] = chrom.RscAlcLst[i] + chrom.Code_RK[i] - floor(chrom.Code_RK[i]);
    }
}

void GnrRscAlcTskSchLstFromCode_RK(chromosome& chrom) {
    vector<int > upr(comConst.NumOfTsk, -1);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }
    //generate resource allocation list and task scheduling order list
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.RscAlcLst[i] = floor(chrom.Code_RK[i]);
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            int decimal = chrom.Code_RK[*lit] - floor(chrom.Code_RK[*lit]);
            if (decimal + PrecisionValue < tmp) {
                tmp = decimal; pit = lit; 
            }
        }
        chrom.TskSchLst[i] = *pit;
        RTI.erase(pit);
        //更新RTI;
        for (int childId : Tasks[chrom.TskSchLst[i]].children) {
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0)   RTI.push_back(childId);
        }
    }
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk, -1);
    chrom.Code_RK.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk, -1);
    chrom.EndTime.resize(comConst.NumOfTsk);
    chrom.StartTime.resize(comConst.NumOfTsk);
    chrom.Code_TD.resize(comConst.NumOfRsc);
    chrom.TskSchPart.resize(comConst.NumOfTsk);
    chrom.RscAlcPart.resize(comConst.NumOfTsk);
    chrom.VTskSchPart.resize(comConst.NumOfTsk, 0);
    chrom.VRscAlcPart.resize(comConst.NumOfTsk, 0);
    chrom.HtUseTime.resize(HstSet.size(), { 9999999999,-1 });
    chrom.IsFrw = true;
}

chromosome GnrChr_HEFT_Baseline_S() {
    chromosome TemChrom;
    IntChr(TemChrom);
    int ScheduleOrder = 0;
    for (int j = 0; j < TskLstInLvl.size(); ++j) {
        if (TskLstInLvl[j].size() < 2) {
            TemChrom.TskSchLst[ScheduleOrder++] = TskLstInLvl[j][0];
            continue;
        }
        vector<int> SonTaskNum;
        for (int i = 0; i < TskLstInLvl[j].size(); ++i)
            SonTaskNum.push_back(Tasks[TskLstInLvl[j][i]].children.size());

        vector<int> ind(TskLstInLvl[j].size());
        IndexSortByValueOnAscend(ind, SonTaskNum);
        for (int i = TskLstInLvl[j].size() - 1; i >= 0; i--) {
            TemChrom.TskSchLst[ScheduleOrder++] = TskLstInLvl[j][ind[i]];
        }
    }
    GnrML_Evl_EFT_S(TemChrom);
    CalculateEnergy1(TemChrom);
    return TemChrom;
}

void UpdateITL(set<double>& ITLofRscId, double& StartTime, double& EndTime) {
    if (ITLofRscId.find(StartTime) != ITLofRscId.end()) {
        ITLofRscId.erase(StartTime);
    }
    else {
        ITLofRscId.insert(StartTime);
    }
    if (ITLofRscId.find(EndTime) != ITLofRscId.end()) {
        ITLofRscId.erase(EndTime);
    }
    else {
        ITLofRscId.insert(EndTime);
    }
}


double GnrML_Evl_EFT_S(chromosome& ch) {

    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;

    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1, TaskId = ch.TskSchLst[i];
        double FinalEndTime = 9999999999, FinalStartTime = 0;

        SelectRsc_EFT_S(ch, ITL, TaskId, RscId, FinalStartTime, FinalEndTime);  //Find the resource that can finish the task earliest
        
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;

        ch.HtUseTime[Rscs[RscId].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[RscId].Hostid].first, ch.StartTime[TaskId]);
        ch.HtUseTime[Rscs[RscId].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[RscId].Hostid].second, ch.EndTime[TaskId]);
        
        UpdateITL(ITL[RscId], FinalStartTime, FinalEndTime);              //update ITL
        
        makespan = XY_MAX(makespan, FinalEndTime);
    }

    ch.MakeSpan = makespan;
    return makespan;
}

double HrsDcd_EFT_S(chromosome& ch) {
    double makespan = 0;
    for (int k = 0; k < HstSet.size(); ++k) {
        ch.HtUseTime[k] = { 9999999999,-1 };
    }
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    vector<int > upr(comConst.NumOfTsk, 0.0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1;
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            double decimal = ch.Code_RK[*lit] - floor(ch.Code_RK[*lit]);
            if (decimal < tmp) {
                tmp = decimal; pit = lit; 
            }
        }
        ch.TskSchLst[i] = *pit;
        RTI.erase(pit);
        double FinalEndTime = InfiniteValue;
        double FinalStartTime = 0;
        SelectRsc_EFT_S(ch, ITL, ch.TskSchLst[i], RscId, FinalStartTime, FinalEndTime);
        ch.StartTime[ch.TskSchLst[i]] = FinalStartTime;
        ch.EndTime[ch.TskSchLst[i]] = FinalEndTime;
        ch.HtUseTime[Rscs[RscId].Hostid].first = XY_MIN(ch.HtUseTime[Rscs[RscId].Hostid].first, ch.StartTime[ch.TskSchLst[i]]);
        ch.HtUseTime[Rscs[RscId].Hostid].second = XY_MAX(ch.HtUseTime[Rscs[RscId].Hostid].second, ch.EndTime[ch.TskSchLst[i]]);
        ch.Code_RK[ch.TskSchLst[i]] = RscId + tmp;
        ch.RscAlcLst[ch.TskSchLst[i]] = RscId;
        UpdateITL(ITL[RscId], FinalStartTime, FinalEndTime); //{update ITL}
        makespan = XY_MAX(makespan, FinalEndTime);
        for (int ChildId : Tasks[ch.TskSchLst[i]].children) {
            upr[ChildId] = upr[ChildId] - 1;
            if (upr[ChildId] == 0)   RTI.push_back(ChildId);
        }
    }
    ch.MakeSpan = makespan;
    return makespan;
}  //ADBRKGA

void SelectRsc_EFT_S(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime) {
    for (int RscIdOfCrnTsk : Tasks[TaskIndex].ElgRsc) {
        double ReadyTime = 0;
        double TransferData = Tasks[TaskIndex].ExternalInputFileSizeSum;
        for (int ParentIndex : Tasks[TaskIndex].parents) { //calculate the ready time and transfer data of the task
            int RscIdOfPrnTsk = ch.RscAlcLst[ParentIndex];
            if (RscIdOfCrnTsk != RscIdOfPrnTsk) {
                TransferData = TransferData + ParChildTranFileSizeSum[ParentIndex][TaskIndex];
            }
            if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]) {
                ReadyTime = ch.EndTime[ParentIndex];
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[RscIdOfCrnTsk].pc + (TransferData + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[RscIdOfCrnTsk].bw;
        double StartTime = FindIdleTimeSlot(ITL[RscIdOfCrnTsk], ExeTime, ReadyTime); //Find an idle time-slot as early as possible from ITL
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + PrecisionValue < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscIndex = RscIdOfCrnTsk;
        }
    }
}

double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime) {
    set<double>::iterator pre = ITLofRscId.begin();
    set<double>::iterator post = ITLofRscId.begin();
    ++post;
    while (post != ITLofRscId.end()) {
        if ((*post - *pre) > ExeTime - PrecisionValue && ReadyTime - PrecisionValue < (*post) - ExeTime) {
            return  XY_MAX(*pre, ReadyTime);
        }
        else {
            ++pre; ++pre; ++post; ++post;
        }
    }
}

//{generate a topological sort randomly}
vector<int> GnrSS_TS() {
    vector<int> SS;
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0) {
            RTI.push_back(i);
        }
    }

    while (!RTI.empty()) {
        int RandVec = rand() % RTI.size();
        int v = RTI[RandVec];
        vector<int>::iterator iter = RTI.begin() + RandVec;
        RTI.erase(iter);
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            --upr[Tasks[v].children[i]];
            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
        }
        SS.push_back(v);
    }
    return SS;
}

vector<int> GnrSS_Lvl() {
    vector<int> ch;
    vector<vector<int>> tem = TskLstInLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        random_shuffle(tem[i].begin(), tem[i].end());   //arrange the tasks in each level
        for (int j = 0; j < tem[i].size(); ++j) {
            ch.push_back(tem[i][j]);
        }
    }
    return ch;
}

void InitProModelOfResAlc(vector<vector<double> >& PMR) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int j : Tasks[i].ElgRsc) {
            PMR[i][j] = 1.0 / Tasks[i].ElgRsc.size();
        }
    }
}


void InitProModelOfTskSch(vector<vector<double>>& PMS, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants) {
    vector<int> STS(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int k = NumOfAncestors[i]; k < NumOfNonDescendants[i]; ++k) {
            PMS[i][k] = 1;
            ++STS[k];
        }
    }//PMS[i][k] represents the probability that the k-th scheduled task is task i

    for (int k = 0; k < comConst.NumOfTsk; ++k) {
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            PMS[i][k] = PMS[i][k] / STS[k];
        }
    }
}

chromosome GnrTskLstOfChr(vector<vector<double> >& PMS, vector<double>& eta_TSO) {
    chromosome chrom;
    IntChr(chrom);
    vector<int > upr(comConst.NumOfTsk, 0);
    list<int> RTI;

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }


    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double sum = 0;

        for (int k : RTI) {
            sum += PMS[k][i] * eta_TSO[k];
        }

        vector<double> SltProb(comConst.NumOfTsk);

        for (int k : RTI) {
            SltProb[k] = PMS[k][i] * eta_TSO[k] / sum;
        }
        double ProbSum = 0;
        double rnd = double(rand() % 100) / 100;
        for (int k : RTI) {
            ProbSum += SltProb[k];
            if (rnd + PrecisionValue < ProbSum) {
                chrom.TskSchLst[i] = k;
                break;
            }
        }

        RTI.erase(find(RTI.begin(), RTI.end(), chrom.TskSchLst[i]));

        for (int k = 0; k < Tasks[chrom.TskSchLst[i]].children.size(); ++k) {
            upr[Tasks[chrom.TskSchLst[i]].children[k]]--;
            if (upr[Tasks[chrom.TskSchLst[i]].children[k]] == 0) {
                RTI.push_back(Tasks[chrom.TskSchLst[i]].children[k]);
            }
        }
    }
    return chrom;
}


void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR) {
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double rnd = double(rand() % 100) / 100;
        double sum = 0;
        for (int j : Tasks[i].ElgRsc) {
            sum += PMR[i][j];
            if (rnd < sum) {
                chrom.RscAlcLst[i] = j;
                break;
            }
        }
    }
}

chromosome GnrPrtByRank_Rnd_S(vector<double>& Rank) {
    chromosome chrom;
    IntChr(chrom);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.RscAlcPart[i] = RandomDouble2(0, comConst.NumOfRsc - 1); //RandomDouble2(0,comConst.NumOfRsc);//rand() % comConst.NumOfRsc + rand() % 1000 / 1000.0 - 0.5;
    }
    RepairMapAndGnrRscAlcLst(chrom); //GnrRscAlcLst(chrom); //
    chrom.TskSchPart = Rank;
    RepairPriorityAndGnrSchOrd(chrom);
    DcdEvl_S(chrom);
    CalculateEnergy1(chrom);
    return chrom;
}

chromosome GnrPrtByRank_EFT_S(vector<double>& Rank) {
    chromosome chrom;
    IntChr(chrom);
    chrom.TskSchPart = Rank;
    RepairPriorityAndGnrSchOrd(chrom);
    GnrML_Evl_EFT_S(chrom);
    CalculateEnergy1(chrom);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.RscAlcPart[i] = chrom.RscAlcLst[i] - 0.5 + (rand() % 10000) / 10000.0;
    }
    return chrom;
}

void RepairPriorityAndGnrSchOrd(chromosome& chrom) {
    vector<int> V, Q;
    vector<int> N(comConst.NumOfTsk, -1);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        N[i] = round(chrom.TskSchPart[i]);
    }
    vector<int> upr(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0) {
            Q.push_back(i);
        }
    }
    int MaxV = -1;
    while (V.size() != comConst.NumOfTsk) {
        for (int i = 0; i < Q.size(); ++i) {
            int TaskId = Q[i];
            int MaxP = -1;
            for (int i1 = 0; i1 < Tasks[TaskId].parents.size(); ++i1) {
                if (MaxP < N[Tasks[TaskId].parents[i1]]) {
                    MaxP = N[Tasks[TaskId].parents[i1]];
                }
            }
            if (N[TaskId] <= MaxP) {
                N[TaskId] = MaxP + 1;
            }
            for (int i1 = 0; i1 < V.size(); ++i1) {
                if (N[TaskId] == N[V[i1]]) {
                    N[TaskId] = MaxV + 1;
                    MaxV += 1;
                    break;
                }
            }
            MaxV = XY_MAX(N[TaskId], MaxV);
            V.push_back(TaskId);
        }
        vector<int> TemQ;
        for (int i = 0; i < Q.size(); ++i) {
            int taskId = Q[i];
            for (int childId : Tasks[taskId].children) {
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    TemQ.push_back(childId);
                }
            }
        }
        Q = TemQ;
    }
    //    chrom.TskSchPart = N;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.TskSchPart[i] = N[i];
    }
    IndexSortByValueOnAscend(chrom.TskSchLst, N);
}

void RepairMapAndGnrRscAlcLst(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = round(ch.RscAlcPart[i]);
        if (RscId < Tasks[i].ElgRsc[0]) { 
            ch.RscAlcPart[i] = Tasks[i].ElgRsc[0];
            ch.RscAlcLst[i] = Tasks[i].ElgRsc[0];
            continue;
        }
        if (RscId > Tasks[i].ElgRsc[Tasks[i].ElgRsc.size() - 1]) { 
            ch.RscAlcPart[i] = Tasks[i].ElgRsc[Tasks[i].ElgRsc.size() - 1];
            ch.RscAlcLst[i] = Tasks[i].ElgRsc[Tasks[i].ElgRsc.size() - 1];
            continue;
        }
        if (find(Tasks[i].ElgRsc.begin(), Tasks[i].ElgRsc.end(), RscId) == Tasks[i].ElgRsc.end()) {
            if (Tasks[i].ElgRsc.size() == 1) {
                ch.RscAlcPart[i] = Tasks[i].ElgRsc[0];
                ch.RscAlcLst[i] = Tasks[i].ElgRsc[0];
            }
            else {
                int TemRscId = FindNearestRscId(i, ch.RscAlcPart[i]);
                ch.RscAlcPart[i] = TemRscId;
                ch.RscAlcLst[i] = TemRscId;
            }
            continue;
        }
        ch.RscAlcLst[i] = RscId;
    }
}

int FindNearestRscId(int TaskId, double value) {
    for (int j = 0; j < Tasks[TaskId].ElgRsc.size() - 1; ++j) {
        if (Tasks[TaskId].ElgRsc[j] < value && value < Tasks[TaskId].ElgRsc[j + 1]) {
            if (Tasks[TaskId].ElgRsc[j + 1] - value < value - Tasks[TaskId].ElgRsc[j]) {
                return Tasks[TaskId].ElgRsc[j + 1];
            }
            else {
                return Tasks[TaskId].ElgRsc[j];
            }
        }
    }
}

void UpdateParticle(chromosome& ch, chromosome& Pbest, chromosome& Gbest, double& runtime, double& SchTime) {
    Parameter_HPSO.InertiaWeight = 0.1 * (1 - (runtime / SchTime)) + 0.9;
    Parameter_HPSO.c1 = 2 * (1 - (runtime / SchTime));
    Parameter_HPSO.c2 = 2 * (runtime / SchTime);
    double r1 = RandomDouble(0, 1);
    double r2 = RandomDouble(0, 1);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.VTskSchPart[i] = Parameter_HPSO.InertiaWeight * ch.VTskSchPart[i] + Parameter_HPSO.c1 * r1 * (Pbest.TskSchPart[i] - ch.TskSchPart[i])
            + Parameter_HPSO.c2 * r2 * (Gbest.TskSchPart[i] - ch.TskSchPart[i]);
        ch.TskSchPart[i] += ch.VTskSchPart[i];

        ch.VRscAlcPart[i] = Parameter_HPSO.InertiaWeight * ch.VRscAlcPart[i] + Parameter_HPSO.c1 * r1 * (Pbest.RscAlcPart[i] - ch.RscAlcPart[i])
            + Parameter_HPSO.c2 * r2 * (Gbest.RscAlcPart[i] - ch.RscAlcPart[i]);
        ch.RscAlcPart[i] += ch.VRscAlcPart[i];
    }
    RepairMapAndGnrRscAlcLst(ch);   //GnrRscAlcLst(ch); //
    RepairPriorityAndGnrSchOrd(ch);
}

double IFBD_S(chromosome& ch) {
    chromosome NewChrom = ch;
    chromosome OldChrom;
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSortByValueOnAscend(ind, OldChrom.EndTime);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.TskSchLst[comConst.NumOfTsk - 1 - i] = ind[i];
        }
        NewChrom.IsFrw = !(OldChrom.IsFrw);
        DcdEvl_S(NewChrom);
        CalculateEnergy1(NewChrom);
    } while (NewChrom.EnergyConsumption + PrecisionValue < OldChrom.EnergyConsumption);
    if (fabs(NewChrom.EnergyConsumption - OldChrom.EnergyConsumption) < PrecisionValue && NewChrom.IsFrw) { 
        ch = NewChrom;
    }
    else { 
        ch = OldChrom;
    }
    return ch.EnergyConsumption;
}

double IFBS_S(chromosome& ch) {
    chromosome NewChrom = ch;
    chromosome OldChrom;
    vector<double> Decimals = GnrDecimalsByAscend();
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSortByValueOnAscend(ind, OldChrom.EndTime);
        //vector<double> Decimals = GnrDecimalsByAscend();
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.Code_RK[ind[comConst.NumOfTsk - 1 - i]] = floor(NewChrom.Code_RK[ind[comConst.NumOfTsk - 1 - i]]) + Decimals[i];
        }
        NewChrom.IsFrw = !(OldChrom.IsFrw);
        NrmDcd_S(NewChrom);
        CalculateEnergy1(NewChrom);
    } while (NewChrom.EnergyConsumption + PrecisionValue < OldChrom.EnergyConsumption);
    if (OldChrom.IsFrw) {
        ch = OldChrom;
    }
    else {
        ch = NewChrom;
    }
    return ch.EnergyConsumption;
} //ADBRKGA

//{ Load Balancing with Communication Reduction Improvement (LBCRI)}
void LBCA_S(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Id(comConst.NumOfRsc, 0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    vector<double> TransferDataSize;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        double TransferDataSize = Tasks[i].ExternalInputFileSizeSum;
        for (int Parent : Tasks[i].parents) {
            int PrnRscId = ch.RscAlcLst[Parent];
            if (RscIndex != PrnRscId) {
                TransferDataSize += ParChildTranFileSizeSum[Parent][i];
            }
        }
        Id[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc + (TransferDataSize + Tasks[i].OFileSizeSum) / VALUE * 8 / Rscs[RscIndex].bw;
        TSK[RscIndex].push_back(i);
    }
    vector<int> ind(comConst.NumOfRsc);
    IndexSortByValueOnAscend(ind, Id);         //sorting according to loads
    int RscWithMinLd = ind[0];          //find out the resource (Rsc) with the lowest load;
    set<int> ST;
    if (fabs(Id[RscWithMinLd]) < PrecisionValue) {
        ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
    }
    else {
        //traverse the tasks allocated to the resource with lowest load and add their parents and children to set ST
        for (int TaskIndex : TSK[RscWithMinLd]) {
            ST.insert(Tasks[TaskIndex].children.begin(), Tasks[TaskIndex].children.end());
            ST.insert(Tasks[TaskIndex].parents.begin(), Tasks[TaskIndex].parents.end());
        }
        //delete the tasks which have been allocated the resource with lowest load
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            ST.erase(TSK[RscWithMinLd][i]);
        }
        //delete the tasks which can not be performed by the resource with lowest load
        for (auto iter = ST.begin(); iter != ST.end();) {
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) ==
                Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if (ST.empty()) {
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s : ST) {
        t.push_back(pair<int, double>(s, Id[ch.RscAlcLst[s]]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    ch.RscAlcLst[t[0].first] = RscWithMinLd;
    DcdEvl_S(ch);
    CalculateEnergy1(ch);
    IFBD_S(ch);
    if (OldCh.EnergyConsumption + PrecisionValue < ch.EnergyConsumption) {
        ch = OldCh;
    }
}

void LBCA_IFBS(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Ld(comConst.NumOfRsc, 0.0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = floor(ch.Code_RK[i]);
        double TransferDataSize = Tasks[i].ExternalInputFileSizeSum;
        for (int Parent : Tasks[i].parents) {
            int PrnRscId = floor(ch.Code_RK[Parent]);
            if (RscIndex != PrnRscId) {
                TransferDataSize += ParChildTranFileSizeSum[Parent][i];
            }
        }
        Ld[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc + (TransferDataSize + Tasks[i].OFileSizeSum) / VALUE * 8 / Rscs[RscIndex].bw;;
        TSK[RscIndex].push_back(i);
    }
    int RscWithMinLd = 0;
    double TemLd = Ld[0];
    for (int j = 1; j < comConst.NumOfRsc; ++j) { //find out the resource (Rsc) with the lowest load -new-xy;
        if (TemLd > Ld[j]) {
            TemLd = Ld[j]; RscWithMinLd = j;
        }
    }
    set<int> ST;
    if (fabs(Ld[RscWithMinLd]) < PrecisionValue) {
        ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
    }
    else {
        //traverse the tasks allocated to the resource with lowest load and add their parents and children to set ST
        for (int TaskIndex : TSK[RscWithMinLd]) {
            ST.insert(Tasks[TaskIndex].children.begin(), Tasks[TaskIndex].children.end());
            ST.insert(Tasks[TaskIndex].parents.begin(), Tasks[TaskIndex].parents.end());
        }
        //delete the tasks which have been allocated the resource with lowest load
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            ST.erase(TSK[RscWithMinLd][i]);
        }
        //delete the tasks which can not be performed by the resource with lowest load
        for (auto iter = ST.begin(); iter != ST.end();) {
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) == Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if (ST.empty()) {//-w
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s : ST) {
        t.push_back(pair<int, double>(s, Ld[floor(ch.Code_RK[s])]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    double decimal = ch.Code_RK[t[0].first] - floor(ch.Code_RK[t[0].first]);
    ch.Code_RK[t[0].first] = RscWithMinLd + decimal;
    NrmDcd_S(ch);
    CalculateEnergy1(ch);
    IFBS_S(ch);
    if (OldCh.EnergyConsumption + PrecisionValue < ch.EnergyConsumption) {
        ch = OldCh;
    }
}

double CalculateEnergy2(chromosome& Chrom) {
    Chrom.EnergyConsumption = 0;
    set<double> AllTime;
    AllTime.insert(Chrom.EndTime.begin(), Chrom.EndTime.end());
    AllTime.insert(Chrom.StartTime.begin(), Chrom.StartTime.end());
    set<double>::iterator it = AllTime.begin(), endIt = AllTime.end();  --endIt;
    for (; it != endIt; ++it) {
        set<double>::iterator nxtIt = it; ++nxtIt;
        if (fabs(*nxtIt - *it) < PrecisionValue)  continue; //avoiding duplicate calculation load at the same time point
        vector<double> HstLd(HstSet.size(), 0);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            int taskId = Chrom.TskSchLst[i];
            int RscId = Chrom.RscAlcLst[taskId];
            if (Chrom.StartTime[taskId] - PrecisionValue < *it && *nxtIt - PrecisionValue < Chrom.EndTime[taskId]) {
                int taskHT = Rscs[RscId].Hostid;
                HstLd[taskHT] += Rscs[RscId].pc / HstSet[taskHT].pc;
            }
        }
        for (int k = 0; k < HstSet.size(); ++k) {
            if (Chrom.HtUseTime[k].first - PrecisionValue < *it && *nxtIt - PrecisionValue < Chrom.HtUseTime[k].second)
                Chrom.EnergyConsumption += (*nxtIt - *it) * CalculatePowerByLoad(HstLd[k], k);
        }
    }
    return Chrom.EnergyConsumption;
}

//计算EC
double CalculateEnergy1(chromosome& Chrom) {
    Chrom.EnergyConsumption = 0;
    vector<vector<int>> tskStOfHst(HstSet.size());
    vector<set<double>> tmStOfHst(HstSet.size());
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        tskStOfHst[Rscs[Chrom.RscAlcLst[i]].Hostid].push_back(i);
        tmStOfHst[Rscs[Chrom.RscAlcLst[i]].Hostid].insert(Chrom.StartTime[i]);
        tmStOfHst[Rscs[Chrom.RscAlcLst[i]].Hostid].insert(Chrom.EndTime[i]);
    }
    for (int k = 0; k < tskStOfHst.size(); ++k) { //calculate the EC of each host
        if (tskStOfHst[k].size() > 0) {
            set<double>::iterator it = tmStOfHst[k].begin(), endIt = tmStOfHst[k].end();
            --endIt;
            for (; it != endIt; ++it) {
                set<double>::iterator nxtIt = it; ++nxtIt;
                if (fabs(*nxtIt - *it) < PrecisionValue) continue;  //avoiding duplicate calculation load at the same time point
                double ld = 0;
                for (int i = 0; i < tskStOfHst[k].size(); ++i) {
                    int tsk = tskStOfHst[k][i];
                    if (Chrom.StartTime[tsk] - PrecisionValue < *it && *nxtIt - PrecisionValue < Chrom.EndTime[tsk]) {
                        ld += Rscs[Chrom.RscAlcLst[tsk]].pc / HstSet[k].pc;
                    }
                }
                Chrom.EnergyConsumption += (*nxtIt - *it) * CalculatePowerByLoad(ld, k);
            }
        }
    }
    return Chrom.EnergyConsumption;
}

double CalculateECByDelta2(chromosome& ch, int& Number) {
    double DeltaEC = 0;
    vector<int>STn; 
    int CurTask = ch.TskSchLst[Number];
    int CurHT = Rscs[ch.RscAlcLst[CurTask]].Hostid;
    for (int i = 0; i < Number; ++i) {
        int taskId = ch.TskSchLst[i];
        int RscId = ch.RscAlcLst[taskId];
        int taskHT = Rscs[RscId].Hostid;
        if ((ch.EndTime[taskId] > ch.StartTime[CurTask] + PrecisionValue && ch.EndTime[taskId] - PrecisionValue < ch.EndTime[CurTask] && taskHT == CurHT) ||
            (ch.StartTime[taskId] > ch.StartTime[CurTask] - PrecisionValue && ch.StartTime[taskId] + PrecisionValue < ch.EndTime[CurTask] && taskHT == CurHT) ||
            (ch.StartTime[taskId] + PrecisionValue < ch.StartTime[CurTask] && ch.EndTime[taskId] > ch.EndTime[CurTask] + PrecisionValue && taskHT == CurHT)) {
            STn.push_back(taskId);
        }
    }
    set<double>AllTime;
    AllTime.insert(ch.StartTime[CurTask]);
    AllTime.insert(ch.EndTime[CurTask]);
    for (int i = 0; i < STn.size(); ++i) {
        if (ch.StartTime[STn[i]] > ch.StartTime[CurTask] + PrecisionValue) {
            AllTime.insert(ch.StartTime[STn[i]]);
        }
        if (ch.EndTime[STn[i]] + PrecisionValue < ch.EndTime[CurTask]) {
            AllTime.insert(ch.EndTime[STn[i]]);
        }
    }
    set<double>::iterator it = AllTime.begin(), endIt = AllTime.end();
    --endIt;

    for (; it != endIt; ++it) {
        set<double>::iterator nxtIt = it;
        ++nxtIt;
        if (fabs(*nxtIt - *it) < 1e-6)  continue;
        double HstLd = 0;
        for (int i = 0; i < STn.size(); ++i) {
            int taskId = STn[i];
            int RscId = ch.RscAlcLst[taskId];
            if (ch.StartTime[taskId] - PrecisionValue < *it && *nxtIt - PrecisionValue < ch.EndTime[taskId]) {
                double TemLd = Rscs[RscId].pc / HstSet[CurHT].pc;
                HstLd += TemLd;
            }
        }
        double CurLd = HstLd + Rscs[ch.RscAlcLst[CurTask]].pc / HstSet[CurHT].pc;
        if (CurLd > 1)
            cout << endl;
        //cout << " CurLd :" << CurLd << endl;
        double TemEc = (*nxtIt - *it) * (CalculatePowerByLoad(CurLd, CurHT) - CalculatePowerByLoad(HstLd, CurHT));
        DeltaEC += TemEc;
    }
    if (ch.HtUseTime[CurHT].second < -PrecisionValue) {// the host is used firstly;
        DeltaEC += (ch.EndTime[CurTask] - ch.StartTime[CurTask]) * CalculatePowerByLoad(0, CurHT);
    }
    else {
        if (ch.StartTime[CurTask] < ch.HtUseTime[CurHT].first - PrecisionValue) {
            DeltaEC += (ch.HtUseTime[CurHT].first - ch.StartTime[CurTask]) * CalculatePowerByLoad(0, CurHT);
        }
        if (ch.EndTime[CurTask] > ch.HtUseTime[CurHT].second + PrecisionValue) {
            DeltaEC += (ch.EndTime[CurTask] - ch.HtUseTime[CurHT].second) * CalculatePowerByLoad(0, CurHT);
        }
    }

    return ch.EnergyConsumption + DeltaEC;
}



double CalculateECByDelta1(chromosome& ch, int& Number) {
    set<double> AllTime;
    vector<int> STn;
    int CurTask = ch.TskSchLst[Number];
    int CurHT = Rscs[ch.RscAlcLst[CurTask]].Hostid;
    for (int i = 0; i < Number; ++i) {
        int taskId = ch.TskSchLst[i];
        int RscId = ch.RscAlcLst[taskId];
        if (Rscs[RscId].Hostid == CurHT) {
            STn.push_back(taskId);
            AllTime.insert(ch.StartTime[taskId]);
            AllTime.insert(ch.EndTime[taskId]);
        }
    }
    //calculate the EC before the current task is assigned;
    double bEC = 0;
    if (!AllTime.empty()) {
        set<double>::iterator it = AllTime.begin(), endIt = AllTime.end();
        --endIt;
        for (; it != endIt; ++it) {
            set<double>::iterator nxtIt = it; ++nxtIt;
            if (fabs(*nxtIt - *it) < PrecisionValue)  continue; //avoiding duplicate calculation load at the same time point
            double HstLd = 0;
            for (int i = 0; i < STn.size(); ++i) {
                int taskId = STn[i];
                int RscId = ch.RscAlcLst[taskId];
                if (ch.StartTime[taskId] - PrecisionValue < *it && *nxtIt - PrecisionValue < ch.EndTime[taskId]) {
                    HstLd += Rscs[RscId].pc / HstSet[CurHT].pc;
                }
            }
            bEC += (*nxtIt - *it) * CalculatePowerByLoad(HstLd, CurHT);
        }
    }
    //calculate the EC after the current task is assigned;
    double aEC = 0;
    AllTime.insert(ch.StartTime[CurTask]);
    AllTime.insert(ch.EndTime[CurTask]);
    STn.push_back(CurTask);
    set<double>::iterator it = AllTime.begin(), endIt = AllTime.end();
    --endIt;
    for (; it != endIt; ++it) {
        set<double>::iterator nxtIt = it; ++nxtIt;
        if (fabs(*nxtIt - *it) < PrecisionValue)  continue; //avoiding duplicate calculation load at the same time point
        double HstLd = 0;
        for (int i = 0; i < STn.size(); ++i) {
            int taskId = STn[i];
            int RscId = ch.RscAlcLst[taskId];
            if (ch.StartTime[taskId] - PrecisionValue < *it && *nxtIt - PrecisionValue < ch.EndTime[taskId]) {
                HstLd += Rscs[RscId].pc / HstSet[CurHT].pc;
            }
        }
        aEC += (*nxtIt - *it) * CalculatePowerByLoad(HstLd, CurHT);
    }
    return ch.EnergyConsumption + aEC - bEC;
}


void UpdatePMR(vector<vector<double>>& PMR, chromosome& bstChrom) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            int count = 0;
            if (bstChrom.RscAlcLst[i] == j) {
                count = 1;
            }
            PMR[i][j] = (1 - Parameter_TSEDA.theta1) * PMR[i][j] + Parameter_TSEDA.theta1 * count;
        }
    }
}

void UpdatePMS(vector<vector<double>>& PMS, chromosome& bstChrom) {

    if (bstChrom.IsFrw) {
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            for (int j = 0; j < comConst.NumOfTsk; ++j) {
                int count = 0;
                if (bstChrom.TskSchLst[i] == j) {
                    count = 1;
                }
                PMS[j][i] = (1 - Parameter_TSEDA.theta2) * PMS[j][i] + Parameter_TSEDA.theta2 * count;
            }
        }
    }
    else {
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            for (int j = 0; j < comConst.NumOfTsk; ++j) {
                int count = 0;
                if (bstChrom.TskSchLst[comConst.NumOfTsk - i - 1] == j) {
                    count = 1;
                }
                PMS[j][i] = (1 - Parameter_TSEDA.theta2) * PMS[j][i] + Parameter_TSEDA.theta2 * count;
            }
        }
    }
}

