"""
"""
#include <fstream>
#include <sstream>
#include "QPHH.h"
#include "GenOperator.h"
#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"
#include <chrono>
#include <cmath>

vector<vector<double>> qTable(NUMSTATES, vector<double>(NUMSACTIONS, 1.0));
map<std::string, double> lastBestFitnessMap;  

void printVector(vector<int> vec) {
	for (int value : vec) {
		cout << value << "  ";
	}
	cout << endl;
}

int getLastParent(chromosome OriginChrom,int task)
{
	vector<int> taskPosition(OriginChrom.TskSchLst.size());
	for (int i = 0; i < OriginChrom.TskSchLst.size(); ++i) {
		taskPosition[OriginChrom.TskSchLst[i]] = i;
	}

	int lastParent;
	if (Tasks[task].parents.size() <= 0)
		lastParent = 0;
	else {
		vector<int> parentLst, parentIndex;
		for (int parent : Tasks[task].parents) {
			parentIndex.push_back(taskPosition[parent]);
		}
		sort(parentIndex.begin(), parentIndex.end());

		lastParent = parentIndex[parentIndex.size() - 1];
	}
	
	return lastParent;
}

int getFirstChild(chromosome OriginChrom, int task)
{
	vector<int> taskPosition(OriginChrom.TskSchLst.size());
	for (int i = 0; i < OriginChrom.TskSchLst.size(); ++i) {
		taskPosition[OriginChrom.TskSchLst[i]] = i;
	}

	int firstChild;
	if (Tasks[task].children.size() <= 0)
		firstChild = OriginChrom.TskSchLst.size() - 1;
	else {
		vector<int> childrenLst, childIndex;

		for (int child : Tasks[task].children) {
			childrenLst.push_back(child);
			childIndex.push_back(taskPosition[child]);
		}

		sort(childIndex.begin(), childIndex.end());

		firstChild = childIndex[0];
	}

	return firstChild;
}



bool isValidTskLst(vector<int> TskSchLst) {
	vector<int> taskPosition(TskSchLst.size());

	for (int i = 0; i < TskSchLst.size(); ++i) {
		taskPosition[TskSchLst[i]] = i;
	}

	for (int i = 0; i < TskSchLst.size(); i++)
	{
		for (int parent : Tasks[i].parents) {
			if (taskPosition[parent] >= taskPosition[i]) {
				return false;
			}
		}

		
		for (int child : Tasks[i].children) {
			if (taskPosition[child] <= taskPosition[i]) {
				return false;
			}
		}
	}

	return true;
}

// Helper function to check if all children of a task are after it in the sequence
bool areAllChildrenAfter(Task task, vector<int> sequence, int taskIndex) {
	for (int chid : task.children) {
		auto it = find(sequence.begin(), sequence.end(), chid);
		if (it == sequence.end()) {
			// Parent not found in sequence
			return false;
		}
		int childIndex = distance(sequence.begin(), it);

		if (childIndex < taskIndex) {
			return false;
		}
	}
	return true;
}


// Helper function to check if all parents of a task are before it in the sequence
bool areAllParentsBefore(Task task, vector<int> sequence, int taskIndex, int& firstParentIndex) {
	for (int parent : task.parents) {
		auto it = find(sequence.begin(), sequence.end(), parent);
		if (it == sequence.end()) {
			// Parent not found in sequence
			return false;
		}
		int parentIndex = distance(sequence.begin(), it);

		if (parentIndex > taskIndex) {
			firstParentIndex = parentIndex;
			return false;
		}
	}
	return true;
}


chromosome alignRandChrom(vector<int> originTskSeq) {
	chromosome	TempRandChrom;
	for (int j = 0; j < originTskSeq.size(); j++) {
		int taskIndex = originTskSeq[j];
		int firstParentIndex;
		if (areAllParentsBefore(Tasks[taskIndex], originTskSeq, j, firstParentIndex)) { 
			continue;
		}
		else {
			originTskSeq.erase(originTskSeq.begin() + j);
			originTskSeq.insert(originTskSeq.begin() + firstParentIndex, taskIndex);
			j--;
		}
	}
	IntChr(TempRandChrom);
	TempRandChrom.TskSchLst = originTskSeq;
	return TempRandChrom;
}




/**********************************************************************
Function:     selectAction
Description:
Input:
Output:
Others:
**********************************************************************/
int SelectAction(int state) {
	int bestAction = 0;
	for (int a = 0; a < NUMSACTIONS; ++a) {
		if (qTable[state][a] > qTable[state][bestAction]) {
			bestAction = a;
		}
	}

	return bestAction;
}

int GetNextState(chromosome CurrentChrom, chromosome NewChrom)
{
	int NextState = -1;

	double logCEC = log(CurrentChrom.EnergyConsumption) + PrecisionValue;
	double logNEC = log(NewChrom.EnergyConsumption) + PrecisionValue;

	double ImproveLv = logCEC - logNEC;

	if (ImproveLv < 0)
		NextState = 0;
	else if (ImproveLv == 0)
		NextState = 1;
	else if (ImproveLv > 0.1 && ImproveLv < 1)//1e−1
		NextState = 2;
	else if (ImproveLv > 0.01 && ImproveLv < 0.1)//1e−2
		NextState = 3;
	else if (ImproveLv > 0.001 && ImproveLv < 0.01)//1e−3
		NextState = 4;
	else if (ImproveLv > 0.0001 && ImproveLv < 0.001)//1e−4
		NextState = 5;
	else if (ImproveLv > 0.00001 && ImproveLv < 0.0001)//1e−5
		NextState = 6;
	else if (ImproveLv > 0.000001 && ImproveLv < 0.00001)//1e−6
		NextState = 7;
	else if (ImproveLv > 0.0000001 && ImproveLv < 0.000001) //1e−7
		NextState = 8;
	else if (ImproveLv > 0.00000001 && ImproveLv < 0.0000001) //1e−8
		NextState = 9;
	else if (ImproveLv > 0.000000001 && ImproveLv < 0.00000001) //1e−9
		NextState = 10;
	else if (ImproveLv > 1)
		NextState = 11;
	else
		NextState = 12;

	return NextState;
}


int GetNextState(chromosome CurrentChrom, chromosome NewChrom, double entropyValue)
{
	cout << "entropyValue: " << entropyValue << endl;
	int NextState = -1;

	double logCEC = log(CurrentChrom.EnergyConsumption) + PrecisionValue;
	double logNEC = log(NewChrom.EnergyConsumption) + PrecisionValue;

	double ImproveLv = logCEC - logNEC;

	if (ImproveLv < 0)
		NextState = 0;
	else if (ImproveLv == 0) 
		NextState = 1;
	else if (ImproveLv > 0.1 && ImproveLv < 1)//1e−1
		NextState = 2;
	else if (ImproveLv > 0.01 && ImproveLv < 0.1)//1e−2
		NextState = 3;
	else if (ImproveLv > 0.001 && ImproveLv < 0.01)//1e−3
		NextState = 4;
	else if (ImproveLv > 0.0001 && ImproveLv < 0.001)//1e−4
		NextState = 5;
	else if (ImproveLv > 0.00001 && ImproveLv < 0.0001)//1e−5
		NextState = 6;
	else if (ImproveLv > 0.000001 && ImproveLv < 0.00001)//1e−6
		NextState = 7;
	else if (ImproveLv > 0.0000001 && ImproveLv < 0.000001) //1e−7
		NextState = 8;
	else if (ImproveLv > 0.00000001 && ImproveLv < 0.0000001) //1e−8
		NextState = 9;
	else if (ImproveLv > 0.000000001 && ImproveLv < 0.00000001) //1e−9
		NextState = 10;
	else if (ImproveLv > 1)
		NextState = 11;
	else
		NextState = 12;

	return NextState;
}

int GetNextState(double originEC, double newEC, double entropyValue)
{
	//cout << "entropyValue: " << entropyValue << endl;
	int NextState = -1;

	double logCEC = log(originEC) + PrecisionValue;
	double logNEC = log(newEC) + PrecisionValue;
	double ImproveLv = logCEC - logNEC;
	
	if (ImproveLv < 0)
		NextState = 0;
	else if (ImproveLv == 0) 
		NextState = 1;
	else if (ImproveLv > 0.1 && ImproveLv < 1)//1e−1
		NextState = 2;
	else if (ImproveLv > 0.01 && ImproveLv < 0.1)//1e−2
		NextState = 3;
	else if (ImproveLv > 0.001 && ImproveLv < 0.01)//1e−3
		NextState = 4;
	else if (ImproveLv > 0.0001 && ImproveLv < 0.001)//1e−4
		NextState = 5;
	else if (ImproveLv > 0.00001 && ImproveLv < 0.0001)//1e−5
		NextState = 6;
	else if (ImproveLv > 0.000001 && ImproveLv < 0.00001)//1e−6
		NextState = 7;
	else if (ImproveLv > 0.0000001 && ImproveLv < 0.000001) //1e−7
		NextState = 8;
	else if (ImproveLv > 0.00000001 && ImproveLv < 0.0000001) //1e−8
		NextState = 9;
	else if (ImproveLv > 0.000000001 && ImproveLv < 0.00000001) //1e−9
		NextState = 10;
	else if (ImproveLv > 1)
		NextState = 11;
	else
		NextState = 12;

	return NextState;
}




bool isValidSwap(const vector<int>& TskSchLst, int idx1, int idx2) {



	vector<int> taskPosition(TskSchLst.size());
	for (int i = 0; i < TskSchLst.size(); ++i) {
		taskPosition[TskSchLst[i]] = i;
	}


	for (int parent : Tasks[idx1].parents) {
		if (taskPosition[parent] >= idx2) {
			return false;
		}
	}

	for (int child : Tasks[idx1].children) {
		if (taskPosition[child] <= idx2) {
			return false;
		}
	}

	for (int parent : Tasks[idx2].parents) {
		if (taskPosition[parent] >= idx1) {
			return false;
		}
	}

	for (int child : Tasks[idx2].children) {
		if (taskPosition[child] <= idx1) {
			return false;
		}
	}
	
	return true;
}

void swapTwoTasksInSchedule0(chromosome& Chrom, set<pair<int, int>>& swappedPairs) {

	unsigned seed = 1997; 
	mt19937 gen(seed); 
	uniform_int_distribution<> dis(0, Chrom.TskSchLst.size() - 1);

	int maxAttempts = Chrom.TskSchLst.size() * (Chrom.TskSchLst.size() - 1) / 2; 
	int attempt = 0;
	chromosome originChrom = Chrom;
	while (attempt < maxAttempts) {
		int idx1 = dis(gen);
		int idx2 = dis(gen);

		if (idx1 == idx2) continue;

		if (idx1 > idx2) {
			swap(idx1, idx2);
		}

		pair<int, int> swapPair = make_pair(idx1, idx2);


		if (swappedPairs.find(swapPair) != swappedPairs.end()) {
			attempt++;
			continue;
		}

		swap(Chrom.TskSchLst[idx1], Chrom.TskSchLst[idx2]);

		alignRandChrom(Chrom.TskSchLst);
		int firstParentIndex;
		if (areAllChildrenAfter(Tasks[Chrom.TskSchLst[idx1]], Chrom.TskSchLst, idx1) &&
			areAllParentsBefore(Tasks[Chrom.TskSchLst[idx1]], Chrom.TskSchLst, idx1, firstParentIndex) &&
			areAllChildrenAfter(Tasks[Chrom.TskSchLst[idx2]], Chrom.TskSchLst, idx2) &&
			areAllParentsBefore(Tasks[Chrom.TskSchLst[idx2]], Chrom.TskSchLst, idx2, firstParentIndex)) {
			swappedPairs.insert(swapPair);
			break;
		}
		else {
			Chrom = originChrom;
			attempt++;
		}

	}
}
void swapTwoTasksInSchedule(chromosome& Chrom, set<pair<int, int>>& swappedPairs) {


	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(0, Chrom.TskSchLst.size() - 1);

	int maxAttempts = Chrom.TskSchLst.size() * (Chrom.TskSchLst.size() - 1) / 2; 
	int attempt = 0;
	chromosome originChrom = Chrom;
	while (attempt < maxAttempts) {
		int idx1 = dis(gen);
		int idx2 = dis(gen);

		if (idx1 == idx2) continue;


		if (idx1 > idx2) {
			swap(idx1, idx2);
		}

		pair<int, int> swapPair = make_pair(idx1, idx2);


		if (swappedPairs.find(swapPair) != swappedPairs.end()) {
			attempt++;
			continue;
		}

		swap(Chrom.TskSchLst[idx1], Chrom.TskSchLst[idx2]);

		alignRandChrom(Chrom.TskSchLst);
		if (isValidTskLst(Chrom.TskSchLst)) {
			swappedPairs.insert(swapPair);
			break;
		}
		else {
			Chrom = originChrom;
			attempt++;
		}
	}
}

bool isValidInsertNew(vector<int> TskSchLst) {

	for (int i = 0; i < TskSchLst.size(); i++)
	{

		for (int parent : Tasks[i].parents) {
			if (TskSchLst[parent] >= TskSchLst[i]) {
				return false;
			}
		}


		for (int child : Tasks[i].children) {
			if (TskSchLst[child] <= TskSchLst[i]) {
				return false;
			}
		}
	}
	return true;
}




chromosome insertTaskInSchedule(chromosome& Chrom) {
	unsigned seed = 1997; 
	mt19937 gen(seed); 
	uniform_int_distribution<> dis(0, Chrom.TskSchLst.size() - 1);

	set<int> triedTasks; 
	vector<int> originTskSchLst = Chrom.TskSchLst;
	int totalTasks = Chrom.TskSchLst.size();

	while (triedTasks.size() < totalTasks) {
		int idx = dis(gen);

		if (triedTasks.find(Chrom.TskSchLst[idx]) != triedTasks.end()) {
			continue; 
		}

		triedTasks.insert(Chrom.TskSchLst[idx]);


		vector <int> tempTskSchLst = Chrom.TskSchLst;
		set<int> generatedPositions; 
		generatedPositions.insert(idx); 
		int attempts = 0;
		const int maxAttempts = Chrom.TskSchLst.size();

		while (attempts < maxAttempts - 1) {

			int task = Chrom.TskSchLst[idx];

			int newPosition = dis(gen);

			if (newPosition == idx || generatedPositions.find(newPosition) != generatedPositions.end()) {
				continue; 
			}
			generatedPositions.insert(newPosition); 


			Chrom.TskSchLst.erase(Chrom.TskSchLst.begin() + idx); 
			Chrom.TskSchLst.insert(Chrom.TskSchLst.begin() + newPosition, task); 

			if (isValidTskLst(Chrom.TskSchLst))
				return Chrom;
			else {
				Chrom.TskSchLst = originTskSchLst;
				attempts++;
			}
				
		}

	}

}


chromosome insertTaskInScheduleFrontHal(chromosome& Chrom) {
	unsigned seed = 1997; 
	mt19937 gen(seed); 
	uniform_int_distribution<> dis(0, Chrom.TskSchLst.size() - 1);

	set<int> triedTasks; 
	vector<int> originTskSchLst = Chrom.TskSchLst;
	int totalTasks = Chrom.TskSchLst.size();

	while (triedTasks.size() < totalTasks) {
		int idx = dis(gen); 

		if (triedTasks.find(Chrom.TskSchLst[idx]) != triedTasks.end()) {
			continue; 
		}

		triedTasks.insert(Chrom.TskSchLst[idx]); 


		vector <int> tempTskSchLst = Chrom.TskSchLst;
		set<int> generatedPositions; 
		generatedPositions.insert(idx); 
		int attempts = 0;
		const int maxAttempts = Chrom.TskSchLst.size() / 2;


		uniform_int_distribution<> frontHalfDis(0, Chrom.TskSchLst.size() / 2 - 1);

		while (attempts < maxAttempts - 1) {


			int task = Chrom.TskSchLst[idx];

			int newPosition = frontHalfDis(gen); 

			if (newPosition == idx || generatedPositions.find(newPosition) != generatedPositions.end()) {
				continue;
			}
			generatedPositions.insert(newPosition); 

	
			Chrom.TskSchLst.erase(Chrom.TskSchLst.begin() + idx); 
			Chrom.TskSchLst.insert(Chrom.TskSchLst.begin() + newPosition, task); 



			if (isValidTskLst(Chrom.TskSchLst))
				return Chrom;
			else {
				Chrom.TskSchLst = originTskSchLst;
				attempts++;
			}

		}
	}

}


chromosome insertTaskInScheduleBackHal(chromosome& Chrom) {
	unsigned seed = 1997; 
	mt19937 gen(seed); 
	uniform_int_distribution<> dis(0, Chrom.TskSchLst.size() - 1);

	set<int> triedTasks; 
	vector<int> originTskSchLst = Chrom.TskSchLst;
	int totalTasks = Chrom.TskSchLst.size();

	while (triedTasks.size() < totalTasks) {
		int idx = dis(gen);

		if (triedTasks.find(Chrom.TskSchLst[idx]) != triedTasks.end()) {
			continue; 
		}

		triedTasks.insert(Chrom.TskSchLst[idx]);

		vector <int> tempTskSchLst = Chrom.TskSchLst;
		set<int> generatedPositions; 
		generatedPositions.insert(idx); 
		int attempts = 0;
		const int maxAttempts = Chrom.TskSchLst.size() / 2;


		uniform_int_distribution<> backHalfDis(Chrom.TskSchLst.size() / 2, Chrom.TskSchLst.size() - 1);

		while (attempts < maxAttempts - 1) {

			int task = Chrom.TskSchLst[idx];

			int newPosition = backHalfDis(gen); 

			if (newPosition == idx || generatedPositions.find(newPosition) != generatedPositions.end()) {
				continue; 
			}
			generatedPositions.insert(newPosition); 

			Chrom.TskSchLst.erase(Chrom.TskSchLst.begin() + idx); 
			Chrom.TskSchLst.insert(Chrom.TskSchLst.begin() + newPosition, task);


			if (isValidTskLst(Chrom.TskSchLst))
				return Chrom;
			else {
				Chrom.TskSchLst = originTskSchLst;
				attempts++;
			}

		}
		//cout << endl;
	}

}

bool compare(int a, int b) {
	return a > b;  
}

void insertTaskToBstPosition(chromosome& Chrom) {
	unsigned seed = 1997; 
	mt19937 gen(seed);
	uniform_int_distribution<> dis(0, Chrom.TskSchLst.size() - 1);

	chromosome OriginChrom = Chrom;

	GnrML_Evl_MEC_S_QPHH(Chrom);
	//cout << "1 " << endl;
	chromosome TempChrom = Chrom;
	int idx = dis(gen); 
	int task = Chrom.TskSchLst[idx]; 

	int lastParent = getLastParent(OriginChrom, task);
	int firstChild = getFirstChild(OriginChrom, task);

	for (int i = lastParent + 1; i < firstChild; i++) 
	{
		if (idx == i) continue;
		Chrom = OriginChrom;

		
		Chrom.TskSchLst.erase(Chrom.TskSchLst.begin() + idx); 
		Chrom.TskSchLst.insert(Chrom.TskSchLst.begin() + i, task); 


		if (isValidTskLst(Chrom.TskSchLst)) {
			GnrML_Evl_MEC_S_QPHH(Chrom);
			if (Chrom.EnergyConsumption + PrecisionValue < TempChrom.EnergyConsumption)
				TempChrom = Chrom;
		}	

	}

	Chrom = TempChrom;
}


chromosome insertTaskAfterNearestParent(chromosome& Chrom) {
	unsigned seed = 1997; 
	mt19937 gen(seed); 
	uniform_int_distribution<> dis(0, Chrom.TskSchLst.size() - 1);

	chromosome OriginChrom = Chrom;

	int idx;
	do {
		idx = dis(gen); 

	} while (Tasks[Chrom.TskSchLst[idx]].parents.empty());

	

	vector<int> taskPosition(Chrom.TskSchLst.size());
	for (int i = 0; i < Chrom.TskSchLst.size(); ++i) {
		taskPosition[Chrom.TskSchLst[i]] = i;
	}

	int task = Chrom.TskSchLst[idx];
	int latestParentPosition = -1;

	for (int parent : Tasks[task].parents) {

		if (taskPosition[parent] > latestParentPosition) {
			latestParentPosition = taskPosition[parent];
			//int minParentNum = parent;
		}
	}
	

	if (latestParentPosition == idx) {
		return OriginChrom;
	}


	Chrom.TskSchLst.erase(Chrom.TskSchLst.begin() + idx); 

	if (latestParentPosition + 1 < Chrom.TskSchLst.size()) {
		Chrom.TskSchLst.insert(Chrom.TskSchLst.begin() + latestParentPosition + 1, task); 
	}
	else {
		Chrom.TskSchLst.push_back(task);
	}
	

	alignRandChrom(Chrom.TskSchLst);

	if (isValidTskLst(Chrom.TskSchLst)) {
		//cout << "True" << endl;
		return Chrom;	
	}
	else {
		cout << "False" << endl;
		return OriginChrom;	
	}

}

void swapTwoTasksInSameLvl(chromosome& Chrom, set<pair<int, int>>& swappedPairs) {

	unsigned seed = 1997; 
	mt19937 gen(seed); 
	uniform_int_distribution<> levelDist(0, TskLstInLvl.size() - 1);
	unordered_set<int> visitedLevels; 
	int layers = TskLstInLvl.size();
	int level;
	for (int attempts = 0; attempts < layers; ++attempts) {
		level = levelDist(gen);
		if (TskLstInLvl[level].size() < 2) {
			continue;
		}

	
		if (visitedLevels.find(level) != visitedLevels.end()) {
			continue;
		}

		visitedLevels.insert(level);

		vector<int> tasksInLvl = TskLstInLvl[level];

		uniform_int_distribution<> taskDist(0, tasksInLvl.size() - 1);

		int idx1 = taskDist(gen);
		int idx2 = taskDist(gen);

		if (idx1 == idx2) continue;

		
		if (idx1 > idx2) {
			swap(idx1, idx2);
		}

		pair<int, int> swapPair = make_pair(idx1, idx2);

	
		if (swappedPairs.find(swapPair) != swappedPairs.end()) {
			continue;
		}
		chromosome originChrom = Chrom;
		swap(Chrom.TskSchLst[idx1], Chrom.TskSchLst[idx2]);

		alignRandChrom(Chrom.TskSchLst);


		if (isValidTskLst(Chrom.TskSchLst)) {
			swappedPairs.insert(swapPair);
			break;
		}
		Chrom = originChrom;
	}

}

void GenPriByRand(vector<double>& randPri, int i) {
	for (int i = 0; i < comConst.NumOfTsk; ++i) {
		randPri[i] = i;
	}
	// Obtain a time-based seed
	unsigned seed = 1997 + 10 * i;

	// Shuffle the sequence
	shuffle(randPri.begin(), randPri.end(), default_random_engine(seed));

}




/**********************************************************************
Function:     getReward
Description:
Input:
Output:
Others:
**********************************************************************/
double GetReward(double originalEC, double newEC)
{
	double reward;
	double logCEC = log(originalEC) + PrecisionValue;
	double logNEC = log(newEC)+ PrecisionValue;
	double ImproveLv = logCEC - logNEC;

	if (ImproveLv < 0) 
		reward = -10;
	else if (ImproveLv >= 1e-1)
		reward = 10;
	else if (ImproveLv >= 1e-2 && ImproveLv < 1e-1) // 1e−2
		reward = 9;
	else if (ImproveLv >= 1e-3 && ImproveLv < 1e-2) // 1e−3
		reward = 8;
	else if (ImproveLv >= 1e-4 && ImproveLv < 1e-3) // 1e−4
		reward = 7;
	else if (ImproveLv >= 1e-5 && ImproveLv < 1e-4) // 1e−5
		reward = 6;
	else if (ImproveLv >= 1e-6 && ImproveLv < 1e-5) // 1e−6
		reward = 5;
	else
		reward = 0;

	return reward;
}


/**********************************************************************
Function:     updateQValue
Description:
Input:
Output:
Others:
**********************************************************************/
void UpdateQValue(int state, int action, double reward, int nextState, double alpha, double gamma)
{
	//cout << "state1: " << state << endl;
	//cout << "nextState1: " << nextState << endl;
	//cout << "action1: " << action << endl;
	//double alpha = 0.5;
	//double gamma = 0.8;

	double maxNextQValue = qTable[nextState][0];
	for (int a = 0; a < NUMSACTIONS; ++a) {
		if (qTable[nextState][a] > maxNextQValue) {
			maxNextQValue = qTable[nextState][a];
		}
	}
	double temp = alpha * (reward + gamma * maxNextQValue - qTable[state][action]);
		qTable[state][action] += temp;

}




chromosome generateRandChrom(int count) {

	vector<int> originTskSeq, BstRandTskSeq;
	chromosome BstRandChrom;
	BstRandChrom.EnergyConsumption = 999999999;
	for (int i = 0; i < count; i++)
	{
		vector<double> Rank_rand(comConst.NumOfTsk, 0);
		vector<double> RandPri(comConst.NumOfTsk, 0);
		chromosome TempRandChrom;
		GenPriByRand(RandPri, i);

		originTskSeq = GetOriginSequence(RandPri); 
		TempRandChrom = alignRandChrom(originTskSeq);
		GnrML_Evl_MEC_S_QPHH(TempRandChrom);

		if (TempRandChrom.EnergyConsumption + PrecisionValue < BstRandChrom.EnergyConsumption)
			BstRandChrom = TempRandChrom;

	}
	return BstRandChrom;
}

chromosome generateRandChrom(int count, int j) {
	
	vector<int> originTskSeq, BstRandTskSeq;
	chromosome BstRandChrom;
	BstRandChrom.EnergyConsumption = 999999999;
	for (int i = 0; i < count; i++)
	{
		vector<double> Rank_rand(comConst.NumOfTsk, 0);
		vector<double> RandPri(comConst.NumOfTsk, 0);
		chromosome TempRandChrom;
		GenPriByRand(RandPri, j);
		
		originTskSeq = GetOriginSequence(RandPri);
		TempRandChrom = alignRandChrom(originTskSeq);
		GnrML_Evl_MEC_S_QPHH(TempRandChrom);

		if (TempRandChrom.EnergyConsumption + PrecisionValue < BstRandChrom.EnergyConsumption)
			BstRandChrom = TempRandChrom;

	}
	return BstRandChrom;
}


pair<int, int> findRandomParentChildPair(const chromosome NewChrom) {
	vector<pair<int, int>> parentChildPairs;

	
	for (int taskId : NewChrom.TskSchLst) {
		const Task& task = Tasks[taskId];
		
		int parent, child;
		if (task.parents.size() == 1 || task.children.size() == 1) {
			if (task.parents.size() == 1) {
				parent = task.parents[0];
				parentChildPairs.emplace_back(parent, taskId);
			}
			else if (task.children.size() == 1) {
				child = task.children[0];
				parentChildPairs.emplace_back(taskId, child);
			}
				
		}
	}



	if (!parentChildPairs.empty()) {

	
		sort(parentChildPairs.begin(), parentChildPairs.end());


		auto last = unique(parentChildPairs.begin(), parentChildPairs.end());

		
		parentChildPairs.erase(last, parentChildPairs.end());

		int seed = 1997;
		srand(static_cast<unsigned>(seed));
		int randomIndex = rand() % parentChildPairs.size();
		return parentChildPairs[randomIndex];
	}

	
	return { -1, -1 };
}

chromosome insertTaskPairsToBstPosition(chromosome& NewChrom) {

	chromosome tempBstChrom = NewChrom;
	GnrML_Evl_MEC_S_QPHH(tempBstChrom);

	pair<int, int> ParChildPair = findRandomParentChildPair(NewChrom);
	if (ParChildPair.first == -1 || ParChildPair.second == -1)
		return NewChrom;

	int parent = ParChildPair.first;
	int child = ParChildPair.second;


	NewChrom.TskSchLst.erase(remove(NewChrom.TskSchLst.begin(), NewChrom.TskSchLst.end(), parent), NewChrom.TskSchLst.end());
	NewChrom.TskSchLst.erase(remove(NewChrom.TskSchLst.begin(), NewChrom.TskSchLst.end(), child), NewChrom.TskSchLst.end());


	chromosome afterDelTskPariChrom = NewChrom;
	int lastParent = getLastParent(tempBstChrom, parent);
	int firstChild = getFirstChild(tempBstChrom, child);


	for (int i = lastParent + 1; i < firstChild - 2; i++) {
		afterDelTskPariChrom = NewChrom;

		afterDelTskPariChrom.TskSchLst.insert(afterDelTskPariChrom.TskSchLst.begin() + i, parent);
		afterDelTskPariChrom.TskSchLst.insert(afterDelTskPariChrom.TskSchLst.begin() + i + 1, child);
		
		if (isValidTskLst(afterDelTskPariChrom.TskSchLst)) {
			GnrML_Evl_MEC_S_QPHH(afterDelTskPariChrom);
			if (afterDelTskPariChrom.EnergyConsumption + PrecisionValue < tempBstChrom.EnergyConsumption)
				tempBstChrom = afterDelTskPariChrom;
		}
	}
	NewChrom = tempBstChrom;
	return tempBstChrom;
}

double getEntropyValue(vector<chromosome> elitePopulation) {

	
	int chromSize = elitePopulation[0].TskSchLst.size();
	vector<vector<int>> taskCount(chromSize, vector<int>(chromSize, 0));
	
	for (const auto& chrom : elitePopulation) {
		for (int pos = 0; pos < chromSize; ++pos) {
			int task = chrom.TskSchLst[pos];
			taskCount[task][pos]++;
		}
	}

	
	vector<vector<double>> taskProbability(chromSize, vector<double>(chromSize, 0.0));
	int eliteSize = elitePopulation.size();
	for (int task = 0; task < chromSize; ++task) {
		for (int pos = 0; pos < chromSize; ++pos) {
			taskProbability[task][pos] = static_cast<double>(taskCount[task][pos]) / eliteSize;
		}
	}

	
	vector<double> singleChromEntropyValue(chromSize, 0.0);
	double sum = 0.0;
	for (int task = 0; task < chromSize; ++task) {
		for (int pos = 0; pos < chromSize; ++pos) {
			if (taskProbability[task][pos] <= 0)
				continue;
			singleChromEntropyValue[task] -=  taskProbability[task][pos] * log2(taskProbability[task][pos]);
		}
		sum += singleChromEntropyValue[task];
	}

	double averageEntropyValue = sum / chromSize;
	return averageEntropyValue;

}
double calculateMedian(vector<chromosome> population) {
	
	sort(population.begin(), population.end());

	int n = population.size();
	if (n % 2 == 0) {
		
		return (population[n / 2 - 1].EnergyConsumption + population[n / 2].EnergyConsumption) / 2.0;
	}
	else {
		
		return  population[n / 2].EnergyConsumption;
	}
}


// 统计频率
vector<unordered_map<int, double>> computePositionFrequencies(const vector<chromosome>& population) {
	size_t sequenceLength = population[0].TskSchLst.size();
	vector<unordered_map<int, double>> positionFrequencies(sequenceLength);

	vector<int> positionCounts(sequenceLength, 0);

	for (const auto& individual : population) {
		for (size_t pos = 0; pos < sequenceLength; ++pos) {
			int element = individual.TskSchLst[pos];
			positionFrequencies[pos][element]++;
			positionCounts[pos]++;
		}
	}

	for (size_t pos = 0; pos < sequenceLength; ++pos) {
		for (auto it = positionFrequencies[pos].begin(); it != positionFrequencies[pos].end(); ++it) {
			it->second /= positionCounts[pos];
		}
	}

	return positionFrequencies;
}
vector<chromosome> generateNewIndividuals(const vector<unordered_map<int, double>>& positionFrequencies, size_t numIndividuals) {
	size_t sequenceLength = positionFrequencies.size();
	vector<chromosome> newIndividuals(numIndividuals);

	random_device rd;
	//int seed = 1997;
	mt19937 gen(rd());

	for (auto& individual : newIndividuals) {
		IntChr(individual); 
		individual.TskSchLst.resize(sequenceLength);

		vector<int> allElements(sequenceLength);
		iota(allElements.begin(), allElements.end(), 0); 

		for (size_t pos = 0; pos < sequenceLength; ++pos) {
			vector<int> elements;
			vector<double> probabilities;

			for (auto it = positionFrequencies[pos].begin(); it != positionFrequencies[pos].end(); ++it) {
				elements.push_back(it->first);
				probabilities.push_back(it->second);
			}

			discrete_distribution<> dist(probabilities.begin(), probabilities.end());


			shuffle(allElements.begin(), allElements.end(), gen); 
			individual.TskSchLst = allElements; 

	
			vector<int> uniqueElements(allElements.begin(), allElements.begin() + sequenceLength);
			individual.TskSchLst = uniqueElements;
		}
		individual = alignRandChrom(individual.TskSchLst);
		GnrML_Evl_MEC_S_QPHH(individual);
	}

	return newIndividuals;
}


void replaceWorstIndividuals(vector<chromosome>& population, vector<chromosome> newIndividuals) {	


	int numNew = newIndividuals.size();
	int numPop = population.size();
	vector<chromosome> mergerIndividuals;


	mergerIndividuals.reserve(numNew + numPop);


	mergerIndividuals.insert(mergerIndividuals.end(), newIndividuals.begin(), newIndividuals.end());
	mergerIndividuals.insert(mergerIndividuals.end(), population.begin(), population.end());

	sort(mergerIndividuals.begin(), mergerIndividuals.end());

	copy(mergerIndividuals.begin(), mergerIndividuals.begin() + numPop, population.begin());
}


chromosome GnrTskLstByReadyTskLst() {
	int seed = 1997;
	chromosome chrom;
	IntChr(chrom);
	vector<int > upr(comConst.NumOfTsk, 0);
	list<int> RTI;


	for (int i = 0; i < comConst.NumOfTsk; ++i) {
		upr[i] = Tasks[i].parents.size();
		if (upr[i] == 0)  RTI.push_back(i);
	}
	srand(seed); 

	for (int i = 0; i < comConst.NumOfTsk; ++i) {
		
		auto it = RTI.begin();
		std::advance(it, rand() % RTI.size());
		int selectedTask = *it;

		chrom.TskSchLst[i] = selectedTask;
		RTI.erase(it);

		for (int child : Tasks[selectedTask].children) {
			upr[child]--;
			if (upr[child] == 0) {
				RTI.push_back(child);
			}
		}
	}

	return chrom;
}


chromosome GnrTskLstByReadyTskLstRand(int i) {
	int seed = 1997 + i;
	chromosome chrom;
	IntChr(chrom);
	vector<int > upr(comConst.NumOfTsk, 0);
	list<int> RTI;


	for (int i = 0; i < comConst.NumOfTsk; ++i) {
		upr[i] = Tasks[i].parents.size();
		if (upr[i] == 0)  RTI.push_back(i);
	}
	srand(seed); 
	for (int i = 0; i < comConst.NumOfTsk; ++i) {

		auto it = RTI.begin();
		std::advance(it, rand() % RTI.size());
		int selectedTask = *it;

		chrom.TskSchLst[i] = selectedTask;
		RTI.erase(it);


		for (int child : Tasks[selectedTask].children) {
			upr[child]--;
			if (upr[child] == 0) {
				RTI.push_back(child);
			}
		}
	}

	return chrom;
}

bool compareTasksPair(const std::pair<int, int>& a, const std::pair<int, int>& b) {
	return b.second < a.second;
}

bool hasDuplicates(const std::vector<int>& vec) {
	std::unordered_set<int> uniqueElements;
	for (int elem : vec) {
		if (!uniqueElements.insert(elem).second) { 
			return true;
		}
	}
	return false;
}


bool compareTaskWithChildAndParentCount(const pair<int, pair<int, int>>& a,
	const pair<int, pair<int, int>>& b) {
	if (a.second.first != b.second.first)
		return a.second.first > b.second.first; 
	return a.second.second < b.second.second; 
}

void IGFiveInsert(chromosome& Chrom, int selectCount) {

	if (hasDuplicates(Chrom.TskSchLst))
		cout << endl;


	random_device rd;
	mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, Chrom.TskSchLst.size() - 1);

	chromosome OriginChrom = Chrom;
	GnrML_Evl_MEC_S_QPHH(Chrom);
	chromosome TempChrom = Chrom; 


	std::vector<int> selectedTasks;
	while (selectedTasks.size() < selectCount) {
		int idx = dis(gen); 
		if (std::find(selectedTasks.begin(), selectedTasks.end(), Chrom.TskSchLst[idx]) == selectedTasks.end()) {
			selectedTasks.push_back(Chrom.TskSchLst[idx]); 
		}
	}


	std::vector<std::pair<int, int>> taskWithChildCount;
	for (int task : selectedTasks) {
		int childCount = Tasks[task].children.size(); 
		taskWithChildCount.push_back({ task, childCount });
	}
	sort(taskWithChildCount.begin(), taskWithChildCount.end(), compareTasksPair);

	for (const auto& taskPair : taskWithChildCount) {
		chromosome temp2;
		IntChr(temp2);
		int task = taskPair.first;
		int idx = std::find(Chrom.TskSchLst.begin(), Chrom.TskSchLst.end(), task) - Chrom.TskSchLst.begin();
		temp2.TskSchLst = Chrom.TskSchLst;
		int lastParent = getLastParent(OriginChrom, task);
		int firstChild = getFirstChild(OriginChrom, task);

		for (int i = lastParent + 1; i < firstChild; i++) {
			if (idx == i) continue;
			Chrom = temp2; 

			Chrom.TskSchLst.erase(Chrom.TskSchLst.begin() + idx);
			chromosome temp1 = Chrom;
			Chrom.TskSchLst.insert(Chrom.TskSchLst.begin() + i, task);

			if (hasDuplicates(Chrom.TskSchLst))
				cout << endl;


			if (isValidTskLst(Chrom.TskSchLst)) {
				GnrML_Evl_MEC_S_QPHH(Chrom); 

				if (Chrom.EnergyConsumption + PrecisionValue < TempChrom.EnergyConsumption) {
					TempChrom = Chrom; 
				}
			}
		}
	}

	Chrom = TempChrom;
}

void logBestFitnessChange(const string& Model, int iteration, const string& filename, const chromosome& CurrentBstChrom) {

	static auto startTime = std::chrono::high_resolution_clock::now(); 

	double& lastBestFitness = lastBestFitnessMap[Model];

	if (CurrentBstChrom.EnergyConsumption + PrecisionValue <  lastBestFitness || lastBestFitness == 0) {
		std::ofstream outFile(filename, std::ios::app);
		if (!outFile) {
			std::cerr << "cannot open: " << filename << std::endl;
			return;
		}

		auto currentTime = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsedTime = currentTime - startTime;

		outFile << Model << "\t" << comConst.NumOfTsk << "\t" <<  iteration << "\t" << elapsedTime.count() << "\t" 
			<< CurrentBstChrom.EnergyConsumption << std::endl;

		lastBestFitness = CurrentBstChrom.EnergyConsumption;
	}
}


void correctDuplicates(vector<int>& child, const vector<int>& otherChild,
	int start, int end) {
	unordered_set<int> usedInChild;
	vector<int> availableNumbers;


	for (int i = 0; i < child.size(); ++i) {
		if (i < start || i > end) {
			usedInChild.insert(child[i]);
		}
	}


	for (int i = 0; i < otherChild.size(); ++i) {
		if ((i < start || i > end) &&
			usedInChild.find(otherChild[i]) == usedInChild.end()) {
			availableNumbers.push_back(otherChild[i]);
			usedInChild.insert(otherChild[i]); 
		}
	}


	int availIndex = 0;
	for (int i = start; i <= end; ++i) {
		if (usedInChild.find(child[i]) != usedInChild.end()) {
			if (availIndex < availableNumbers.size()) {
				child[i] = availableNumbers[availIndex++];
			}
			else {
				
				cerr << "Error: Not enough unique elements for correction" << endl;
			}
		}
		usedInChild.insert(child[i]);
	}
}


void Cross_DP_ModifyOnlyFirst(chromosome& ch1, const chromosome& ch2) {
	
	if (ch1.TskSchLst.size() != ch2.TskSchLst.size() || ch1.TskSchLst.size() < 2) {
		cerr << "Parents must have the same size and at least 2 elements" << endl;
		return;
	}

	
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dist(0, ch1.TskSchLst.size() - 1);

	int point1 = dist(gen);
	int point2 = dist(gen);

	if (point1 > point2) {
		swap(point1, point2);
	}


	vector<int> temp = ch1.TskSchLst;

	for (int i = point1; i <= point2; ++i) {
		temp[i] = ch2.TskSchLst[i]; 
	}

	correctDuplicates(temp, ch2.TskSchLst, point1, point2);


	alignRandChrom(temp);

	ch1.TskSchLst = temp;
}



void Crossover_QPHH(chromosome& ch1, const chromosome ch2) {


	chromosome y1 = ch1;
	chromosome y2 = ch2;

	//Cross_SP(x1, x2);
	Cross_DP_ModifyOnlyFirst(y1, y2);
}


void localSearch_Cross(vector <chromosome> sortTempElitePopulation, const chromosome& Chrom, int savePopSize, vector <chromosome>& BestPopulationAfterLS) {

	//cout << "before: ";
	for (int i = 0; i < savePopSize; i++)
	{
		IntChr(BestPopulationAfterLS[i]);
		BestPopulationAfterLS[i].TskSchLst = sortTempElitePopulation[savePopSize - 1 -i].TskSchLst;


	}
	//cout << "LS_Cross" << endl;

	vector<future<void>> futures;
	for (int i = 0; i < BestPopulationAfterLS.size(); i++) {
		futures.push_back(async(launch::async, [&, i] {

			Crossover_QPHH(BestPopulationAfterLS[i], Chrom);
			GnrML_Evl_MEC_S_QPHH(BestPopulationAfterLS[i]);
			}));
	}

	for (auto& f : futures) {
		f.get();
	}
}
void localSearch(vector <chromosome> sortTempElitePopulation, int selectCount, int savePopSize, vector <chromosome>& BestPopulationAfterLS) {

	for (int i = 0; i < savePopSize; i++)
	{
		IntChr(BestPopulationAfterLS[i]);
		//BestPopulationAfterLS[i].TskSchLst = sortTempElitePopulation[savePopSize - 1 -i].TskSchLst;
		BestPopulationAfterLS[i].TskSchLst = sortTempElitePopulation[i].TskSchLst;

	}
	//cout << "LS" << endl;

	vector<future<void>> futures;
	for (int i = 0; i < BestPopulationAfterLS.size(); i++) {
		futures.push_back(async(launch::async, [&, i] {
			IGFiveInsert(BestPopulationAfterLS[i], selectCount);
			}));
	}

	for (auto& f : futures) {
		f.get();
	}
}

chromosome runQPHH(string Model,string ECfileName,string XmlFile, string RscAlcFile, double& SchTime, int& iteration, double PopSize, int selectCount, double epslion, double gamma) {

	for (int i = 0; i < NUMSTATES; ++i) {
		for (int j = 0; j < NUMSACTIONS; ++j) {
			qTable[i][j] = 1.0;
		}
	}
	clock_t start = clock();
	ReadFile(XmlFile, RscAlcFile);
	ConfigParameter_QPHH();
	CalculateLevelList();


	vector<double> Rank_b(comConst.NumOfTsk, 0);
	vector<double> ww(comConst.NumOfTsk, 0);
	W_Cal_Average_S(ww);

	Calculate_Rank_b_S(Rank_b, ww);


	chromosome BstChrom = GnrChr_HMEC_S(Rank_b);
	chromosome Chrom_HEFT = GnrChr_HEFT_S(Rank_b);

	double originECSum = 0.0;
	vector<chromosome> population;
	int populationSize = 500;
	for (int i = 0; i < populationSize; i++)
	{
		chromosome RandChrom = GnrTskLstByReadyTskLstRand(i);
		GnrML_Evl_MEC_S_QPHH(RandChrom);
		population.push_back(RandChrom);
	}

	 double originECMedian = calculateMedian(population);

	sort(population.begin(), population.end());
	int eliteNum = PopSize;
	vector<chromosome> elitePopulation(population.begin(), population.begin() + eliteNum);

	chromosome BstRandChrom = elitePopulation[0];

	if (Chrom_HEFT.EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption) {
		BstChrom = Chrom_HEFT;
	}

	if (BstRandChrom.EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption) {
		BstChrom = BstRandChrom;
	}


	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distribAction(0, NUMSACTIONS - 1);
	uniform_int_distribution<> distribState(0, NUMSTATES - 1);
	uniform_int_distribution<> distribPop(0, population.size() - 1);



	int CurrentState = distribState(gen);
	double RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;

	chromosome CurrentBstChrom;

	CurrentBstChrom = BstChrom;


	vector<set<pair<int, int>>> swappedPairs(elitePopulation.size());


	vector<chromosome> tempElitePopulation(elitePopulation.size());
	for(int i = 0; i < elitePopulation.size(); i++)
	{
		IntChr(tempElitePopulation[i]);
		tempElitePopulation[i].TskSchLst = elitePopulation[i].TskSchLst;
	}

	int notImproveCount = 0;

	int action = -1;
	double LSRate = 1;

	while (RunTime + PrecisionValue < SchTime) {

		double originBstEC = BstChrom.EnergyConsumption;


		double SelectRand = ((double)rand() / (RAND_MAX)) / 1.0;
		double newECSum = 0.0;
		if (SelectRand < epslion)
			action = distribAction(gen);
		else
			action = SelectAction(CurrentState);

		int index = 0;

		switch (action) {
		case 0: {
			//cout << "action = " << action << endl;
			vector<future<void>> futures;
			for (int i = index; i < tempElitePopulation.size(); i++) {
				futures.push_back(async(launch::async, [&, i] {
					swapTwoTasksInSchedule(tempElitePopulation[i], swappedPairs[i]);
					swapTwoTasksInSchedule(tempElitePopulation[i], swappedPairs[i]);
					swapTwoTasksInSchedule(tempElitePopulation[i], swappedPairs[i]);
					GnrML_Evl_MEC_S_QPHH(tempElitePopulation[i]);
					}));
			}

			for (auto& f : futures) {
				f.get();
			}
			break;
		}
		case 1: {
			vector<future<void>> futures;

			for (int i = index; i < tempElitePopulation.size(); i++) {
				futures.push_back(async(launch::async, [&, i] {
					insertTaskInSchedule(tempElitePopulation[i]);
					insertTaskInScheduleFrontHal(tempElitePopulation[i]);
					insertTaskInScheduleBackHal(tempElitePopulation[i]);
					GnrML_Evl_MEC_S_QPHH(tempElitePopulation[i]);
					}));
			}

			for (auto& f : futures) {
				f.get();
			}
			break;
		}

		case 2: {
			vector<future<void>> futures;

			for (int i = index; i < tempElitePopulation.size(); i++) {
				futures.push_back(async(launch::async, [&, i] {
					insertTaskToBstPosition(tempElitePopulation[i]);
					}));
			}

			for (auto& f : futures) {
				f.get();
			}
			break;
		}
		case 3: {
			vector<future<void>> futures;

			for (int i = index; i < tempElitePopulation.size(); i++) {
				futures.push_back(async(launch::async, [&, i] {
					insertTaskAfterNearestParent(tempElitePopulation[i]);
					GnrML_Evl_MEC_S_QPHH(tempElitePopulation[i]);
					}));
			}

			for (auto& f : futures) {
				f.get();
			}
			break;
		}
		case 4: {
			vector<future<void>> futures;

			for (int i = index; i < tempElitePopulation.size(); i++) {
				futures.push_back(async(launch::async, [&, i] {
					swapTwoTasksInSameLvl(tempElitePopulation[i], swappedPairs[i]);
					GnrML_Evl_MEC_S_QPHH(tempElitePopulation[i]);
					}));
			}

			for (auto& f : futures) {
				f.get();
			}
			break;
		}

		case 5: {
			vector<future<void>> futures;

			for (int i = index; i < tempElitePopulation.size(); i++) {
				futures.push_back(async(launch::async, [&, i] {
					insertTaskPairsToBstPosition(tempElitePopulation[i]);
					}));
			}

			for (auto& f : futures) {
				f.get();
			}
			break;
		}

		}

		vector <chromosome> sortTempElitePopulation(tempElitePopulation.begin(), tempElitePopulation.end());
		sort(sortTempElitePopulation.begin(), sortTempElitePopulation.end());
		chromosome NewBstChrom = sortTempElitePopulation[0];

		int NextState = GetNextState(CurrentBstChrom, NewBstChrom);

		double reward = GetReward(CurrentBstChrom.EnergyConsumption, NewBstChrom.EnergyConsumption);

		double curRunTimeAlpha = (double)(clock() - start) / CLOCKS_PER_SEC;// 计算运行时间
		double alpha = 1 - (0.9 * curRunTimeAlpha) / SchTime;

		UpdateQValue(CurrentState, action, reward, NextState,alpha, gamma);

		CurrentState = NextState;
		
		CurrentBstChrom = NewBstChrom;


		if (CurrentBstChrom.EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption) {
			BstChrom = CurrentBstChrom;
		}

		int savePopSize = eliteNum / 2;
		vector <chromosome> BestPopulationAfterLS(savePopSize);


		double LSRand = ((double)rand() / (RAND_MAX)) / 1.0;

		if (LSRand < LSRate) {
			localSearch_Cross(sortTempElitePopulation, BstChrom, savePopSize, BestPopulationAfterLS);
		}
		else {
			localSearch(sortTempElitePopulation, selectCount, savePopSize, BestPopulationAfterLS);
		}
		LSRate *= 0.95;

		sort(BestPopulationAfterLS.begin(), BestPopulationAfterLS.end());
		chromosome BstChromAfterLS = BestPopulationAfterLS[0];


		vector<chromosome> TemElitePopulationAfterLS;
		TemElitePopulationAfterLS.reserve(sortTempElitePopulation.size() + BestPopulationAfterLS.size());
		TemElitePopulationAfterLS.insert(TemElitePopulationAfterLS.end(), sortTempElitePopulation.begin(), sortTempElitePopulation.end());
		TemElitePopulationAfterLS.insert(TemElitePopulationAfterLS.end(), BestPopulationAfterLS.begin(), BestPopulationAfterLS.end());
		sort(TemElitePopulationAfterLS.begin(), TemElitePopulationAfterLS.end());

		if (TemElitePopulationAfterLS[0].EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption) {
			BstChrom = TemElitePopulationAfterLS[0];
		}

		vector<chromosome> sortTempElitePopulationTemp(TemElitePopulationAfterLS.begin(), TemElitePopulationAfterLS.begin() + sortTempElitePopulation.size());


		vector<chromosome> TemElitePopulation(sortTempElitePopulationTemp.size());
		for (int i = 0; i < sortTempElitePopulationTemp.size(); i++)
		{
			IntChr(TemElitePopulation[i]);
			TemElitePopulation[i].TskSchLst = sortTempElitePopulationTemp[i].TskSchLst;
		}

		tempElitePopulation = TemElitePopulation;

		++iteration;

		RunTime = (double)(clock() - start) / CLOCKS_PER_SEC;
	}
	SchTime = (double)(clock() - start) / CLOCKS_PER_SEC;
	return BstChrom;
}
