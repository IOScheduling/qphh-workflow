#pragma once


#ifndef FRAME_COMMON_H
#define FRAME_COMMON_H


#include "classAndVarDefine.h"
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <map>
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <numeric>
#include <random>
#include <list>
#include <unordered_set>
#include <thread>
#include <unordered_map>
#include <future>


#define PrecisionValue 1e-6
#define InfiniteValue 9999999999
#define VALUE 1000000
#define XY_MIN(i, j)   (((i) > (j)) ? (j) : (i))
#define XY_MAX(i, j)   (((i) < (j)) ? (j) : (i))
#define XY_SWAP(x, y, type) {type tmp = (x); (x) = (y); (y) = (tmp);}
#endif //FRAME_COMMON_H
