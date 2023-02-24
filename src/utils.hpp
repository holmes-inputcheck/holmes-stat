#ifndef EMP_ZK_UTILS_H
#define EMP_ZK_UTILS_H

#include "emp-zk-holmes/emp-zk.h"
#include <cstdio>
#include <iostream>
#include "emp-tool/emp-tool.h"
#if defined(__linux__)
#include <sys/time.h>
	#include <sys/resource.h>
#elif defined(__APPLE__)
#include <unistd.h>
#include <sys/resource.h>
#include <mach/mach.h>
#endif

using namespace emp;
using namespace std;

int get_num_range_bits(uint64_t B_low, uint64_t B_high) {
	return ceil(log2(B_high - B_low + 1));
}

#endif //EMP_ZK_UTILS_H
