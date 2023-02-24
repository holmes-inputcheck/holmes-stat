#!/usr/bin/env python
from common import *
import scipy.stats

subprocess.run(["cmake", "."])
subprocess.run(["make"])

print("running the bench_mean_check")

fixed_point_shift = 10000
print("running the case for entries with a fixed-point shift of " + str(fixed_point_shift))

time.sleep(5)

def copy_p_value_to_csv(num_records):
    while not os.path.isfile("benchmark_result.txt"):
        time.sleep(0.01)
    f = open("benchmark_result.txt", "r")
    res = f.readlines()
    f.close()
    f = open("t_test_p_graph_dataset1.csv", "w")
    first = True
    f.write("EntriesChanged,DegreeOfFreedom,Tscore,Pvalue\n")
    for line in res:
        if first:
            first = False
            continue 
        sep = line.split(',')
        sep[-1] = sep[-1].strip('\n')
        # print(sep[-2])

        # compute t-test statistic
        p_value = scipy.stats.t.sf(abs(float(sep[-1])), df=float(sep[-2])) * 2
        sep.append(str(p_value) + "\n")
        f.write(",".join(sep))
    f.close()
    print(" done")
    os.remove("benchmark_result.txt")

num_records = 41118
write_configure_info(str(num_records) + " " + str(fixed_point_shift))
subprocess.run(["bin/graph_t_test_dataset1"])
copy_p_value_to_csv(num_records)