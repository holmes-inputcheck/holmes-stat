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
    f = open("f_test_p_graph.csv", "w")
    first = True
    f.write("EntriesChanged,Fscore,Pvalue\n")
    for line in res:
        if first:
            first = False
            continue 
        sep = line.split(',')
        sep[-1] = sep[-1].strip('\n')

        # compute p value for the F-test
        p_value = (1-scipy.stats.f.cdf(abs(float(sep[-1])), num_records - 1, num_records - 1)) * 2

        # correct p value greater than 1
        p_value = min(p_value, 1)
        
        sep.append(str(p_value) + "\n")
        f.write(",".join(sep))
    f.close()
    print(" done")
    os.remove("benchmark_result.txt")


num_records = 40000
write_configure_info(str(num_records) + " " + str(fixed_point_shift))
subprocess.run(["bin/graph_f_test"])
copy_p_value_to_csv(num_records)