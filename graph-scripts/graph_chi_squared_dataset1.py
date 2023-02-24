#!/usr/bin/env python
from common import *
import scipy.stats
import numpy as np

subprocess.run(["cmake", "."])
subprocess.run(["make"])

print("running the graph chi squared for dataset 1")

k = 200
num_records = 41118

time.sleep(5)

def find_unnormalized_p_value(num_records):
    while not os.path.isfile("benchmark_result.txt") or not os.path.isfile("chi_expectation.txt"):
        time.sleep(0.01)
    f = open("benchmark_result.txt", "r")
    res = f.readlines()
    f.close()
    
    expec = []

    fe = open("chi_expectation.txt", "r")
    fer = fe.readlines()
    fe.close()
    fer_first = True
    for line in fer:
        # preprocess csv
        if fer_first:
            fer_first = False
            continue
        sep = line.strip('\n')

        # fetch expectation
        expec.append(float(sep))

    f = open("chi_test_p_graph_dataset1.csv", "w")
    first = True
    f.write("EntriesChanged,ActualPvalue,UnnormalPvalue,JLPvalue\n")

    list_iters = []
    list_jlchi = []
    list_unchi = []
    list_actualchi = []

    for line in res:
        # preprocess csv
        if first:
            first = False
            continue 
        sep = line.split(',')
        sep[-1] = sep[-1].strip('\n')

        # parse entries
        num_iters = int(sep[0])
        actualchi = float(sep[1])
        unchi = float(sep[2])
        jlchi = float(sep[3])

        # save entries
        list_iters.append(num_iters)
        list_actualchi.append(actualchi)
        list_unchi.append(unchi)
        list_jlchi.append(jlchi)

    # Number of nonzero features in public dataseet
    nonzero_D = 643

    prob = [x / num_records for x in expec]
    cov = [[0 for j in range(len(prob))] for i in range(len(prob))]
    mean = [0 for i in range(len(prob))]

    # compute covariance
    for i in range(len(prob)):
        for j in range(len(prob)):
            cov[i][j] = prob[i] * prob[j]
  
    # compute variance diagonals
    for i in range(len(prob)):
        cov[i][i] = prob[i] * (1 - prob[i])

    cov = np.array(cov)

    rng = np.random.default_rng()
    num_samples = 10000
    multi = rng.multivariate_normal(mean, cov, size=num_samples)

    stat_values = []
    for sample in multi:
        stat_values.append(np.sum(np.square(sample))) # multiply the array by num_records in original
    stat_values.sort()

    print("---")
    print("0% quartile: ", stat_values[0])
    print("25% quartile: ", stat_values[num_samples // 4])
    print("50% quartile: ", stat_values[2 * num_samples // 4])
    print("75% quartile: ", stat_values[3 * num_samples // 4])
    print("100% quartile: ", stat_values[-1])

    for i in range(len(list_actualchi)):
        actual_p_value = scipy.stats.chi2.sf(list_actualchi[i], df=nonzero_D - 1)

        closest = min(stat_values, key=lambda x:abs(x-list_unchi[i]))
        desired_index = stat_values.index(closest)
        unnormal_p_value = 1 - desired_index / num_samples

        jl_closest = min(stat_values, key=lambda x:abs(x-list_jlchi[i]))
        jl_desired_index = stat_values.index(jl_closest)
        jl_p_value = 1 - jl_desired_index / num_samples
        
        f.write(str(list_iters[i]) + "," + str(actual_p_value) + "," + str(unnormal_p_value) + "," + str(jl_p_value) + "\n") 

    f.close()
    print(" done")
    os.remove("benchmark_result.txt")

write_configure_info(str(k) + " " + str(num_records))
subprocess.run(["bin/graph_chi_squared_dataset1"])
find_unnormalized_p_value(num_records)