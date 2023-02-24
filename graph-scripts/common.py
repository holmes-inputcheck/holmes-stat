import subprocess
import os
import time

def configure_network():
    subprocess.run(["./network_setup.sh"])

def write_configure_info(msg):
    f = open("benchmark_input.txt", "w")
    f.write(msg)
    f.close()

def copy_benchmark_result_to_log(summary):
    while not os.path.isfile("benchmark_result.txt"):
        time.sleep(0.01)
    f = open("benchmark_result.txt", "r")
    res = f.readlines()
    f.close()
    f = open("all_results.txt", "a")
    f.write("\r\n" + summary + "\r\n")
    f.write("\r\n".join(res) + "\r\n")
    f.close()
    print(summary + " done")
    os.remove("benchmark_result.txt")