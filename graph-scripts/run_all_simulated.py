import subprocess

subprocess.run("python3 graph-scripts/graph_t_test.py \
    && python3 graph-scripts/graph_z_test.py \
    && python3 graph-scripts/graph_f_test.py \
    && python3 graph-scripts/graph_chi_squared.py", shell=True)