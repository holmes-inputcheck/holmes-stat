import subprocess

subprocess.run("python3 graph-scripts/graph_t_test_dataset1.py \
    && python3 graph-scripts/graph_z_test_dataset1.py \
    && python3 graph-scripts/graph_f_test_dataset1.py \
    && python3 graph-scripts/graph_chi_squared_dataset1.py", shell=True)