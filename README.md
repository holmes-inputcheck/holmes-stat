## Setup for the graphs and prerequisites:

1. [10 minutes] Ensure that the prerequisites are built.

Install [TexLive/MacTeX](https://tug.org/texlive/)

Install python3.
```
pip3 install numpy
pip3 install matplotlib
pip3 install scipy
pip3 install pandas
pip3 install latex
```

2. [10 minutes] Either install the AWS instance that has HOLMES pre-built [here](https://github.com/holmes-inputcheck/holmes-library), or you can build two of HOLMES’ dependencies, specifically FLINT and OpenSSL [here](https://github.com/holmes-inputcheck/holmes-library#requirements).

3. [1 minute] Clone the stat graph repository and cd to holmes-stat:
```
git clone https://github.com/holmes-inputcheck/holmes-stat.git && cd holmes-stat
```

## Execution for the graphs

1. [1-2 hours] Run the graphs for the simulated dataset with:
```
python3 graph-scripts/run_all_simulated.py
```

2. [1-2 hours] Run the graphs for the real dataset with:
```
python3 graph-scripts/run_all_dataset.py
```

## Retrieving graph results

1. [1 minute] Plot the graph for the simulated dataset with:
```
​​python3 graph-scripts/graph_p_value.py
```

2. [1 minute] Plot the graph for the real dataset with:
```
python3 graph-scripts/graph_p_value_dataset.py
```
