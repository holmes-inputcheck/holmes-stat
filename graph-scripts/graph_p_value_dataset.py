import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath,mathptmx,boldmath}']

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
#mpl.rcParams['xtick.labelsize'] = 20
#mpl.rcParams.update({'figure.max_open_warning': 0})
#mpl.rcParams.update({"font.size": 34,"figure.autolayout": True})
plt.rcParams["font.family"] = "Helvetica Neue"


columns = ["EntriesChanged", "Pvalue"]
jl_columns = ["EntriesChanged", "ActualPvalue", "UnnormalPvalue", "JLPvalue"]
dt = pd.read_csv("t_test_p_graph_dataset1.csv", usecols=columns)
dz = pd.read_csv("z_test_p_graph_dataset1.csv", usecols=columns)
df = pd.read_csv("f_test_p_graph_dataset1.csv", usecols=columns)
dc = pd.read_csv("chi_test_p_graph_dataset1.csv", usecols=jl_columns)

num_entries = 40000
threshold = 17500
x = np.arange(0, threshold, 0.1)
y5 = [0.05 for y in x]

half_linspace = [i for i in range(0, threshold)]
for i in range(threshold, num_entries):
    half_linspace.append(None)
half_linspace = np.array(half_linspace)



ax = plt.axes()
#ttest = ax.plot(dt.EntriesChanged, dt.Pvalue, 'k', label="T-test", linewidth=2)
ttest = ax.plot(dt.EntriesChanged, dt.Pvalue, color='blue', linewidth=2, label="T-test")
#ztest = ax.plot(dz.EntriesChanged, dz.Pvalue, 'b--', label="Z-test")
ztest = ax.plot(dz.EntriesChanged, dz.Pvalue, color='orange', linestyle=':', linewidth=2, label="Z-test")
#ftest = ax.plot(df.EntriesChanged, df.Pvalue, 'y', label="F-test")
ftest = ax.plot(df.EntriesChanged, df.Pvalue, color='yellow', linewidth=2, label="F-test")
#normaltest = ax.plot(dc.EntriesChanged, dc.ActualPvalue, 'r', label="Normalized Chi-square")
normaltest = ax.plot(dc.EntriesChanged, dc.ActualPvalue, color='pink', linewidth=2, label="Normalized Chi-square")

#unnormaltest = ax.plot(dc.EntriesChanged, dc.UnnormalPvalue, 'm-.', label="Unnormalized Chi-square", linewidth=1.5)
unnormaltest = ax.plot(dc.EntriesChanged, dc.UnnormalPvalue, color='red', linewidth=2, label="Unnormalized Chi-square")
#jltest = ax.plot(dc.EntriesChanged, dc.JLPvalue, 'b-.', label="JL w/ Legendre PRF Chi-square", linewidth=1.5)
jltest = ax.plot(dc.EntriesChanged, dc.JLPvalue, color='green', linestyle=':', linewidth=2, label="JL w/ Legendre PRF Chi-square")

plt.ylabel(r"\textbf{P-Value}",fontsize=10, weight='bold')
np.arange
plt.yticks(np.array([0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]), fontsize=10)
plt.xlabel(r"\textbf{Ratio of data inputs poisoned (\%)}", fontsize=10)

statvalue = ax.plot(x, y5, color='black', linestyle='--', label="P-Value 0.05")
plt.legend()

plt.xticks(np.arange(0, threshold + 1, 2000))
ax.xaxis.set_major_formatter(mtick.PercentFormatter(num_entries))

fig = plt.gcf()
fig.set_size_inches(6,4, forward=True)

plt.show()
