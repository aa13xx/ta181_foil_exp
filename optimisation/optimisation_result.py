import pandas
import matplotlib.pyplot as plt

df = pandas.read_csv("optimisation/optimization_results.csv", usecols=["Beam_Current", "Counting_Time", "Collection_Time", "Distance", "Irradiation_Time", "Total_Time"])

plt.figure(0)
plt.plot(df.Total_Time, color="steelblue", alpha=0.8)
plt.hlines(y=366.12, color="grey", ls =':', label=f"366.12s", xmin = 0, xmax=250)

plt.legend()
plt.xlim(0,250)
plt.ylim(0,15000)
plt.xlabel('Iteration')
plt.ylabel('Total Time (s)')
plt.minorticks_on()
#plt.grid(True)
plt.tight_layout()
plt.savefig(f"optimisation_iteration_result.png")