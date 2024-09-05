import pandas
import matplotlib.pyplot as plt

df = pandas.read_csv("optimisation/optimization_results_233.csv", usecols=["Beam_Current", "Cooling_Time", "Counting_Time", "Distance", "Irradiation_Time", "Total_Time"])
df2 = pandas.read_csv("optimisation/optimization_results_233_2.csv", usecols=["Beam_Current", "Cooling_Time", "Counting_Time", "Distance", "Irradiation_Time", "Total_Time"])

plt.figure(0)
plt.plot(df.Total_Time, color="steelblue", alpha=0.8)
plt.plot(df2.Total_Time, color="steelblue", alpha=0.5)
plt.hlines(y=250, color="grey", ls =':', label=f"250s", xmin = 0, xmax=250)

plt.legend()
plt.xlim(0,15)
plt.ylim(0,5000)
plt.xlabel('Iteration')
plt.ylabel('Total Time (s)')
plt.minorticks_on()
#plt.grid(True)
plt.tight_layout()
plt.savefig(f"optimisation_iteration_result_233.png")