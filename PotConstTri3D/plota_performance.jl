
using PyPlot

plt=PyPlot

x = [624,9048,9358,9772,12128,17520,17980]
y1 = [2.05165308,336.799750236,372.826819311,415.375716659,617.49832289,1086.066939104,1140.771201388]
y2 = [0.732925931,15.394120959,17.500025624,19.551630967,36.065009333,75.679239542,134.270465291]
plt.close("all")
plt.figure()
plt.plot(x,y1,linestyle="-",marker="x")
plt.axis("tight") # Fit the axis tightly to the plot
plt.title("Time to assembly matrix A and vector b")
plt.xlabel("Nuber of nodes")
plt.ylabel("Time (seconds)")
plt.grid("on")

plt.figure()
plt.plot(x,y2,linestyle="-",marker="x")
plt.axis("tight") # Fit the axis tightly to the plot
plt.title("Time to solve the linear system A x = b")
plt.xlabel("Nuber of nodes")
plt.ylabel("Time (seconds)")
plt.grid("on")


plt.show()
