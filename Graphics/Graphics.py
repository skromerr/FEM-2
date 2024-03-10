import matplotlib.pyplot as plt
import numpy as np

x, y, Az = [], [], []

with open("results.txt") as file:
    for line in file:
        xT, yT, AzT = line.split()
        x.append(float(xT))
        y.append(float(yT))
        Az.append(float(AzT))

x = list(dict.fromkeys(x))
y = list(dict.fromkeys(y))

levels = np.linspace(min(Az), max(Az), 8)



# plot
fig, ax = plt.subplots()

X, Y = np.meshgrid(x, y)
Az = np.reshape(Az, (len(y), len(x)))
colorBar = plt.contourf(X, Y, Az, levels=100, cmap='jet')
plt.colorbar(colorBar, ax=ax, format='%.0e')
plt.show()
