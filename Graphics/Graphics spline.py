import matplotlib.pyplot as plt
import numpy as np

x, y, xReal, yReal = [], [], [], []

with open("resultSpline.txt") as file:
    for line in file:
        xT, yT = line.split()
        x.append(float(xT))
        y.append(float(yT))

with open("mu") as file:
    for line in file:
        yT, xT = line.split()
        xReal.append(float(xT))
        yReal.append(float(yT))


# plot
fig, ax = plt.subplots()

plt.plot(x, y, 'y')
plt.plot(xReal, yReal, '--w')
plt.title('Mu by |B|')
plt.xlabel('|B|')
plt.ylabel('Mu')
ax.set_facecolor('blue')
plt.show()
