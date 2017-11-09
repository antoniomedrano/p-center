import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0,10,.1)
y = np.sin(t)
z = np.cos(t)

plt.figure(figsize=(10,6))
plt.plot(t, y, 'r-', label='y') # y
plt.plot(t, z, 'b-', label='z') # z
plt.plot(5, .5, 'go')
circle1 = plt.Circle((6, .6), 0.5, color='k', fill=False)
plt.gcf().gca().add_artist(circle1)
circle2 = plt.Circle((5, .5), 0.5, color='b', fill=False)
plt.gcf().gca().add_artist(circle2)
plt.plot([.3, .4], [.3, .4], color='green', marker='o',
     markerfacecolor='blue', markersize=12)
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()
plt.axis('equal')
plt.show()