import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0,10,.1)
y = np.sin(t)
z = np.cos(t)

plt.figure(figsize=(10,6))
plt.plot(t, y, 'r-', label='y') # y
plt.plot(t, z, 'b-', label='z') # z
circle2 = plt.Circle((5, .5), 0.5, color='b', fill=False)
plt.gcf().gca().add_artist(circle2)
plt.xlabel('Time')
plt.legend()
plt.axis('equal')
plt.show()