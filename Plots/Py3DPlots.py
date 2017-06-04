import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie


d = load_movie()

print("Plotting...")

plt.plot(d['bx'],d['by'])
plt.show()
