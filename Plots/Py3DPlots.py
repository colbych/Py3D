import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie

d = load_movie()

plt.plot(d['bx'],d['by'])
