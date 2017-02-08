# Reference: http://matplotlib.org/examples/pylab_examples/triplot_demo.html

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


v = []
y = []
triangles = []


with open('waveSpeedCheck2.dat', 'r') as f:
    noPoints = int(f.readline())
    for i in range(noPoints):
        y_v = map(float, f.readline().strip().split())
        y.append(y_v[0])
        v.append(y_v[1])



v = [b for (a,b) in sorted(zip(y,v))]
y.sort()

plt.plot(y, v, 'o-', label='The wave speed along the Initial Wavefront')
plt.legend()
plt.show()