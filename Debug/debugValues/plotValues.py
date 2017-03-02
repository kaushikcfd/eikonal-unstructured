# Reference: http://matplotlib.org/examples/pylab_examples/triplot_demo.html

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})

x = []
y = []
T = []
triangles = []

#-----------------------------------------------------------------------------
# Reading values from the files
#-----------------------------------------------------------------------------

with open('square-refined-Nodes.dat', 'r') as f:
    noNodes = int(f.readline())
    for i in range(noNodes):
        x_y = map(float, f.readline().strip().split())
        x.append(x_y[0])
        y.append(x_y[1])

with open('square-refined-Elements.dat', 'r') as f:
    noElements = int(f.readline())
    for i in range(noElements):
        node1_2_3 = map(int, f.readline().strip().split())
        triangles.append(node1_2_3)

with open('square-refined.dat', 'r') as f:
    noValues = int(f.readline())
    for i in range(noValues):
        t = float(f.readline())
        T.append(t)


#-----------------------------------------------------------------------------
# Defining the triangulation, with the help of the values that have been read.
#-----------------------------------------------------------------------------

triang = tri.Triangulation(x, y, triangles)


#-----------------------------------------------------------------------------
# Plotting the figure
#-----------------------------------------------------------------------------
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontour(triang, T, 10)
plt.colorbar()
#plt.title("Planar Wavefront on Unstructured Grid")
plt.title(r"$F = 2 - \frac{1}{2}\cos^2\left(\pi\left(y-\frac{1}{2}\right)\right)$- Rectangular Structured Grid")

plt.show()
