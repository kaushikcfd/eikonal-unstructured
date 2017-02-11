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

with open('square-structured-Nodes.dat', 'r') as f:
    noNodes = int(f.readline())
    for i in range(noNodes):
        x_y = map(float, f.readline().strip().split())
        x.append(x_y[0])
        y.append(x_y[1])

with open('square-structured-Elements.dat', 'r') as f:
    noElements = int(f.readline())
    for i in range(noElements):
        node1_2_3 = map(int, f.readline().strip().split())
        triangles.append(node1_2_3)

with open('dahiya_paper-structured.dat', 'r') as f:
    noValues = int(f.readline())
    for i in range(noValues):
        t = float(f.readline())
        T.append(t)


#-----------------------------------------------------------------------------
# Defining the triangulation, with the help of the values that have been read.
#-----------------------------------------------------------------------------

triang = tri.Triangulation(x, y, triangles)

#-----------------------------------------------------------------------------
# Refine data
#-----------------------------------------------------------------------------
refiner = tri.UniformTriRefiner(triang)
tri_refi, z_test_refi = refiner.refine_field(T, subdiv=3)

#-----------------------------------------------------------------------------
# Plotting the figure
#-----------------------------------------------------------------------------
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontour(triang, T, 10)
#plt.tricontour(tri_refi, z_test_refi, 10)
plt.colorbar()
plt.title("Planar Wavefront with no sweeping")

plt.show()
