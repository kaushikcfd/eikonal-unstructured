# Reference: http://matplotlib.org/examples/pylab_examples/triplot_demo.html

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


x = []
y = []
triangles = []


with open('square_Nodes.dat', 'r') as f:
    noNodes = int(f.readline())
    for i in range(noNodes):
        x_y = map(float, f.readline().strip().split())
        x.append(x_y[0])
        y.append(x_y[1])

with open('square_Elements.dat', 'r') as f:
    noElements = int(f.readline())
    for i in range(noElements):
        node1_2_3 = map(int, f.readline().strip().split())
        triangles.append(node1_2_3)

plt.figure()
plt.gca().set_aspect('equal')
#plt.axis([-1.2, 1.2, -1.2, 1.2])
plt.triplot(x, y, triangles, '-', color='black')
plt.title('Triplot of user-specified triangulation')
plt.xlabel('X')
plt.ylabel('Y')

i = 0
with open('States.dat', 'r') as f:
    for i in range(noNodes):
        state = map(int, f.readline().strip().split())
        if(state[0]==0):
            plt.scatter(x[i], y[i], color='green', marker='s', s=80)
        elif(state[0]==1):
            plt.scatter(x[i], y[i], color='blue', marker='o', s=80)
        else:
            plt.scatter(x[i], y[i], color='red', marker='^', s=80)

plt.show()
