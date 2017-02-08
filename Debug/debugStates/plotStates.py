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
#plt.axis([-.2, 1.3, -.2, 1.3])
plt.triplot(x, y, triangles, '-', color='black')
plt.title('The program is currently stuck here.')
plt.xlabel('X')
plt.ylabel('Y')

i = 0
with open('States2.dat', 'r') as f:
    for i in range(noNodes):
        state = map(int, f.readline().strip().split())
        if(state[0]==0):
            al = plt.scatter(x[i], y[i], color='green', marker='s', s=80)
        elif(state[0]==1):
            nb = plt.scatter(x[i], y[i], color='blue', marker='o', s=80)
        else:
            far = plt.scatter(x[i], y[i], color='red', marker='d', s=50)

plt.legend((al, nb, far), ('Alive($\mathcal{A}$)', 'Narrow-Band($\mathcal{N}$)', 'Far-Away($\mathcal{F}$)'), scatterpoints=1, loc='upper right', fontsize=8)
plt.show()
