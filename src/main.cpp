#include <iostream>
#include "Mesh2D.h"

#define FILE_PATH "arbitrary_unstructured.msh"

void displayArray(float *x, unsigned n)
{
    unsigned i;
    for(i=0; i<n; i++)
        printf("%d\t%6.2f\n",i+1, x[i]);
}

void displayMatrix(int **A, unsigned m, unsigned n)
{
    unsigned i,j;
    for(i=0; i<m; i++)
    {
        for(j=0;j<n;j++)
            printf("%d\t",A[i][j]);
        printf("\n");
    }
}

int main() {
    Mesh2D mesh;
    mesh.readFromFile(FILE_PATH);
    mesh.write("arbit_Nodes.dat", "arbit_Elements.dat");
    return 0;
}
