#include "EikonalSolver.h"


EikonalSolver::EikonalSolver(Mesh2D* _mesh, function<float(float, float)> _F, function<float(float, float)> _v1, function<float(float, float)> _v2) {
    mesh = _mesh;
    F = _F;
    v1 = _v1;
    v2 = _v2;
    
    // The functions have been intialized, now using these functions to set the values for each nodes.
   	
   	float x, y;
   	Node** nodes = mesh->nodes;
    int noNodes = mesh->getNoOfNodes();
    for( int i = 0; i < noNodes; i++) {
        x = nodes[i]-> getX();
        y = nodes[i]-> getY();
        nodes[i]->setF(F(x, y)); // Setting the wave speed for the $i^th$ node.
        nodes[i]->setv1(v1(x, y)); // Setting the medium speed (x-direction) for the $i^th$ node.
        nodes[i]->setv2(v2(x, y)); // Setting the medium speed (y-direction) for the $i^th$ node.
        if(nodes[i]->getState() == NARROW_BAND)
            narrowBandHeap.push(nodes[i]);
    }
}

void EikonalSolver::plotStates(string outputFile="States.dat") {
    FILE* pFile;
    pFile = freopen(outputFile.c_str(), "w", stdout);
    Node** nodes = mesh->nodes;
    int noNodes = mesh->getNoOfNodes();
    for( int i = 0; i < noNodes; i++) {
        if((nodes[i]->getState())==ALIVE){
            fprintf(pFile, "%d\n", 0);
        }
        else if((nodes[i]->getState())==NARROW_BAND){
            fprintf(pFile, "%d\n", 1);
        }
        else{
            fprintf(pFile, "%d\n", 2);
        }
    }
    fclose(pFile);
}
