#ifndef GMSH_INTERPRETER_EIKONALSOLVER_H
#define GMSH_INTERPRETER_EIKONALSOLVER_H

#include <queue>
#include <functional>
#include <string>
#include "Mesh2D.h"
#include "Node.h"

using namespace std;

/**
 * This class is declared to specialy handle the priority_queue.
 * This serves as an input while declaring the cardinal priority_queue, which is being used as a heap.
 */
class Compare {
public:
    bool operator()(Node* n1, Node *n2) {
        return ((n1->getT()) > (n2->getT()));
    }
};


/**
 * Purpose of this class is to solve to the eikonal PDE over unstructures domains
 */
class EikonalSolver {
private:
    Mesh2D *mesh;

    // Functions that define the problem
    function<float(float, float)> F; /// The wave speed in the medium
    function<float(float, float)> v1; /// The medium speed in the x-direction
    function<float(float, float)> v2; /// The medium speed in the y-direction

    priority_queue<Node*, vector<Node*>, Compare> narrowBandHeap;/// The heap controlling the narrow band points

public:
    EikonalSolver(Mesh2D* _mesh, function<float(float, float)> _F, function<float(float, float)> _v1, function<float(float, float)> _v2);
    void readInitialFront(string inputFile);
    void plot(string outputFile);

};

#endif
