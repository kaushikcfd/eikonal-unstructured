#ifndef GMSH_INTERPRETER_EIKONALSOLVER_H
#define GMSH_INTERPRETER_EIKONALSOLVER_H

#include <queue>
#include <functional>
#include <string>
#include <set>
#include <queue>
#include <cmath>

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
    
    // Debugging data-structures
    vector<float> y_cord, v_cord;// This will be used to check the value of the wavespeed that will be assigned.
    // Debuggin functions
    void checkWaveSpeed();// This will print the wavespeed values to a file.

    // Private functions
    float solution(float F, float vx, float vy, float a, float b, float c, float d); // Computes the solution to the quadratic eqn.
    void scheme(Node* n); /// This is the function which actually loops through all the neighboring elements and computes the value of T.
    void refreshHeap(); /// This function would help us construct the heap again from scratch.
    float calculateCharacteristic(float T, Node* n, float a, float b, float c, float d);/// This is the function which calculates the causality direction, which is represented as Phi in the paper of reference.
    bool checkCausality(Node *n0, float thetaStart, float thetaEnd, float phi);/// This is the function which tests whether the causality conditions is met. 

public:
    EikonalSolver(Mesh2D* _mesh, function<float(float, float)> _F, function<float(float, float)> _v1, function<float(float, float)> _v2);
    void plotStates(string outputFile);
    int solve();

};

#endif
