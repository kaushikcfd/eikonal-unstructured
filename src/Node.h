//
// Created by kaushikcfd on 7/1/17.
//
#include "Element.h"
#include <iostream>
#include <vector>

#ifndef GMSH_INTERPRETER_NODE_H
#define GMSH_INTERPRETER_NODE_H
using namespace std;

class Node {
private:
    int index;
    float x; /// The x-coordinate of the node
    float y; /// The y-coordinate of the node
    int noOfNbgElements; /// The number of neighbouring elements for the  node
    vector<Element*> nbgElements; /// This is used to facilitate the operations such as pushing Elemetns, etc.
    vector<float> nbgThetaStart;
    vector<float> nbgThetaEnd;

    int tagState; /// This denotes the state of the node which would help in the FMM algorithm
    int tagAccept; /// This tag helps to located whether the node is accepted or not

public:
    // Constructor
    Node();

    // Updating the value of the data variables
    void updateCoords( float _x, float _y );
    void pushNbgElement( Element* );
    void updateState(int _state);
    void updateAccept(int _accept);

    // Operator
    bool operator==(Node &rhs) const;

    // Index related functions
    void setIndex(int _index);
    int getIndex();

    // Externing the values to the user
    float getX(); /// For the user to read the value of x-coordinate
    float getY(); /// For the user to read the value of y-coordinate
    int getNoOfNbgElements(); /// For the user to read the value of number of nbg Elements
    Element** getNbgElements(); /// For the user to read the value of indices of the nbg Elements

    ~Node();
};

#endif
