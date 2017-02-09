//
// Created by kaushikcfd on 7/1/17.
//
#include "Element.h"
#include <iostream>
#include <vector>
#include <limits>

#ifndef GMSH_INTERPRETER_NODE_H
#define GMSH_INTERPRETER_NODE_H

#define INF 34028234663852885981170418348451692544.00 // This is the maximum value taken by a float divided by 10

#define ALIVE 100
#define NARROW_BAND 200
#define FAR_AWAY 300

#define ACCEPTED 1000
#define AMBIGUOUS 2000
#define NOT_ACCEPTED 3000

/**
 * The different possible accept solutions are explained below:
 * ACCEPTED: The `T` value is != INF, and the value has come from accepted nodes.
 * AMBIGUOUS: The `T` value is !=INF, and the value has come from atleast one NOT_ACCEPTET/AMBIGUOUS node(s)
 * NOT_ACCEPTED: The `T` value is = INF, this is equivalent to a FAR_AWAY node.
 */ 

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
    

    float T; /// The minimum arrival time; this is specific to the question
    float F; /// The wave speed in the given medium, at this node
    float v1; /// The medium speed in the x-direction 
    float v2; /// The medium speed in the y-direction

    int timesRecomputed; /// This variable stores the number of times a variables value is recomputed.

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
    float getX() const ; /// For the user to read the value of x-coordinate
    float getY(); /// For the user to read the value of y-coordinate
    int getNoOfNbgElements(); /// For the user to read the value of number of nbg Elements
    Element** getNbgElements(); /// For the user to read the value of indices of the nbg Elements

    // ThetaStart and ThetaEnd related routines
    float* getNbgThetaStart();
    float* getNbgThetaEnd();

    // `T` related routines
    void setT(float _T); // For the user to change the value in `T`
    float getT(); // For the user to get `T`

    // `F, v` related routines
    void setF(float _F);
    void setv1(float _v1);
    void setv2(float _v2);
    float getF();
    float getv1();
    float getv2();
    
    // Getting the tags i.e. tagState & tagAccept
    int getState();
    int getAccept();
    int getTimesRecomputed();

    // Destructor
    ~Node();
};

#endif
