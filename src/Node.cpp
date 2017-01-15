//
// Created by kaushikcfd on 7/1/17.
//

#include "Node.h"
#include <cmath>
#include <algorithm>

#define ANGLE(a) (a >= 0 ? (a) : (360.0+a))

Node::Node() {
    x = 0.0;
    y = 0.0;
    noOfNbgElements = 0;
    tagState = 0;
    tagAccept = 0;
}

void Node::updateCoords( float _x, float _y) {
    x = _x;
    y = _y;
}

void Node::pushNbgElement( Element* _nbgElement ) {
    nbgElements.push_back(_nbgElement);

    // Computing the thetaStart and thetaEnd for the current entered element

    float thetaStart, thetaEnd; // Dummy variables to store the values, finally pushed to the data-members
    float x1, y1, x2, y2; // In order to extract the coords of the other two nodes of the same element


    _nbgElement->assigningOtherCoords(*this, x1, y1, x2, y2);
    thetaStart = ANGLE(atan2(y1 - y, x1- x)*180.0/3.1415926);
    thetaEnd = ANGLE(atan2(y2 - y, x2- x)*180.0/3.1415926);
    if(!((((thetaEnd - thetaStart) > 0) && ((thetaEnd - thetaStart) < 180)) ||(((thetaEnd - thetaStart) > -360) && ((thetaEnd - thetaStart) < -180)) ))
        swap(thetaEnd, thetaStart);

    nbgThetaStart.push_back(thetaStart);
    nbgThetaEnd.push_back(thetaEnd);

    noOfNbgElements++;

    return ;
}

void Node::updateState(int _state) {
  tagState = _state ;
}

void Node::updateAccept(int _accept) {
    tagAccept = _accept;
}


bool Node::operator==(Node &rhs) const {
    return ((x == rhs.getX())&&(y == rhs.getY()));
}

float Node::getX() {
    return x;
}

float Node::getY() {
    return y;
}

int Node::getNoOfNbgElements() {
    return noOfNbgElements;
}

Element** Node::getNbgElements() {
    return &nbgElements[0];
}

void Node::setIndex(int _index) {
    index = _index;
}

int Node::getIndex() {
    return index;
}

Node::~Node() {
    nbgElements.clear();
    nbgThetaStart.clear();
    nbgThetaEnd.clear();
}
