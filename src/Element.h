//
// Created by kaushikcfd on 7/1/17.
//


#ifndef GMSH_INTERPRETER_ELEMENT_H
#define GMSH_INTERPRETER_ELEMENT_H
#include <iostream>
class Node; // Forward declaration
/**
 * Currently the element is only able to handle triangular elements
 * */

class Element {
private:
    int index;
    Node** nodes; // Currently, there is only support for triangular elements

    int whichNodeOfElement(const Node &_node);

public:
    Element();
    void setNode1(Node* _node);
    void setNode2(Node* _node);
    void setNode3(Node* _node);
    void setIndex(int _index);
    int getIndex();
    Node* getNode(int node_index);
    void assigningOtherCoords(const Node &_node, float &x1, float &y1, float &x2, float &y2);

    ~Element();

};
#endif //GMSH_INTERPRETER_ELEMENT_H
