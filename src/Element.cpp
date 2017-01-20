//
// Created by kaushikcfd on 7/1/17.
//

#include "Element.h"
#include "Node.h"
#include <cstdio>

Element::Element() {
    nodes = new Node*[3]; /// Assigning only three nodes because, currently considering only unstructured mesh with triangular elements
}

void Element::setNode1(Node* _node) {
    nodes[0] = _node;
    return ;
}

void Element::setNode2(Node* _node) {
    nodes[1] = _node;
    return ;
}

void Element::setNode3(Node* _node) {
    nodes[2] = _node;
    return ;
}

Node* Element::getNode (int node_index) {
    return (nodes[node_index]); /// Returned the address of required Node
}

void Element::setIndex(int _index) {
    index = _index;
}

int Element::getIndex() {
    return index;
}

int Element::whichNodeOfElement(const Node* const _node) {
    for(int i = 0; i < 3; i++) {
        if(_node == nodes[i])
            return i;
    }
    printf("Not found Node at the end of the element. Exit!\n");
    return -1;
}

void Element::assigningOtherCoords(const Node* const _node, float &x1, float &y1, float &x2,
                                  float &y2) {
    int count = 0;
    int node_indexInElement = whichNodeOfElement(_node);
    for (int i=0; i<3; i++) {
        if(node_indexInElement != i){
            if(count == 0){
                x1 = nodes[i]->getX();
                y1 = nodes[i]->getY();
                ++count;
            }
            else {
                x2 = nodes[i]->getX();
                y2 = nodes[i]->getY();
                return ;
            }
        }
    }
}

/**
  * @brief      Used to get the other two nodes of the elements and the inputs &n1. &n2 would be updated accordingly.
  *
  * @param[in]  _node  The node whose neighbouring nodes are to be found out.
  * @param[out] n1     The adress of the first node would be stored in this pointer.
  * @param[out]	n2	   The adress of the other node would be stored in this pointer.
  */
void Element::assigningOtherNodes(Node* _node, Node* n1, Node* n2) {
	int count = 0;
    int node_indexInElement = whichNodeOfElement(_node);
    for (int i=0; i<3; i++) {
        if(node_indexInElement != i){
            if(count == 0){
                n1 = nodes[i];
                ++count;
            }
            else {
                n2 = nodes[i];
                return ;
            }
        }
    }
}

Element::~Element() {
    delete [] nodes;
}
