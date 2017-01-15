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

int Element::whichNodeOfElement(const Node &_node) {
    for(int i = 0; i < 3; i++) {
        if(&_node == nodes[i])
            return i;
    }
    printf("Not found Node at the end of the element. Exit!\n");
    return -1;
}

void Element::assigningOtherCoords(const Node &_node, float &x1, float &y1, float &x2,
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

Element::~Element() {
    delete [] nodes;
}
