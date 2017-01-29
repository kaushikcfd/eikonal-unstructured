//
// Created by kaushikcfd on 13/12/16.
//

#include "Mesh2D.h"

Mesh2D::Mesh2D() {
    noNodes = 0;
    noElements = 0;
}

/**
 * This function is used to allocate memory for the `nodes` array once the number of nodes is known.
 */
void Mesh2D::allocateNodes() {
    nodes = new Node*[noNodes];
    for(int i=0; i<noNodes; i++)
        nodes[i] = new Node;
    return ;
}

/**
 * This function is used to allocate memort for the `elements` array once the number of elements is known.
 */
void Mesh2D::allocateElements() {
    elements = new Element*[noElements];
    for(int i=0; i<noElements; i++)
        elements[i] = new Element;
    return ;
}


void Mesh2D::readFromFile(string filename) {
    FILE* pFile = freopen(filename.c_str(), "r", stdin);
    string s;
    int i, node_no, elm_no, elm_type, tag1, tag2, tag3;
    float node_x, node_y, node_z;// This will always be zeros, as this is the 2D case
    int node_index1, node_index2, node_index3;

    for(i=0; i<4; i++)// Omitting the top 4 lines, which describe the Mesh Format
        getline(cin, s);

    scanf("%u", &noNodes);// Reading the number of nodes involved in the mesh
    allocateNodes();

    for(i=0; i<noNodes; i++) {
        scanf("%d %f %f %f", &node_no, &node_x, &node_y, &node_z);// Reading the x,y,z of each node.
        nodes[i]->updateCoords(node_x, node_y);
        nodes[i]->setIndex(i);
    }

    for(i=0; i<3; i++)// Omitting the 2 lines, which describe the syntax of a .msh file
        getline(cin, s);

    scanf("%u", &noElements);// Reading the number of elements that are to be read.

    allocateElements();

    for(i=0;i<noElements;i++) {
        scanf("%d %d",&elm_no, &elm_type);
        if(elm_type == 2){
            // Setting the index of the element
            elements[i]->setIndex(i);
            // Reading the three nodes of each elements.
            scanf("%d %d %d %d %d %d", &tag1, &tag2, &tag3, &node_index1, &node_index2, &node_index3);
            elements[i]->setNode1(nodes[node_index1-1]);
            elements[i]->setNode2(nodes[node_index2-1]);
            elements[i]->setNode3(nodes[node_index3-1]);

            for (int j=0; j<3; j++) {
                elements[i]->getNode(j)->pushNbgElement(elements[i]);
            }
        }

        else{
            // Omitting the line, as other elements are of no use.
            getline(cin, s);
            --i;
            --noElements;
        }
    }
    fclose(pFile);
    return ;
}

int Mesh2D::getNoOfNodes(){
    return noNodes;
}
int Mesh2D::getNoOfElements(){
    return noElements;
}

void Mesh2D::write(string nodeFile, string elementFile) {


    int i;
    FILE* pFile = freopen(nodeFile.c_str(), "w", stdout);
    fprintf(pFile, "%d\n", noNodes);
    for (i=0; i<noNodes; i++) {
        fprintf(pFile, "%6.3f\t%6.3f\n", nodes[i]->getX(), nodes[i]->getY());
    }
    fclose(pFile);

    pFile = freopen(elementFile.c_str(), "w", stdout);
    fprintf(pFile, "%d\n", noElements);
    for(i=0; i<noElements; i++) {
        fprintf(pFile, "%d\t%d\t%d\n", elements[i]->getNode(0)->getIndex(), elements[i]->getNode(1)->getIndex(), elements[i]->getNode(2)->getIndex());
    }
    fclose(pFile);

    return ;
}

Mesh2D::~Mesh2D() {
    for (int i=0; i < noNodes; i++)
        delete nodes[i];
    delete [] nodes;
    for (int i=0; i < noElements; i++)
        delete elements[i];
    delete [] elements;
}
