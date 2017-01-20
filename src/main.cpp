/**
 * This file is the main driver file.
 * I plan to add all the settings for the for the functions in this file.
 * The definition of the functions for F(x,y), v1(x, y), v2(x, y) would all be present here.
 *  
 */

#include "main.h"
#define FILE_PATH "square-structured.msh"

/**
 * This function serves as an input function for wave speed. In all the papers this is denoted as `F(x, y)`
 * @param  x The x-coordinate for the function
 * @param  y The y-coordinate for the function
 * @return   Speed of the wave in the medium of interest
 */
float waveSpeed(float x, float y) {
    return 1.0;
}

/**
 * This function is used to give the input function for the speed in the x-direction. In the paper this is denoted by v1(x, y)
 * @param  x The x-coordinate for the function
 * @param  y The y-coordinate for the function
 * @return   The medium speed in the x-direction
 */
float mediumSpeed1(float x, float y) {
    return 0.0;
}

/**
 * This function is used to give the input function for the speed in the y-direction. In the paper this is denoted by v2(x, y)
 * @param  x The x-coordinate for the function
 * @param  y The y-coordinate for the function
 * @return   The medium speed in the y-direction
 */
float mediumSpeed2(float x, float y) {
    return 0.0;
}

// The above functions define the problem.


//----- Functions to assign the initial wavefront

/**
 * @brief      To check whether a satisifies the condition of the intial wavefront.
 *
 * @param[in]  n	The node which is to be checked whether it satisifies the constraints of the initial wavefront.     
 *
 * @return     True if alive, False otherwise.
 */
bool isAlive_planar(Node* n) {
	/// The function is constructed in such a way that for now it is being asumed that we have a planar wavefront at x=0;
	return ((n->getX()) == 0); 
}

void setAlivePoints(Mesh2D* mesh, function<bool(Node*)> isAlive) {
	int noNodes = mesh->getNoOfNodes(), noNbgElements, i;
	Node **nodes = mesh->nodes;
	Element** nbgElements;
	Node *n1=NULL, *n2=NULL;//The other two nodes which have to be made narrow band points
	// TODO: Start from here, all you need to do is get the nbgElements and then in each of the nbgElements check for the 
	for(i=0; i<noNodes; i++) {
		if(isAlive(nodes[i])){
			nodes[i]->updateState(ALIVE);// Updated the state
			nodes[i]->setT(0.0);// Update the time to be zero ;-)
			nodes[i]->updateAccept(ACCEPTED);			              
            //--- Now, adding the neighbouring nodes of the Alive nodes to the narrow band.
             noNbgElements = nodes[i]->getNoOfNbgElements();
             nbgElements = nodes[i]->getNbgElements();
             for(int j=0; j<noNbgElements; j++) {
             	nbgElements[j]->assigningOtherNodes(nodes[i], n1, n2);
             	if(!isAlive(n1)){
             		n1->updateState(NARROW_BAND);
             		n1->setT(100.0);//// TODO: Find a way to calculate this.
             	}
             	if(!isAlive(n2)){
             		n2->updateState(NARROW_BAND);
             		n2->setT(100.0);//// TODO: Find a way to calculate this.             		
             	}
             }
		}
	}
	return;
}

/**---------Starting the driver, and calling all the functions-------**/
int main() {
    Mesh2D mesh;
    mesh.readFromFile(FILE_PATH);

    /**--------------------Initializing-------------------------**/

     


    return 0;
}































/**Use this code to check whether the mesh is generated correctly**/

/*FILE* pFile = freopen("nbgTheta.dat", "w", stdout);

    for(int i=0; i < mesh.getNoOfNodes(); i++ ) {
        fprintf(pFile, "Node %d\n", i);
        Node *currentNode = mesh.nodes[i];
        for(int j = 0; j < currentNode->getNoOfNbgElements(); j++ ) {
            fprintf(pFile, "%6.3f\t", currentNode->getNbgThetaStart()[j]);
        }
        fprintf(pFile, "\n");
        for(int j = 0; j < currentNode->getNoOfNbgElements(); j++ ) {
            fprintf(pFile, "%6.3f\t", currentNode->getNbgThetaEnd()[j]);
        }
        fprintf(pFile, "\n------------------------------------------------------\n");
    }
    fclose(pFile);
    mesh.write("Nodes.dat", "Elements.dat");
    */
   // The above operation was to read from the given .msh file & then output the neighbouring elements and the angles formed by the elements w.r.t. a node
