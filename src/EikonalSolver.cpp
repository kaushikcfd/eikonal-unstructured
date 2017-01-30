#include "EikonalSolver.h"


EikonalSolver::EikonalSolver(Mesh2D* _mesh, function<float(float, float)> _F, function<float(float, float)> _v1, function<float(float, float)> _v2) {
    mesh = _mesh;
    F = _F;
    v1 = _v1;
    v2 = _v2;
    
    // The functions have been intialized, now using these functions to set the values for each nodes.

    float x, y;
    Node** nodes = mesh->nodes;
    int noNodes = mesh->getNoOfNodes();
    for( int i = 0; i < noNodes; i++) {
        x = nodes[i]-> getX();
        y = nodes[i]-> getY();
        nodes[i]->setF(F(x, y)); // Setting the wave speed for the $i^th$ node.
        nodes[i]->setv1(v1(x, y)); // Setting the medium speed (x-direction) for the $i^th$ node.
        nodes[i]->setv2(v2(x, y)); // Setting the medium speed (y-direction) for the $i^th$ node.
        if(nodes[i]->getState() == NARROW_BAND)
            narrowBandHeap.push(nodes[i]);
    }
}

void EikonalSolver::plotStates(string outputFile="States.dat") {
    FILE* pFile;
    pFile = freopen(outputFile.c_str(), "w", stdout);
    Node** nodes = mesh->nodes;
    int noNodes = mesh->getNoOfNodes();
    for( int i = 0; i < noNodes; i++) {
        if((nodes[i]->getState())==ALIVE){
            fprintf(pFile, "%d\n", 0);
        }
        else if((nodes[i]->getState())==NARROW_BAND){
            fprintf(pFile, "%d\n", 1);
        }
        else{
            fprintf(pFile, "%d\n", 2);
        }
    }
    fclose(pFile);
}


/**
 * @brief This function is used to refresh the whole heap. This maybe a time consuming process but yes, there is no escaping!
 */
void EikonalSolver::refreshHeap() {
    queue<Node*> q;

    // Populating the queue `q`, with the elements of the narrowBandHeap.
    while(!narrowBandHeap.empty()){
        q.push(narrowBandHeap.top());
        narrowBandHeap.pop();
    }

    // Now populating the narrowBandHeap, with the queue `q`.
    while(!q.empty()) {
        narrowBandHeap.push(q.front());
        q.pop();
    }

    // Successfully rearranged the priority queue according to the priorities. Exiting!
    return ;
}

void EikonalSolver::scheme(Node* n) {
    int noNbgElements;
    int initialAccept = n->getAccept(); /// Initial state might be required to know the operations to be performed on the queue.
    float initialSolution = n->getT();

    Element** nbgElements;
    Node *n1=NULL, *n2=NULL;//The other two nodes of the element which is of consideration.

    noNbgElements = n->getNoOfNbgElements(); // Getting the number of nbg. Elements of the node.
    nbgElements = n->getNbgElements(); // Getting all the information about the nbg. Elements.
    
    vector<float> possibleSolutions(noNbgElements, INF);/// This stores the time computed by the solution using that particular element.
    vector<int> solutionStates(noNbgElements, 0);/// This stores the state of the solution that has been computed using that element. The legend for the solutionStates has been mentioned ahead.
    /**
     * `i` is the element from which the solution has been computed.
     * If possibleSolution[i] == 0 : The solution has the value INF, and must not be used in any computation ahead.
     * If possibleSolution[i] == 1 : One of the nodes of the element is FAR AWAY and this must not be acepted but the
     * If possibleSolution[i] == 2: Both the nodes have meaningful value, and this node can be accepted. 
     */

    float m11, m12, m21, m22;
    float x, y, a1 , a2, b1, b2; // Using the notation used by Dahiya & Baskar, Characteristic Fast Marcing Method on Triangular Grids for the generalized eikonal equation in moving media.
    float lengthAX, lengthBX; // These are the lengths of AX and BX. 

    float a, b, c, d; // The coefficients to the quadratic equation.
    float F, v1, v2; // Properties of the node
    float t;

    /**Initializing the variables x and y**/
    x  = n->getX();
    y  = n->getY();
    /**Initializing the properties of the node.**/
    F  = n->getF();
    v1 = n->getv1();
    v2 = n->getv2(); 

    /**Getting the parameters to set up the scheme.**/
    a1 = n1->getX(); a2 = n1->getY();
    b1 = n2->getX(); b2 = n2->getY();

    lengthAX = sqrt((x-a1)*(x-a1) + (y-a2)*(y-a2));
    lengthBX = sqrt((x-b1)*(x-b1) + (y-b2)*(y-b2));

    m11 = (x - a1) / lengthAX;
    m12 = (y - a2) / lengthAX;
    m21 = (x - b1) / lengthBX;
    m22 = (x - b2) / lengthBX;

    /**Looping through all the Nbg. Elements**/
    for(int j=0; j<noNbgElements; j++) {
        nbgElements[j]->assigningOtherNodes(n, n1, n2); // n1, n2 contains the information about the other two nodes of the element.
        if((n1->getAccept() == NOT_ACCEPTED ) && (n2->getAccept() == NOT_ACCEPTED )){
            solutionStates[j] = 0; // As the information is not taken from the correct nodes.
        }
        else if((n1->getAccept() == NOT_ACCEPTED ) && (n2->getAccept() != NOT_ACCEPTED )) {
            a =  - (m12/((m22*m11 - m12*m21)*(lengthBX))) ;
            b =    (((n2->getT())*m12)/((m22*m11 - m12*m21)*(lengthBX))) ;
            c =    (m11/((m22*m11 - m12*m21)*(lengthBX))) ;
            d =  - (((n2->getT())*m11)/((m22*m11 - m12*m21)*(lengthBX))) ;
            
            // Solving the quadratic solution as the coefficients are known.
            t = solution(F, v1, v2, a, b, c, d); /// The solution of the current triangle has been found.
            possibleSolutions[j] = t; // Noting the solution obtained.
            
            // Setting the solution state approppriately.
            solutionStates[j] = 1;
        }
        else if((n1->getAccept() != NOT_ACCEPTED ) && (n2->getAccept() == NOT_ACCEPTED )) {
            // compute the value but note down.
            a =  (m22/((m22*m11 - m12*m21)*(lengthAX)));                 
            b = -(((n1->getT())*m22)/((m22*m11 - m12*m21)*(lengthAX)));
            c = -(m21/((m22*m11 - m12*m21)*(lengthAX)));                
            d =  (((n1->getT())*m21)/((m22*m11 - m12*m21)*(lengthAX))); 
            
            // Calculating the solution to the quadratic equation.
            t = solution(F, v1, v2, a, b, c, d); /// The solution of the current triangle has been found.
            possibleSolutions[j] = t; // 
            
            // Setting the solution state appropriately 
            solutionStates[j] = 1;
        }
        else {

            a =  (m22/((m22*m11 - m12*m21)*(lengthAX)))                - (m12/((m22*m11 - m12*m21)*(lengthBX))) ;
            b = -(((n1->getT())*m22)/((m22*m11 - m12*m21)*(lengthAX))) + (((n2->getT())*m12)/((m22*m11 - m12*m21)*(lengthBX))) ;
            c = -(m21/((m22*m11 - m12*m21)*(lengthAX)))                + (m11/((m22*m11 - m12*m21)*(lengthBX))) ;
            d =  (((n1->getT())*m21)/((m22*m11 - m12*m21)*(lengthAX))) - (((n2->getT())*m11)/((m22*m11 - m12*m21)*(lengthBX))) ;
            
            // Calculating the solution to the quadratic equation.
            t = solution(F, v1, v2, a, b, c, d); /// The solution of the current triangle has been found.
            possibleSolutions[j] = t; // Noting down the solution computed using the particular element. 
            
            // Setting the solution states with the conditions mentioned alongside.
            if((n1->getAccept() == ACCEPTED ) && (n2->getAccept() == ACCEPTED )) 
                solutionStates[j] = 2; // State = 2; as everything is perfectly fine!
            else
                solutionStates[j] = 1;// One of the node may be AMBIGUOUS and hence this node cannot be directly accepted!
            
            // Here add the cardinality check, which may again affect the solutionState[j], Not doing it for now!
        }
    }
    /// Now, we are checking for the column which correctly satisfies all the conditions and its minimum is taken to compute the value of `T`
    float minTime = INF;
    int indexOfMinElement = -1; // This stores the element corresponding to the minimum.

    for(int j = 0; j < noNbgElements; j++ ){
        if( (solutionStates[j]==2) && (possibleSolutions[j] < minTime) ) {
            // Here in the condition, we can expect one more condition about cardinality
            minTime = possibleSolutions[j];
            indexOfMinElement = j;
        }
    }

    

    if(indexOfMinElement > -1){
        // This means that we had one element with the minimum time which satisifies all the conditions and hence it must be accepted.
        n->updateState(NARROW_BAND); // Adding the node to the Narrow Band 
        n->updateAccept(ACCEPTED);/// Accepting the node.
        n->setT(minTime);/// Updating the `T` of the node.
    
        if(initialAccept != ACCEPTED){
            // Inititally the `Node n` was not a part of the NarrowBand and hence this must be pushed into the narrowBandHeap
            narrowBandHeap.push(n);  
        }
        else {
            // This must a narrowBand Node, and hence the edits have been done in the heap. Now the heap must be updated to sort once again.
            if(minTime!=initialSolution)// If the new Solution conicided with the earlier Solution, no need to make changes in the narrowBandHeap.
                refreshHeap();
        }
    }
    else {
        /// The node cannot be accepted, but if atleast we could set the `T` value from one of the solutions, even though it may be not accepted, then it is fine.
        for(int j = 0; j < noNbgElements; j++){
            if( (solutionStates[j]==1) && (possibleSolutions[j] < minTime) ) {
                // Here in the condition, we can expect one more condition about cardinality
                minTime = possibleSolutions[j];
            }
        }
        n->updateState(NARROW_BAND); // Making the node Narrow_Band
        n->updateAccept(AMBIGUOUS);
        n->setT(minTime);
    }

    
    return ;
}

/**
 * @brief This is the function which calls the schemeF(), and returns the computed new time. This calls the functions scheme, and does most of the computations
 * @param n The `Node* n` contains the address of the node whose arrival time has to be computed.
 */

void EikonalSolver::recompute(Node* n) {


    return ;
}


/**
 * @brief Solves the quadratic equation for a node, and calculates the Tij at that point.
 
 * @param F The wave speed at the node.
 * @param vx The speed of the medium in the x-direction.
 * @param vy The speed of the medium in the y-direction.
 * @param a, b, c, d : These paramaeters have been explained in the paper D. Dahiya, et. al.
 * @return The computed wave approach time at the node
 */

float EikonalSolver::solution(float F, float vx, float vy, float a, float b, float c, float d)
{
    double v[2],FF;
    double Sol[2];
    double dummy1,dummy2,T;
    double coef1,coef2,coef3,discr;

    v[0]=vx;v[1]=vy;
    coef1=0.0;coef2=0.0;coef3=0.0;
    FF=F*F;

    dummy1=FF-v[0]*v[0];
    dummy2=FF-v[1]*v[1];

    coef1=dummy1*a*a-2.0*v[0]*v[1]*a*c+dummy2*c*c;
    coef2=2.0*(dummy1*a*b+v[0]*a-v[0]*v[1]*(a*d+b*c)+v[1]*c+dummy2*c*d);
    coef3=dummy1*b*b+2.0*v[0]*b-2.0*v[0]*v[1]*b*d+2.0*v[1]*d+dummy2*d*d-1.0;

    discr=coef2*coef2-4.0*coef1*coef3;
    if(discr<0)
    {
//      printf("Error: Solution becomes complex at (i,j)= (%d,%d). \n",i,j);
//      exit(0);
        return(INF);
    }

    if(coef1==0)
    {
//        printf("coef1 = 0 at (%d, %d) with dummy1 = %lf, dummy2 = %lf  a=%lf c=%lf.\n",i,j,dummy1,dummy2,a,c);
        exit(0);
        return(INF);
    }

    discr=sqrt(discr);
    T=2.0*coef1;

    Sol[0]=(-1*coef2+discr)/T;
    Sol[1]=(-1*coef2-discr)/T;
    /**Choosing the maximum**/
    T=((Sol[0]>Sol[1])?Sol[0]:Sol[1]);

    return(T);
} // end of solution

/**
 * @brief This is the main function which has to be called which ensures that all the rest of the functions are called appropriately through this function.
 * @return The state of the computation. `0` indicates  that everything was computed stably.
 */

int EikonalSolver::solve() {

    Node* currentNode;/// This is the node which currently has the minimum arrival time in the given heap.

    Node *n1, *n2; /// The other two nodes of a given neighboring element.

    set<Node*> set_N, set_F; /// The neighboring nodes of the currentNode which are needed to be recomputed are either belonging to NarrowBand or are belonging to Far-Away points.
    /// If they are belonging to the `set N`  they are pushed into it, and likewise for `set F`.
    /// And at once all the nodes are scanned through these sets and their `T` is recomputed.

    int noNbgElements;
    Element** nbgElements; 

    // Taking the minimum of all the current nodes in the narrow band
    while(!narrowBandHeap.empty()){

        currentNode = narrowBandHeap.top();// Now it is ensured that the node with the minimum is considered.

        narrowBandHeap.pop(); // Removed from the heap.

        currentNode->updateState(ALIVE);/// Marking the node as alive.
        
        noNbgElements = currentNode->getNoOfNbgElements();
        nbgElements = currentNode->getNbgElements();

        for(int j=0; j<noNbgElements; j++) {
            nbgElements[j]->assigningOtherNodes(currentNode, n1, n2);
            /// We are currently looping over all the nodes which are the neighboring nodes of THE TOP NODE of the heap.

            // Now, the aim is to bifurcate all the neighboring nodes so that they can be treated separately.  
            
            /**Observing the state of `n1` **/
            if((n1->getState()) == NARROW_BAND){
                set_N.insert(n1);
            }
            else if((n1->getState()) == FAR_AWAY){
                set_F.insert(n1);
            }

            /**Observing the state of `n2` **/
            if((n2->getState()) == NARROW_BAND){
                set_N.insert(n2);
            }
            else if((n2->getState()) == FAR_AWAY){
                set_F.insert(n1);
            }            
        } /// Now we have marked out all our nodes which are eligible to be recomputed.

        while(!set_N.empty()){
            // Run recompute for each of the Narrow Band points
            // This should update the T which is present in them already
            Node* n =  *set_N.begin();
            // Write the code here to run the re-compute.
            set_N.erase(set_N.begin());/// As the corresponding node has been filled in with the updated value of T, so now removing it from the set.
        }

        while(!set_F.empty()) {
            // Main purpose of this loop is to ensure that the node which was initially far away, should now be pushed into the Narrow Band.
            Node* n = *set_F.begin();
            // Write the code here to compute the `T` for the node.
            // First check whether the recompute was successful, only then add the node to the narrow band heap.
            n->updateState(NARROW_BAND); // Assuming everything went fine.

            set_F.erase(set_F.begin());/// This step should always be there, as the node should always be removed from the Far Away region.
        }
    }

    return 0;
}
