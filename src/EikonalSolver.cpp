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
            nodes[i]->updateAccept(ACCEPTED); // This is the condition which is demanded from the boundary condition.
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


float EikonalSolver::calculateCharacteristic(float T, Node* n, float a, float b, float c, float d){
    float phi = 0.0;

    // Collecting information about the given node.
    float F = n->getF();
    float v1 = n->getv1();
    float v2 = n->getv2();

    // Calculating the parameter `G` mentioned in the paper.
    float G = 1 - v1 * (a*T + b) - v2 * (c*T + d);
    
    // Calculating the dyerator and the dxomiator of dy/dx.
    float dy = F * F * (c*T + d) + v2 * G;
    float dx = F * F * (a*T + b) + v1 * G;

    phi = atan2(dy, dx)*(180.0/3.1415926);// Converted the angle from radians to degrees.
    if(phi<0)// Changing the range of phi from (-180, 180] to [0. 360).
        phi+=360;
    
    return phi;
}

bool EikonalSolver::checkCausality(Node *n0, float thetaStart, float thetaEnd, float phi) {
    bool result = false;// Initially initialized to `false` so that it gets toggled once the condition is satisified.
    
    // `phi` is the angle which is the characteristic direction. We need to check the element for if the negative of the
    // characteristic direction lies in the bounds of the given element.
    
    phi+=180;// This gives no problems if the angle is contained within the first 2 quadrants.
    if(phi>=360)
        phi-=360;// This preserved the angle made by the negative of the characteristic angle.

    if(thetaStart > thetaEnd){
        if(((phi>=thetaStart)&&(phi<360))||((phi>=0)&&(phi<thetaEnd)))
            result = true;
    }
    else {
        if(((phi>=thetaStart)&&(phi<thetaEnd)))
            result = true;
    }

    return result;
}

void EikonalSolver::scheme(Node* n) {
    int noNbgElements;
    int initialAccept = n->getAccept(); /// Initial state might be required to know the operations to be performed on the queue.
    float initialSolution = n->getT();

    Element** nbgElements;
    Node *n1=NULL, *n2=NULL;//The other two nodes of the element which is of consideration.
    float *nbgThetaStart = n->getNbgThetaStart(), *nbgThetaEnd = n->getNbgThetaEnd(); // Initialized the arrays with the neighboring thetaStart and thetaEnd, which is the major advantage of self-writen library for mesh interpretation.
    float thetaStart, thetaEnd;
    float phi; // This is the variable which would be used for storing the characteristic angle.

    noNbgElements = n->getNoOfNbgElements(); // Getting the number of nbg. Elements of the node.
    nbgElements = n->getNbgElements(); // Getting all the information about the nbg. Elements.
    
    vector<float> possibleSolutions(noNbgElements, INF);/// This stores the time computed by the solution using that particular element.
    
    vector<bool> calculationIndicator(noNbgElements, false); /// This indicator turns true when a solution is calculated.
    vector<bool> acceptanceIndicator(noNbgElements, false);/// This is an array which stores the information, whether the solution which is arising from the neighboring solution is coming only from accepted nodes.
    vector<bool> causalityIndicator(noNbgElements, false);/// This is the array which stores the information, about whether the characterisitic angle is arising from this direction.
    

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
           // This means both the nodes have their `T` value as `inf`, hence solution cannot be calculated. 
           // So, making all the indicators as FALSE.
           calculationIndicator[j] = false; // Since, no solution is calculated.
           acceptanceIndicator[j] = false; // Since, both of the nodes around are not accepted.
           causalityIndicator[j] = false; // Obviously, since no characteristic direction has been calculated.           
        }
        else if((n1->getAccept() == NOT_ACCEPTED ) && (n2->getAccept() != NOT_ACCEPTED )) {
            a =  - (m12/((m22*m11 - m12*m21)*(lengthBX))) ;
            b =    (((n2->getT())*m12)/((m22*m11 - m12*m21)*(lengthBX))) ;
            c =    (m11/((m22*m11 - m12*m21)*(lengthBX))) ;
            d =  - (((n2->getT())*m11)/((m22*m11 - m12*m21)*(lengthBX))) ;
            
            // Solving the quadratic solution as the coefficients are known.
            t = solution(F, v1, v2, a, b, c, d); /// The solution of the current triangle has been found.
            possibleSolutions[j] = t; // Noting the solution obtained.
            
            // Calculating the characteristic angle arising due to the solution.
            phi = calculateCharacteristic(t, n, a, b, c, d);

           // Assigning the values to thetaStart and thetaEnd.
           
            thetaStart = nbgThetaStart[j];
            thetaEnd = nbgThetaEnd[j];
            

            // Assigning the indicators appropriately.
            calculationIndicator[j] = true; // As we have some finite value of `T` stored here.
            acceptanceIndicator[j] = false; // Obviously as one of the solutions is coming from FAR_AWAY node.
            causalityIndicator[j] = checkCausality(n, thetaStart, thetaEnd, phi);

        }
        else if((n1->getAccept() != NOT_ACCEPTED ) && (n2->getAccept() == NOT_ACCEPTED )) {
            // compute the value but note down.
            a =  (m22/((m22*m11 - m12*m21)*(lengthAX)));                 
            b = -(((n1->getT())*m22)/((m22*m11 - m12*m21)*(lengthAX)));
            c = -(m21/((m22*m11 - m12*m21)*(lengthAX)));                
            d =  (((n1->getT())*m21)/((m22*m11 - m12*m21)*(lengthAX))); 


            // Solving the quadratic solution as the coefficients are known.
            t = solution(F, v1, v2, a, b, c, d); /// The solution of the current triangle has been found.
            possibleSolutions[j] = t; // Noting the solution obtained.
            
            // Calculating the characteristic angle arising due to the solution.
            phi = calculateCharacteristic(t, n, a, b, c, d);

           // Assigning the values to thetaStart and thetaEnd.
           
            thetaStart = nbgThetaStart[j];
            thetaEnd = nbgThetaEnd[j];
            

            // Assigning the indicators appropriately.
            calculationIndicator[j] = true; // As we have some finite value of `T` stored here.
            acceptanceIndicator[j] = false; // Obviously as one of the solutions is coming from FAR_AWAY node.
            causalityIndicator[j] = checkCausality(n, thetaStart, thetaEnd, phi);

            
        }
        else {

            a =  (m22/((m22*m11 - m12*m21)*(lengthAX)))                - (m12/((m22*m11 - m12*m21)*(lengthBX))) ;
            b = -(((n1->getT())*m22)/((m22*m11 - m12*m21)*(lengthAX))) + (((n2->getT())*m12)/((m22*m11 - m12*m21)*(lengthBX))) ;
            c = -(m21/((m22*m11 - m12*m21)*(lengthAX)))                + (m11/((m22*m11 - m12*m21)*(lengthBX))) ;
            d =  (((n1->getT())*m21)/((m22*m11 - m12*m21)*(lengthAX))) - (((n2->getT())*m11)/((m22*m11 - m12*m21)*(lengthBX))) ;


            // Solving the quadratic solution as the coefficients are known.
            t = solution(F, v1, v2, a, b, c, d); /// The solution of the current triangle has been found.
            possibleSolutions[j] = t; // Noting the solution obtained.
            
            // Calculating the characteristic angle arising due to the solution.
            phi = calculateCharacteristic(t, n, a, b, c, d);

           // Assigning the values to thetaStart and thetaEnd.
           
            thetaStart = nbgThetaStart[j];
            thetaEnd = nbgThetaEnd[j];

            // Assigning the indicators appropriately.
            calculationIndicator[j] = true; // As we have some finite value of `T` stored here.
            
            // Setting the acceptance indicator using the below if-else condition. 
            // The node is marked `true` for acceptanceIndicator iff both the neighboring nodes are accepted.
            if((n1->getAccept() == ACCEPTED ) && (n2->getAccept() == ACCEPTED ))
                acceptanceIndicator[j] = true;
            else
                acceptanceIndicator[j] = false;
            
            causalityIndicator[j] = checkCausality(n, thetaStart, thetaEnd, phi);// Using the function to check the causality condition.
            
                        
            
        }
    }
    /// Now, we are checking for the column which correctly satisfies all the conditions and its minimum is taken to compute the value of `T`
    float minTime = INF;
    int indexOfMinElement = -1; // This stores the element corresponding to the minimum.

    for(int j = 0; j < noNbgElements; j++ ){
        if( (acceptanceIndicator[j] == true) && (causalityIndicator[j] == true) && (possibleSolutions[j] < minTime) ) {
            // These statements will be entered only if both the neighboring nodes are ac
            minTime = possibleSolutions[j];
            indexOfMinElement = j;
        }
    }

    

    if(indexOfMinElement > -1){
        // This means that we had one element with the minimum time which satisifies all the conditions and hence it must be accepted.
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
            if( (possibleSolutions[j] < minTime) ) {// WARNING: WE NEED SOME CONDITION HERE. ASK Baskar Sir.
                // Here in the condition, we can expect one more condition about cardinality
                minTime = possibleSolutions[j];
            }
        }
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

    set<Node*> set_new; /// The neighboring nodes of the currentNode which are needed to be recomputed are either belonging to NarrowBand or are belonging to Far-Away points, and at once all the nodes are scanned through these sets and their `T` is recomputed.

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
            /// We are currently looping over all the nodes which are the neighboring nodes of THE TOP NODE of the heap.
            nbgElements[j]->assigningOtherNodes(currentNode, n1, n2);/// n1, n2 now contain the address of the current neighboring element.
            
            /// Inserting both the nodes into the set_new. They are only inserted into the set_new if they are not ALIVE.
            if(n1->getState()!=ALIVE)
                set_new.insert(n1);
            if(n2->getState()!=ALIVE)
                set_new.insert(n2);
        
        } /// Now we have marked out all our nodes which are eligible to be recomputed.

        while(!set_new.empty()){
            // Run recompute for each of the Narrow Band points
            // This should update the T which is present in them already
            Node* n =  *set_new.begin();
            scheme(n);
            n->updateState(NARROW_BAND);
            // Write the code here to run the re-compute.
            set_new.erase(set_new.begin());/// As the corresponding node has been filled in with the updated value of T, so now removing it from the set.
        }

    }

    return 0;
}
