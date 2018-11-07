//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) INF/UFG 2018. All rights reserved.                                                       
//                                                                                                        
// License is granted to copy, to use, and to make and to use derivative                                  
// works for research and evaluation purposes, provided that INF/UFG is                                   
// acknowledged in all documentation pertaining to any such copy or                                       
// derivative work. INF/UFG grants no other licenses expressed or                                             
// implied.                                                                                               
//                                                                                                        
// INF/UFG MAKES NO REPRESENTATIONS CONCERNING EITHER THE                                                 
// MERCHANTABILITY OF THIS SOFTWARE OR THE SUITABILITY OF THIS SOFTWARE                                     
// FOR ANY PARTICULAR PURPOSE.  The software is provided "as is" without                                  
// express or implied warranty of any kind.                                                               
//																											
// These notices must be retained in any copies of any part of this											
// software. 																								
//															    
// Contributed by Pinto, L. L., Cardoso, K. V., Maculan, N. 	
// Designed by Fernandes, K. C. C.								
//															    
// The source code in this file regards the paper, described below, submitted to C&OR.				    
// 															 
// Title: An Exact and Polynomial Approach for a Bi-Objective Integer Programming Problem Regarding Network   
//        Flow Routing.										        
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////Random Topology/////////////////////////////////////////////////////////
#include <map>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <fstream>
#include <iostream>
#include <cstdlib>  
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include <limits.h>
#include <ctime>  

#include <cstring>
   
using namespace std;
   

  typedef struct {
  int source;
  int target;
  int data[100][100];
  map<pair<int, int>, IloBoolVar> variables;
   } Flow;
  int main (int argc, char **argv) {
  ////////////////////////////////////////////////////////////////
  double inicio = clock(); 
  ///////////////////////////////////////////////////////////////  
  // Print FILE 
   ofstream outputFile;
   outputFile.open("lpex1.txt");
   ofstream inputFile;
   inputFile.open("LogFile.txt");
   std::ofstream LogFile("LogFile.txt");
  ////////////////////////////////////////////////////////////////////
  //Initialization
  ///////////////////////////////////////////////////////////////////
  int instancia = 100; //numbers of instances
  int numVertex = 100; //numbers of nodes
  int numFlow = 80;  //number of nodes of network  
   vector<vector<int> > graph(numVertex + 1); 
   Flow *flows = new Flow[numFlow+1];
  int mo = 3; // initial number of nodes and mo is equal also the number of edge added on each new node
  int t = numVertex - mo; //number of nodes that we will put
  int numEdge= 2*mo*t; //edge number of network
  int numvar= numEdge*numFlow;//variables number

/////////////////////////////////////////////////////////////////////////////
//Build a matrix cost of edges
///////////////////////////////////////////////////////////////////////////////
  float**  costEdge = new float*[numVertex];
    for (int i = 0;i < numVertex;i++ ) {
         costEdge[i] = new float[numVertex];
   }
  
 	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build a matrix X where each row is the solution vector of each iteration
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int** X = new int*[numFlow];
   for (int i = 0;i < numFlow;i++ ) {
       X[i] = new int[numvar];
    }
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Build matrix VO that in each row representing the objective value vector
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int** VO = new int*[numFlow];
   for (int i = 0;i < numFlow;i++ ) {
       VO[i] = new int[2];
   }
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Build matrix PO that in each row representing the objective value vector of Pareto-optimal solution 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int** PO = new int*[numFlow];
   for (int i = 0;i < numFlow;i++ ) {
       PO[i] = new int[2];
   }

   
 //////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build jumps matrix numFlow x numFlow that it insertes the number of jumps of each flow in each iteration
////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int** salto = new int*[numFlow];
   for (int i = 0;i < numFlow;i++ ) {
       salto[i] = new int[numFlow];
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int *variable = new int[numvar]; //soluction vector with numvar elements
	int *sumFlow = new int[numvar/numFlow]; //vector of sum of flows in the numEdge edges in each iteration
	int *bottleneck = new int[numFlow]; // vector with the bottleneck in each iteration
	int *ori = new int[numFlow];//vector with the source of each flow
    int *dest = new int[numFlow]; //vector with destination of each flow

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Random Graph
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int *conect = new int[numVertex]; //number of edges insert in each node
	int *nodeAdj = new int[numVertex];//vector that has the number of adjacency nodes for each node
	int sumConectTotal= 0;//total de conecção
	int *position = new int[numVertex+1]; // conection total until node i

//Adjacency matrix of network
     int** inc=new int*[numVertex];// matriz de adjacência
     for(int i = 0; i < numVertex; i++) {
       inc[i]=new int[numVertex];
     }
	
    int** loc = new int*[numVertex];
    for (int i = 0;i < numVertex;i++ ) {
       loc[i] = new int[numVertex];
    }

	float** prob = new float*[numVertex];
   for (int i = 0;i < numVertex;i++ ) {
       prob[i] = new float[numVertex];
   }

	for(int l = 0; l < numVertex; ++l){
        nodeAdj[l] = 0;
     }
	  
    for(int l = 0; l < numVertex; ++l){
        conect[l] = 0;
     }
 
	for (int i=0; i< numVertex; i++){
        for (int j = 0; j < numVertex;j++){
                  inc[i][j] = 0;
	    }
    }
	  
  srand(30);
//connect the node (m_o + 1) to m nodes of existing nodes (here m = m_o)
	for(int i = 1; i <= mo; i++){
		   graph[mo+1].push_back(i);
		   graph[i].push_back(mo+1);
		   inc[mo][i-1] = 1;
		   inc[i-1][mo] = 1;
	       nodeAdj[i-1] = nodeAdj[i-1] + 1;
		   nodeAdj[mo] = nodeAdj[mo];
	       conect[i-1] = conect[i-1] + 1;
		   conect[mo] = conect[mo] + 1;
		   
	}
	int tempo=0;
	position[0] = 0;
	position[1]= conect[0];
	for(int i = 2; i <= mo+1; i++){
		position[i] = position[i-1]+conect[i-1];
	}
	sumConectTotal= position[mo+1];

//calculate the probability and the minimum distance between the connectivitiess
    for(int i = 1; i <= mo+1; i++){
		prob[tempo][i-1]=  float(conect[i-1])/ (sumConectTotal);
	}
	
// insert new nodes

	int vmo; 
	int vpi = 0;
	int novo = 0;
	for(int i = mo+2; i <= numVertex; i++){
		vmo = 0;
		while (vmo < mo){
			novo = (rand()%(sumConectTotal)) + 1;
			for(int j = 1; j <= i-1; ++j){
				if(inc[j-1][i-1] == 0 && inc[i-1][j-1] ==0 & novo <= position[j] && novo > position[j-1]){//j != vp[vmo-2]){
					graph[i].push_back(j);
					graph[j].push_back(i);
					vmo= vmo+1;
					conect[i-1]= conect[i-1]+1;
					conect[j-1]= conect[j-1]+1;
					inc[j-1][i-1] = 1;
		            inc[i-1][j-1] = 1; 
		            nodeAdj[i-1] = nodeAdj[i-1] + 1;
					nodeAdj[j-1] = nodeAdj[j-1] + 1;
					
				}
			}
		}
// updating the positions 
		position[0] = 0;
	    position[1]= conect[0];
	    for(int j = 2; j <= i; j++){
		   position[j] = position[j-1]+conect[j-1];
	    }
		sumConectTotal = position[i];
		// calculate the probbility 
		tempo = tempo+1;
		for(int j = 1; j <= i; j++){
		    prob[tempo][j-1]= float (conect[j-1])/(sumConectTotal);
	    }
		
	}
/*
//Print the graph
	for(int i=1;i <= numVertex; i++){
		for(int j=1;j<= numVertex; j++){
			loc[i-1][j-1]=0;
		}
	}
	
	for(int i=1;i <= numVertex; i++){
		for(int j=1;j<= numVertex; j++){
			if (inc[i-1][j-1] == 1 && loc[i-1][j-1] == 0){
				outputFile<< i << "--" << j << ";" << endl ;
				loc[i-1][j-1] = 1;
				loc[j-1][i-1] = 1;
			}
	     }
     	}
	}

//Calculate maximum number adjacency nodes
  int maxnodeAdj = 0;
  int posmaxEdge = 0;
  for (int l=0; l< numVertex; ++l){
            if (nodeAdj[l] > maxnodeAdj){ 
			   maxnodeAdj = nodeAdj[l];
			   posmaxEdge = l+1;}
   } 
  int posminEdge =0;
  int minnodeAdj = INT_MAX;
  for (int l=0; l< numVertex; ++l){
            if (nodeAdj[l] < minnodeAdj){ 
			   minnodeAdj = nodeAdj[l];
			   posminEdge = l+1;}
   }
   */
/*
  cout <<" Maximum number of  adjacent node = "  << (maxnodeAdj) << endl;
  outputFile <<" Maximum number of  adjacent node ="  << (maxnodeAdj) << endl;
  cout <<" A node with adjacency node maximum   = "  << posmaxEdge << endl;
  outputFile <<" A node with adjacency node maximum = "  << posmaxEdge << endl;
  
  cout <<" Minimum number of  adjacent node = "  << (minnodeAdj) << endl;
  outputFile <<" Minimum number of  adjacent node ="  << (minnodeAdj) << endl;

  cout <<" A node with adjacency node minimum = "  << posminEdge << endl;
  outputFile <<" A node with adjacency node minimum = "  << posminEdge << endl;
    
    */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//The random cost of edges
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int j = 0; j < numVertex; j++) {
      for (int k = j; k < numVertex; k++) {
		  if(inc[j][k]==1 ){
			      costEdge[j][k]= 1;
				  //costEdge[j][k]= 1 + ( rand() % 4 );
				  costEdge[k][j]= costEdge[j][k];
			  }
		  if(inc[j][k]==0){
				  costEdge[j][k]= INT_MAX;
				  costEdge[k][j]= INT_MAX;
			  }
		  
		  }
   }

  /*
  for (int j = 1; j <= numVertex; j++) {
      for (int k = 1;  k <= numVertex; k++) {
		  if(inc[j-1][k-1]==1 ){
	 		  outputFile <<" The edge cost_"<<j<<"_"<<k<< "=" << costEdge[j-1][k-1] << endl;}
	  }
   }
  */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Solving the model while P_epsilon is veasible in each instance ////////////////////////    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int contador;    // Declare the number of instances
   contador=0;                 
      
   while (contador <= instancia)      
  {
      double inicio = clock();
	  contador++;  
///////////////////////////////////////////////////////////////////////////////////////////////////////////
   int epsilon = numFlow;//initial value of epsilon
   int cont=0; //number of times that P_epsilon is feasible
   int numdelete = 0; //number of dominated solutions
   int numPareto = 0; //number of Pareto-optimal solutions
   int prod = 0;
   int pos = 0;
   int cost = 0;
   
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for(int i = 0; i < numFlow; ++i){
	  for(int j = 0; j < numFlow; ++j){
	   salto[i][j] = 0; 
     }
  }            
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Building flows settings  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
//source vector of each flow 
   for(int i = 0; i < numFlow; ++i){
	   ori[i] = 0; 
     }
//destination vector of each flow 
   for(int i = 0; i < numFlow; ++i){
	     dest[i] = 0; 
   }
//Multiple Sources and multiple destinations   
   for(int i = 1;i <= numFlow;i++){
	       memset(flows[i].data, 0, sizeof(flows[i].data));
  }
  
  srand(time(0)); 
  for(int i = 1;i <= numFlow;i++){
	  while( ori[i-1] == dest[i-1] ){    
	       ori[i-1] = (rand()%(numVertex)) + 1; //gera um inteiro de 1 a numVertex
		   dest[i-1] = (rand()%(numVertex)) + 1; }
	   flows[i].source = ori[i-1];
	   flows[i].target = dest[i-1];}
  
//Multiple sources and single destination
  /*
  for(int i = 1;i <= numFlow;i++){
	       memset(flows[i].data, 0, sizeof(flows[i].data));
  }
  int d = (rand()%(numVertex)) + 1;
  for(int i = 1;i <= numFlow;i++){
      ori[i]=d;
      dest[i]=d;
	while(ori[i] == dest[i]){
		ori[i] = (rand()%(numVertex)) + 1; }
		 flows[i].source = ori[i];
	       flows[i].target = dest[i];}
  */
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
 inputFile <<"%%%%%%%%%%%%%%%%%%%%%%%% Grid topology – input flow setting %%%%%%%%%%%%%%%%%%%%%"  << endl;
 inputFile <<"Instance = "<< contador  << endl;
 inputFile <<"Number of nodes = "<< numVertex  << endl;
 inputFile <<"Number of fluxos = "<< numFlow  << endl;
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 for(int i = 1;i <= numFlow;i++){
      inputFile <<"flows["<<i<<"].source = "<< flows[i].source<<";" <<endl;
      inputFile <<"flows["<<i<<"].target = " << flows[i].target<<";"<< endl;
      cout <<"Flow"<< i << "- source ="<< flows[i].source << " - destination =" << flows[i].target << endl;
      //outputFile <<"Fluxo"<< i << "- origem ="<< flows[i].source << " - destino =" << flows[i].target <        
 }
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Building of Loop for resolution P_epsilon
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  do{     
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Building model
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  IloEnv   env;
  IloModel model(env);
  IloBoolVarArray variables(env);
  IloRangeArray constraints(env);
  CPX_ON;//conserva a memória quando possível
  CPX_NODESEL_DFS; //busca em profundidade
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//The building of the variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
 for (int i = 1; i <= numFlow; i++) {
    map<pair<int, int>, bool> added;
    for (int j = 1; j <= numVertex; j++) {
      for (int k = 0; k < graph[j].size(); k++) {
        stringstream varName;
        varName << "x_" << i << "_" << j << "_" << graph[j][k];
        IloBoolVar newVar = IloBoolVar(env, 0, 1, varName.str().c_str());
        flows[i].variables[make_pair(j, graph[j][k])] = newVar;
        variables.add(newVar);
		
		}
    }
  }
  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
//Building the constraint 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
 for (int i = 1; i <= numFlow; i++) {

//Constraint 1: 
	stringstream outConstraintName;
    stringstream inConstraintName;
    outConstraintName << "Out_flow_" << i;
    inConstraintName << "In_flow_" << i;

    IloExpr outConstraint(env);
    IloExpr inConstraint(env);

     for (int j = 0; j < graph[flows[i].source].size(); j++) {
          outConstraint += flows[i].variables[make_pair(flows[i].source, graph[flows[i].source][j])];
          outConstraint -=   flows[i].variables[make_pair(graph[flows[i].source][j], flows[i].source)];

    }

    for (int j = 0; j < graph[flows[i].target].size(); j++) {
         inConstraint += flows[i].variables[make_pair(flows[i].target, graph[flows[i].target][j])];       
         inConstraint -= flows[i].variables[make_pair(graph[flows[i].target][j], flows[i].target)];

    }
//Constraint 2: 
// this loop below excludes the source and destination vertices
     for (int j = 1; j <= numVertex; j++) {
       if (j == flows[i].source || j == flows[i].target)
        continue;

	  
      IloExpr pathConstraint(env);
      stringstream pathConstraintName;

      pathConstraintName << "Path_" << i << "_" << j;

      for (int k = 0; k < graph[j].size(); k++) {
        pathConstraint += flows[i].variables[make_pair(graph[j][k], j)];
        pathConstraint -= flows[i].variables[make_pair(j, graph[j][k])];
      }

      constraints.add(IloRange(env, 0, pathConstraint, 0, pathConstraintName.str().c_str()));
       }

      constraints.add(IloRange(env, 1, outConstraint, 1, outConstraintName.str().c_str()));
      constraints.add(IloRange(env, -1, inConstraint, -1, inConstraintName.str().c_str()));
       }
  
 //Constraint 3:  
  for (int i = 1; i <= numVertex; i++) {
    for (int j = 0; j < graph[i].size(); j++) {//aquizero <
      
	  IloExpr alphaConstraint(env);
      stringstream alphaConstraintName;

      alphaConstraintName << "Epsilon_" << i << "_" << graph[i][j];
	  

      for (int k = 1; k <= numFlow; k++) {
        alphaConstraint += flows[k].variables[make_pair(i, graph[i][j])];
		
      }


      constraints.add(IloRange(env, 0, alphaConstraint, epsilon, alphaConstraintName.str().c_str()));

    }
  }
   
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //objective function of P_epsilon
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   IloExpr obj(env);
   stringstream objName;
   objName << "Objective_" << obj << endl;
   for (int i = 1; i <= numVertex; i++) {
	   for (int j = 0; j < graph[i].size(); j++) {
   			for (int k = 1; k <= numFlow; k++) {
					obj += costEdge[i-1][graph[i][j]-1]*flows[k].variables[make_pair(i, graph[i][j])];
				 	
	                   		 
			}
			 
	}
   }
 
  
//Building the model P_epsilon  
   model.add(variables);
   model.add(IloMinimize(env, obj));
   model.add(constraints);
  //Print input 
   if(cont == 0){
     cout<<"Number of nodes = "<< numVertex << endl;
     cout<<"Number of edges = "<< numEdge << "=" << numEdge << endl;
     cout<<"Number of flows = "<< numFlow << endl;
     cout<<"Number of variables = "<< numvar << endl;
    
     //outputFile <<"Number of nodes = "<< numVertex << endl;
     //outputFile <<"Number of edges = "<< numEdge << endl;
     //outputFile <<"Number of flows = "<< numFlow << endl;
   //outputFile <<"Number of variables = "<< numvar << endl;
    }
   //Print the modelo
   /*if(cont == 0){
       cout << variables << endl;
       outputFile << variables << endl;
       cout << model << endl;
       outputFile << model << endl;}
    */
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Optimize the problem P_epsilon
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      IloCplex cplex(model);
	  //cplex.setOut(LogFile); 
   
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if ( cplex.solve() ) {
	   cont = cont + 1;}
   else {
	  env.error() << "Failed to optimize LP" << endl;
	  cout << "The P-Alpha is infeasible to iterations >=  " << cont+1 << endl;  
	  //outputFile << "The P-Alpha is infeasible to iterations >=  " << cont+1 << endl;
	  cout <<"The quantity of times that P_epsilon was solved = " << cont+1 << endl;
      cout << "The quantity of solutions dominated = " << numdelete << endl;
      cout << "The quantity of optimal-Pareto solutions, i.e., |X| = " << numPareto << endl;
      //outputFile <<"The quantity of times that P_epsilon was solved = " << cont+1 << endl;
      //outputFile << "The quantity of solutions dominated = " << numdelete << endl;
      //outputFile << "The quantity of optimal-Pareto solutions, i.e., |PO| = " << numPareto << endl;
      double fim = clock();
      double  tempo = fim - inicio;
	  //outputFile << "The time of execution until here is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;	
	  cout << "The time of execution until here is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;
	  ////////cout << numPareto << " " << numdelete << " " << numPareto+ numdelete << " " << cont+1  << " " << tempo / (double)CLOCKS_PER_SEC << " " << VO[0][0] << " " << VO[cont-1][0] <<  endl;
      ///////outputFile << numPareto << " " << numdelete << " " << numPareto+ numdelete << " " << cont+1  << " " << tempo / (double)CLOCKS_PER_SEC << " " << VO[0][0] << " " << VO[cont-1][0] <<  endl;
	  int primeiro = 0;
	  for(int i=0;i<=cont;i++){
		  if( PO[i][0] != INT_MAX && primeiro == 0 ){
			  cout << numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont+1  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " <<PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
		      outputFile << numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont+1  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " <<PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
	          primeiro = primeiro + 1;
			  /*for (int k = 0; k < numFlow; k++ ){
				  outputFile <<salto[0][k] << "   " << salto[cont-1][k] << endl;}*/
		  }
	  }
	  
//Clean 
	  cplex.clearModel();
	  cplex.clear();
	  variables.end();
      constraints.endElements();
      constraints.end();
      cplex.end();
      model.end();
	  env.end();
	  break;
   }
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
//Print the model results
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   env.out() << "Solution status = " << cplex.getStatus() << endl;
   env.out() << "Solution value  = " << cplex.getObjValue() << endl;
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
//Updating the epsilon value
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////      
   
   for(int i = 0; i < numvar; ++i){
	    variable[i] = cplex.getValue(variables[i]);
   }
   
   for(int i = 0; i < numEdge; ++i){
	    sumFlow[i] = 0;
   }

   for(int j = 0; j < numEdge; ++j){
     for(int l = 0; l < numFlow; ++l){
		 prod = l*numEdge;
		 pos = j + (l*numEdge);
		 sumFlow[j] += variable[pos];
     }
   }
//vector with the bottleneck in each iteration after to solve P_epsilon
    for (int i = 0; i < numFlow; ++i){
        bottleneck[i] = 0;
    }
 
    for (int l=0; l< numEdge; ++l){
             if (sumFlow[l] > bottleneck[cont-1]) 
			   bottleneck[cont-1] = sumFlow[l];
    } 

    cout <<"The bottleneck in the iteration "  << cont << " = " << bottleneck[cont-1] << endl;
    //outputFile <<"The bottleneck in the iteration "  << cont << " = " << bottleneck[cont-1] << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    epsilon = bottleneck[cont-1] - 1;  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Store the total of jumps of each flow from the current iteration 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 for(int k = 0; k < numFlow; ++k){
		 for(int i = 0; i < numEdge; ++i){
		 salto[cont-1][k] = salto[cont-1][k] +variable[i + k*numEdge]; 
     }
  }

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////   //Store the solutions (variables) of current iteration  in the row (cont-1) of matrix X with (numFLow)columns
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
   for(int j = 0; j < numvar; ++j){
        X[cont-1][j] = variable[j]; 
   }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
//Store the solutions (variables) of current iteration  in the row (cont-1) of matrix X with (numFLow)columns
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    VO[cont-1][0] = bottleneck[cont-1];
    VO[cont-1][1] = cplex.getObjValue();
	VO[cont][0] = INT_MAX;
    VO[cont][1] = INT_MAX;

	PO[cont-1][0] = VO[cont-1][0];
    PO[cont-1][1] = VO[cont-1][1];
	PO[cont][0] = INT_MAX;
    PO[cont][1] = INT_MAX;
    
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Cleanup
///////////////////////////////////////////////////////////////////////////////////////////////////////
	cplex.clearModel();
	cplex.clear();
	cplex.clearModel();
	variables.end();
    constraints.endElements();
    constraints.end();
    cplex.end();
    model.end();
	env.end();
	
	cout << "The objective vector [f1(x)  f2(x)] in the iteration " << cont << " =  [" << VO[cont-1][0] << "  " << VO[cont-1][1]  << "]" <<  endl; 
    //outputFile << "The objective vector [f1(x)  f2(x)] in the iteration " << cont << " =  [" << VO[cont-1][0] << "  " << VO[cont-1][1]  << "]" <<  endl; 
    cout << "The value of epsilon in the next iteration , ie , iteration " << cont+1 << " = " << epsilon << endl;
    //outputFile << "The value of alpha in the next iteration , ie , iteration  " << cont+1 << "= " << alpha << endl;
  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
//Deleting the dominated solutions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
// when solution of iteration i dominates the solution of iteration i-1 the line i-1 of X is replaced by another one
// with null elements .  
	for ( int i = 0; i < numvar; i++ ){
         if( cont > 1 && VO[cont-1][1] == VO[cont-2][1])
             X[cont-2][i] = 0;
   }
 
  // when solution of iteration i dominates the solution of iteration i-1 the line i-1 of VO is replaced by another one
  // with infinite elements.
   for ( int i = 0; i < 2; i++ ){
	   if( cont > 1 && VO[cont-1][1] == VO[cont-2][1]){
             PO[cont-2][i] = INT_MAX;
			 
	   }
   }

//calculate the dominated solution number
   if( cont > 1 && VO[cont-1][1] == VO[cont-2][1]){
	   numdelete = numdelete+1;
    }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << " The Matrix VO presents in each row the objective vector of generated solutions, i.e, " << endl;
   //outputFile << " The Matrix PO presents in each row the objective vector of generated solutions, i.e, " << endl;
   
   for (int i = 0; i <= cont; i++ ){
	 cout <<  "VO_" << i << " * " << VO[i][0] << "  " << VO[i][1] << endl;
     //outputFile <<  "VO_" << i << " * " << VO[i][0] << "  " << VO[i][1]  << endl;
    }
	   
    cout << " The Matrix PO presents in each row the objective vector of Pareto-opitmal solutions, i.e, " << endl;
     //outputFile << " The Matrix PO presents in each row the objective vector of Pareto-opitmal solutions, i.e, " << endl;
   
   for (int i = 0; i <= cont; i++ ){
	 cout <<  "PO_" << i << " * " << PO[i][0] << "  " << PO[i][1] << endl;
     //outputFile <<  "PO_" << i << " * " << PO[i][0] << "  " << PO[i][1]  << endl;
   }
 
//Calculate the cardinality of complete minimal  set of Pareto-optimal solutions 
    numPareto = cont - numdelete; 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (epsilon == 0 ){
		cout << "The P-epsilon reached the last iteration = " << cont << endl;///atingiu o valor minimo de alpha em cont
		//outputFile << "The P-epsilon reached the last iteration = " << cont << endl;
        cout <<"The quantity of times that P_epsilon was solved = " << cont << endl;
        cout << "The quantity of solutions dominated = " << numdelete << endl;
        cout << "The quantity of optimal-Pareto solutions, i.e., |X| = " << numPareto << endl;
        //outputFile <<"The quantify of times that P_alpha was solved = " << cont << endl;
        //outputFile << "The quantity of solutions dominated = " << numdelete << endl;
        //outputFile << "The quantity of optimal-Pareto solutions, i.e., |X| = " << numPareto << endl;
        
        double fim = clock();
        double  tempo = fim - inicio;
	    ////outputFile << "The time of execution is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;	
	    //cout << "The time of execution is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;
		//cout << numPareto << " " << numdelete << " " << numPareto+ numdelete << " " << cont  << " " << tempo / (double)CLOCKS_PER_SEC << " " << VO[0][0] << " " << VO[cont-1][0] <<  endl;
        //outputFile << numPareto << " " << numdelete << " " << numPareto+ numdelete << " " << cont  << " " << tempo / (double)CLOCKS_PER_SEC << " " << PO[0][0] << " " << PO[cont-1][0] <<  endl;  
        int primeiro = 0;
	    for(int i=0;i<=cont;i++){
		  if( PO[i][0] != INT_MAX && primeiro == 0 ){
			  cout << numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " <<PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
		      outputFile<< numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " <<PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
			  primeiro = primeiro + 1;
			  /*for (int k = 0; k < numFlow; k++ ){
				  outputFile <<salto[0][k] << "   " << salto[cont-1][k] << endl;}*/
		  }
	  }
   }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
     double fim = clock();
     double  tempo = fim - inicio;
	 //outputFile << "The running time until here is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;	
	 cout << "The running time  until here is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;
	 
	 
	
	} while (epsilon != 0);
  

  } 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   system ("pause");
   outputFile.close();
       
     delete [] nodeAdj;
     delete [] ori;
     delete [] dest;
     delete [] flows;
	 delete [] variable;
     delete [] sumFlow;
     delete [] bottleneck;
	 delete [] conect;
	 delete [] position;
	
	 for(int i=0 ; i < numFlow;i++) {
      delete [] X[i];}
      delete [] X;
 
    for(int i=0 ; i < numFlow;i++) {
       delete [] PO[i];}
       delete [] PO;
   
	 
	for(int i=0 ; i < numFlow;i++) {
     delete [] VO[i];}
     delete [] VO;
	 
	 
	 for(int i=0 ; i < numVertex;i++) {
      delete [] inc[i];}
      delete [] inc;

     for(int i=0 ; i < numVertex;i++) {
      delete [] loc[i];}
      delete [] loc;
	 
	 for(int i=0 ; i < numVertex;i++) {
      delete [] prob[i];}
      delete [] prob;

    for(int i=0 ; i < numFlow;i++) {
     delete [] salto[i];}
     delete [] salto;

	  for(int i=0 ; i < numVertex;i++) {
      delete [] costEdge[i];}
      delete [] costEdge;


   return 0;
} 