//----------------------------------------------------------------------
// NAME: S. Tornquist
// FILE: adjacency_matrix.cpp
// DATE: Spring 2021
// DESC: function implementations of adjacency_matrix.h
//----------------------------------------------------------------------


#ifndef ADJACENCY_MATRIX_CPP
#define ADJACENCY_MATRIX_CPP

#include "graph.h"
#include "adjacency_matrix.h"

  // default constructor
  AdjacencyMatrix :: AdjacencyMatrix(int vertex_count){
     number_of_edges = 0;
     number_of_vertices = vertex_count;
     adj_matrix = new int*[number_of_vertices*number_of_vertices]; //creating new matrix of int pointers
     for(int i=0; i<(number_of_vertices*number_of_vertices); i++){ // sets every pointer to null ptr to indicate empty
          adj_matrix[i] = nullptr;
     }
  };

  // destructor
  AdjacencyMatrix :: ~AdjacencyMatrix(){
       make_empty(); //clearing matrix
  };

  // copy constructor
  AdjacencyMatrix :: AdjacencyMatrix(const AdjacencyMatrix& rhs){
     number_of_vertices = -1; //used to indicate the lhs being passed to the assignment operator is new and does not need to be cleared
     *this = rhs; // defer to assignment operator
  };

  // assignment operator
  AdjacencyMatrix& AdjacencyMatrix :: operator=(const AdjacencyMatrix& rhs){
       if(this != &rhs){ // list1= list1
          if(number_of_vertices != -1){ //had to make sure the garbage wasnt triggering make empty when object was uninitialized
               make_empty(); //clearing lhs if it needs to be cleared
          }
          number_of_vertices = rhs.number_of_vertices;
          number_of_edges = rhs.number_of_edges;
          adj_matrix = new int*[number_of_vertices*number_of_vertices];
          for(int i=0; i<(number_of_vertices*number_of_vertices); i++){ // sets every pointer to null ptr to indicate empty
               adj_matrix[i] = nullptr;
          }
          for(int i=0; i < number_of_vertices*number_of_vertices; i++){ // moving old values in reverse order of the chains as to keep them in their orginal order
               if(rhs.adj_matrix[i] != nullptr){
                    adj_matrix[i] = new int;
                    *adj_matrix[i] = *rhs.adj_matrix[i]; 
               }
          }    
     }
     //return lhs(this)
     return *this;
  };

  // add or update an edge between existing vertices in the graph
  void AdjacencyMatrix :: set_edge(int v1, int edge_label, int v2){
     if((0 <= v1) && (v1 < number_of_vertices) && (0 <= v2) && (v2 < number_of_vertices)){ // making sure it is valid vertices
          if(adj_matrix[index(v1,v2)] == nullptr){ //new edge
               adj_matrix[index(v1,v2)] = new int; //creaitn new int to point to
               *adj_matrix[index(v1,v2)] = edge_label; //adding value
               number_of_edges++;
          }
          else{ //updating old edge
               *adj_matrix[index(v1,v2)] = edge_label;
          }
          
     }
  };

  // check if the given edge is in the graph
  bool AdjacencyMatrix :: has_edge(int v1, int v2) const{
     if((0 <= v1) && (v1 < number_of_vertices) && (0 <= v2) && (v2 < number_of_vertices)){ // making sure it is valid vertices
          if(adj_matrix[index(v1,v2)] != nullptr){
               return true;
          }
     }
     return false; //otherwise not found
  };  

  // returns true if there is an edge between v1 and v2 in the graph
  // the output parameter returns the edge label
  bool AdjacencyMatrix :: get_edge(int v1, int v2, int& edge) const{
       if((0 <= v1) && (v1 < number_of_vertices) && (0 <= v2) && (v2 < number_of_vertices)){// making sure it is valid vertices
          if(adj_matrix[index(v1,v2)] != nullptr){
               edge = *adj_matrix[index(v1,v2)]; //same at has_edge but also return the value of the edges to the return variable
               return true;
          }
     }
     return false;
  }; 
  
  // get all vertices on an outgoing edge from given vertex
  void AdjacencyMatrix :: connected_to(int v1, std::list<int>& vertices) const{
       if((0 <= v1) && (v1 < number_of_vertices)){ // making sure the vertices exist go to list of outcoming vertices
          for(int i =0; i < number_of_vertices; i++){
               if(adj_matrix[index(v1,i)] != nullptr){ //iterating down the row adding all non nullptr's indicating an edge
                    vertices.push_back(i); //adding
               }
          }
     }
  };

  // get all vertices on an incoming edge to given vertex
  void AdjacencyMatrix :: connected_from(int v2, std::list<int>& vertices) const{
       if((0 <= v2) && (v2 < number_of_vertices)){ // making sure the vertices exist go to list of outcoming vertices
          for(int i =0; i < number_of_vertices; i++){
               if(adj_matrix[index(i,v2)] != nullptr){//iterating down the column adding all non nullptr's indicating an edge
                    vertices.push_front(i);
               }
          }
     }
  };

  // get all vertices adjacent to a vertex, that is, connected to or connected
  // from the vertex (may return duplicate vertices)
  void AdjacencyMatrix :: adjacent(int v, std::list<int>& vertices) const{
       connected_from(v,vertices); //iterating down the column adding all non nullptr's indicating an edge

       connected_to(v,vertices); //iterating down the row adding all non nullptr's indicating an edge
       // this method does not prevent overlaps
  };

  // get number of nodes in the graph
  int AdjacencyMatrix :: vertex_count() const{
       return number_of_vertices;
  }; 

  // get number of edges in the graph
  int AdjacencyMatrix :: edge_count() const{
       return number_of_edges;
  }; 

  
  // helper function to delete adj_list
  void AdjacencyMatrix :: make_empty(){
       for(int i =0; i < number_of_vertices*number_of_vertices; i++){
               if(adj_matrix[i] != nullptr){
                    delete adj_matrix[i];//deleting all ptr that that are allocated
               }
          }
          delete[] adj_matrix; // deleting array of ptrs
          number_of_vertices =0;
          number_of_edges =0;
  };

  int AdjacencyMatrix :: index(int x, int y) const{
     int ind = (x*number_of_vertices) + y; //calculating index for 1D array given row and column indexes of VxV matrix
     return ind;
  };


//hw5
//removing all edges btween src and end
bool AdjacencyMatrix ::remove_edge(int src,int end){
     if((0 <= src) && (src < number_of_vertices) && (0 <= end) && (end < number_of_vertices)){// making sure it is valid vertices
          if(adj_matrix[index(src,end)] != nullptr){
               delete adj_matrix[index(src,end)]; //same at has_edge but also return the value of the edges to the return variable
               adj_matrix[index(src,end)] = nullptr;
               number_of_edges--;
               return true;
          }
     }
     return false;
};


#endif
