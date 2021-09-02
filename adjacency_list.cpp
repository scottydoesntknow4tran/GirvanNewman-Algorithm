//----------------------------------------------------------------------
// NAME: S. Tornquist
// FILE: adjacency_list.cpp
// DATE: Spring 2021
// DESC: fucntion implementations of adjacency_list.h
//----------------------------------------------------------------------


#ifndef ADJACENCY_LIST_CPP
#define ADJACENCY_LIST_CPP

#include "graph.h"
#include "adjacency_list.h"

  // default constructor
AdjacencyList :: AdjacencyList(int vertex_count){
     number_of_edges = 0;
     number_of_vertices = vertex_count;
     adj_list_in  = new Node*[number_of_vertices]; //creating tables
     adj_list_out  = new Node*[number_of_vertices];
     for(int i=0; i<number_of_vertices; i++){ // sets every pointer to null ptr to indicate empty
          adj_list_out[i] = nullptr;
          adj_list_in[i] = nullptr;
     }
  };

  // destructor
AdjacencyList :: ~AdjacencyList(){
     make_empty(); //clears the lists
};

  // copy constructor
AdjacencyList :: AdjacencyList(const AdjacencyList& rhs){
     number_of_vertices = -1; //used to indicate the lhs being passed to the assignment operator is new and does not need to be cleared
     *this = rhs; // defer to assignment operator
};

  // assignment operator
AdjacencyList& AdjacencyList :: operator=(const AdjacencyList& rhs){
     if(this != &rhs){ // list1= list1
          if(number_of_vertices != -1){ //had to make sure the garbage wasnt triggering make empty when object was uninitialized
               make_empty();
          }
          number_of_vertices = rhs.number_of_vertices;
          number_of_edges = 0;
          adj_list_in  = new Node*[number_of_vertices];
          adj_list_out  = new Node*[number_of_vertices];
          for(int i=0; i<number_of_vertices; i++){ // sets every pointer to null ptr to indicate empty
               adj_list_out[i] = nullptr;
               adj_list_in[i] = nullptr;
          }
          for(int i=0; i < number_of_vertices; i++){ // moving old values in reverse order of the chains as to keep them in their orginal order
               if(rhs.adj_list_out[i] != nullptr){
                    Node* ptr = rhs.adj_list_out[i];
                    while(ptr != nullptr){// goes to last value in chain
                         set_edge(i,ptr->edge, ptr->vertex);
                         ptr = ptr->next;
                    }
               }
          }    
     }
     //return lhs(this)
     return *this;
};

  // add or update an edge between existing vertices in the graph
void AdjacencyList :: set_edge(int v1, double edge_label, int v2){
     if((0 <= v1) && (v1 < number_of_vertices) && (0 <= v2) && (v2 < number_of_vertices)){ // making sure the vertices exist
          //checking for existence of edge already
          Node* ptr = adj_list_out[v1]; 
          bool found = false;
          while(ptr != nullptr){
               if(ptr->vertex == v2){
                    //found, check for update
                    //if(ptr->edge != edge_label){
                         ptr->edge = edge_label; //updating
                         found = true;
                         break;
                    //}
               }
               ptr = ptr->next;
          }

          //if edge already exists, update edge value
          if(found == true){ // updating other list
               ptr = adj_list_in[v2]; 
               while(ptr != nullptr){
                    if(ptr->vertex == v1){
                         //found, check for update
                         //if(ptr->edge != edge_label){
                              ptr->edge = edge_label; //updating
                              break;
                         //}
                    }
               ptr = ptr->next;
               }
          }
          
          //if the node was not found, create it
          if(found == false){ // if the edge doesnt already exist
               //adding edge to out list at v1 location
               ptr = adj_list_out[v1]; 

               // creating new node and filling data
               Node* itr = new Node; 
               itr->edge = edge_label;
               itr->vertex = v2;

               //setting new node's next to whatever the head was pointing to and moving head to new node(inserting at the front)
               itr->next = adj_list_out[v1];
               adj_list_out[v1] = itr;


               //Repeating for in list
               //adding edge to out list at v1 location
               ptr = adj_list_in[v1]; 

               // creating new node and filling data
               itr = new Node; 
               itr->edge = edge_label;
               itr->vertex = v1;

               //setting new node's next to whatever the head was pointing to and moving head to new node(inserting at the front)
               itr->next = adj_list_in[v2];
               adj_list_in[v2] = itr;

               number_of_edges++; //updating number of edges
          }
     }

};

  // check if the given edge is in the graph
bool AdjacencyList :: has_edge(int v1, int v2) const{
     if((0 <= v1) && (v1 < number_of_vertices) && (0 <= v2) && (v2 < number_of_vertices)){ // making sure the vertices exist
          Node* ptr = adj_list_out[v1]; //go to list of outcoming vertices from v1
          while(ptr !=  nullptr){ // loopr through linked list, adding to return list
               if(ptr->vertex == v2){
                    return true; //if connection is found
               }
               ptr = ptr->next;
          }
     }
     return false; //if it is not found
};  

  // returns true if there is an edge between v1 and v2 in the graph
  // the output parameter returns the edge label
bool AdjacencyList :: get_edge(int v1, int v2, double& edge) const{
     if((0 <= v1) && (v1 < number_of_vertices) && (0 <= v2) && (v2 < number_of_vertices)){ // making sure the vertices exist
          Node* ptr = adj_list_out[v1]; //go to list of outcoming vertices from v1
          while(ptr !=  nullptr){ // loopr through linked list, adding to return list
               if(ptr->vertex == v2){
                    edge = ptr->edge;
                    return true; //if connection is found
               }
               ptr = ptr->next;
          }
     }
     edge = 0;
     return false; //if it is not found
}; 
  
  // get all vertices on an outgoing edge from given vertex
void AdjacencyList :: connected_to(int v1, std::list<int>& vertices) const{
     if((0 <= v1) && (v1 < number_of_vertices)){ // making sure the vertices exist
          Node* ptr = adj_list_out[v1]; //go to list of outcoming vertices
          while(ptr !=  nullptr){ // loop through linked list, adding to return list
               vertices.push_front(ptr->vertex);
               ptr = ptr->next;
          }
     }
};

  // get all vertices on an incoming edge to given vertex
void AdjacencyList :: connected_from(int v2, std::list<int>& vertices) const{
     if((0 <= v2) && (v2 < number_of_vertices)){ // making sure the vertices exist
          Node* ptr = adj_list_in[v2]; //go to list of incoming vertices
          if(ptr != nullptr){
               while(ptr !=  nullptr){ // loop through linked list, adding to return list
                    vertices.push_front(ptr->vertex);
                    ptr = ptr->next;
               }
          }
     }
};

  // get all vertices adjacent to a vertex, that is, connected to or connected
  // from the vertex (may return duplicate vertices)
void AdjacencyList :: adjacent(int v, std::list<int>& vertices) const{
     //outgoing edges
     connected_to(v, vertices);

     //incoming edges
     connected_from(v, vertices);
     // this code does not prevent duplicates
};

  // get number of nodes in the graph
int AdjacencyList :: vertex_count() const{
     return number_of_vertices;
}; 

  // get number of edges in the graph
int AdjacencyList :: edge_count() const{
     return number_of_edges;
}; 

  
  // helper function to delete adj_list
void AdjacencyList :: make_empty(){
     if(adj_list_out != nullptr){
          for(int i=0; i<number_of_vertices; i++){// deletes all old nodes 
               if(adj_list_out[i]!=nullptr){
                    Node* ptr = adj_list_out[i]; //set to head of list
                    Node* itr = ptr->next;
                    while(itr!=nullptr){ //iterate down list deleting all nodes
                         delete ptr;
                         ptr = itr;
                         itr = itr->next;
                    }
               delete ptr;
               }
          }
          delete[] adj_list_out;//deleting array
     }
     // repeating for list in
     if(adj_list_in != nullptr){
          for(int i=0; i<number_of_vertices; i++){// deletes all old nodes 
               if(adj_list_in[i]!=nullptr){
                    Node* ptr = adj_list_in[i]; //set to head of list
                    Node* itr = ptr->next;
                    while(itr!=nullptr){ //iterate down list deleting all nodes
                         delete ptr;
                         ptr = itr;
                         itr = itr->next;
                    }
               delete ptr;
               }
          }
          delete[] adj_list_in; //deleting array
          number_of_edges =0;
          number_of_vertices=0;
     }
};

//hw5
//removing all edges btween src and end
bool AdjacencyList ::remove_edge(int src,int end){
     if((0 <= src) && (src < number_of_vertices) && (0 <= end) && (end < number_of_vertices)){ // making sure the vertices exist
          bool found = false;

          //checking first link
          //checking out list from src
          Node* ptr = adj_list_out[src]; //go to list of outcoming vertices from src
          if(ptr != nullptr){
               if(ptr->vertex == end){ //if the next vertex is the end vertice
                    adj_list_out[src] = ptr->next; //moving link
                    delete ptr;//deleting itr
                    found =true;
               }
          }
          //checking in list from end
          ptr = adj_list_in[end]; //go to list of outcoming vertices from v1
          if(ptr != nullptr){
               if(ptr->vertex == src){ //if the next vertex is the end vertice
                    adj_list_in[end] = ptr->next; //moving link
                    delete ptr;//deleting ptr 
                    found =true;
               }
          }

          //checking other links
          //checking out list from src
          ptr = adj_list_out[src]; //go to list of outcoming vertices from v1
          if(ptr!= nullptr){
               Node* itr = ptr->next;
               while(itr !=  nullptr){ // looping through linked list, adding to return list
                         if(itr->vertex == end){ //if the next vertex is the end vertice
                              ptr->next = itr->next; //setting ptr to itrs next
                              delete itr;//deleting itr
                              found ==true;
                              break;
                         }
                    ptr = ptr->next;
                    itr = itr->next;
               }
          }
          //checking in list from end
          ptr = adj_list_in[end]; //go to list of outcoming vertices from v1
          if(ptr!= nullptr){
               Node* itr = ptr->next;
               while(itr !=  nullptr){ // looping through linked list, adding to return list
                         if(itr->vertex == src){ //if the next vertex is the end vertice
                              ptr->next = itr->next; //setting ptr to itrs next
                              delete itr;//deleting itr
                              found = true;
                              break;
                         }
                    ptr = ptr->next;
                    itr = itr->next;
               }
          }
          if(found){
               number_of_edges--;
               return true;
          }
          return false;
     }
     return false;
};


#endif

