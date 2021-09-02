//----------------------------------------------------------------------
// Name: 
// File: project_tests.cpp
// Date: Spring 2021
// Desc: Unit tests for list and matrix graph implementations
//----------------------------------------------------------------------


#include <iostream>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <gtest/gtest.h>
#include "graph.h"
#include "adjacency_list.h"
#include "adjacency_matrix.h"


using namespace std;


//----------------------------------------------------------------------
// Helper functions for testing
//----------------------------------------------------------------------

void print_graph(const Graph& g)
{
  int n = g.vertex_count();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (g.has_edge(i,j)) {
        double edge_label;
        g.get_edge(i, j, edge_label);
        cout << "(" << i << "," << edge_label << "," << j << ")" << endl;
      }
    }
  }
}


// Project tests

//----------------------------------------------------------------------
// Adjacency List Tests
//----------------------------------------------------------------------

//testing helper algorithm for number of shorstest path from a source vertex

TEST(AdjacencyListTest, BasicNumberOfShortestPaths) {
  AdjacencyList g(6);
  g.set_edge(0, 0, 1);
  g.set_edge(0, 0, 3);
  g.set_edge(1, 0, 2);
  g.set_edge(1, 0, 4);
  g.set_edge(2, 0, 5);
  g.set_edge(3, 0, 4);
  g.set_edge(4, 0, 5);
  ASSERT_EQ(7, g.edge_count());
  ASSERT_EQ(6, g.vertex_count());
  Map num_shortest_paths;
  Map distance;
  std::list<std::pair<int,std::pair<int,int>>> ordered_by_distance;
  g.number_of_shortest_paths(0, num_shortest_paths, distance, ordered_by_distance);
  // check size
  ASSERT_EQ(g.vertex_count(), num_shortest_paths.size());
  ASSERT_EQ(g.vertex_count(), distance.size());
  ASSERT_EQ(g.vertex_count(), ordered_by_distance.size());
  //checking number of shortest paths
  ASSERT_EQ(1, num_shortest_paths[0]);
  ASSERT_EQ(1, num_shortest_paths[1]);
  ASSERT_EQ(1, num_shortest_paths[2]);
  ASSERT_EQ(1, num_shortest_paths[3]);
  ASSERT_EQ(2, num_shortest_paths[4]);
  ASSERT_EQ(3, num_shortest_paths[5]);
  //checking distance
  ASSERT_EQ(0, distance[0]);
  ASSERT_EQ(1, distance[1]);
  ASSERT_EQ(2, distance[2]);
  ASSERT_EQ(1, distance[3]);
  ASSERT_EQ(2, distance[4]);
  ASSERT_EQ(3, distance[5]);
  //checking ordering by distance

  //ASSERT_EQ(ordered_by_distance.front(), distance[5]);
}


TEST(AdjacencyListTest, BasicEdgeScoreTest) {
  AdjacencyList g(6);
  g.set_edge(0, 0, 1);
  g.set_edge(0, 0, 3);
  g.set_edge(1, 0, 2);
  g.set_edge(1, 0, 4);
  g.set_edge(2, 0, 5);
  g.set_edge(3, 0, 4);
  g.set_edge(4, 0, 5);
  ASSERT_EQ(7, g.edge_count());
  ASSERT_EQ(6, g.vertex_count());
  AdjacencyList g1(6);
  g.girvan_newman(g1);
  // check size
  ASSERT_EQ(g1.edge_count(), 3);
  //checking the edges
  ASSERT_EQ(g1.has_edge(0,3) or g1.has_edge(3,0), true);
  ASSERT_EQ(g1.has_edge(1,4) or g1.has_edge(4,1), true);
  ASSERT_EQ(g1.has_edge(2,5) or g1.has_edge(5,2), true);
  
  //checking communities
  Map communities;
  g1.bfs_connected_components(communities);
  int community1 = communities[0];
  int community2 = communities[1];
  int community3 = communities[2];
  ASSERT_NE(community1, community2);
  ASSERT_NE(community1, community3);
  ASSERT_NE(community1, community3);

  ASSERT_EQ(community1, communities[3]);

  ASSERT_EQ(community2, communities[4]);

  ASSERT_EQ(community3, communities[5]);
}

TEST(AdjacencyListTest, MediumGrivanNewmanTest) {
  AdjacencyList g(7);
  g.set_edge(0, 0, 1);
  g.set_edge(1, 0, 2);
  g.set_edge(2, 0, 0);
  g.set_edge(1, 0, 3);
  g.set_edge(3, 0, 5);
  g.set_edge(3, 0, 4);
  g.set_edge(4, 0, 5);
  g.set_edge(5, 0, 6);
  g.set_edge(6, 0, 3);
  ASSERT_EQ(9, g.edge_count());
  ASSERT_EQ(7, g.vertex_count());
  AdjacencyList g1(7);
  g.girvan_newman(g1);
  // check size
  ASSERT_EQ(g1.edge_count(), 8);
  //checking the edge values
  ASSERT_EQ(g1.has_edge(0,1) or g1.has_edge(1,0), true);
  ASSERT_EQ(g1.has_edge(1,2) or g1.has_edge(2,1), true);
  ASSERT_EQ(g1.has_edge(2,0) or g1.has_edge(0,2), true);
  ASSERT_EQ(g1.has_edge(4,3) or g1.has_edge(3,4), true);
  ASSERT_EQ(g1.has_edge(5,4) or g1.has_edge(4,5), true);
  ASSERT_EQ(g1.has_edge(3,5) or g1.has_edge(5,3), true);
  ASSERT_EQ(g1.has_edge(6,5) or g1.has_edge(5,6), true);
  ASSERT_EQ(g1.has_edge(3,6) or g1.has_edge(6,3), true);
  //checking communities
  Map communities;
  g1.bfs_connected_components(communities);
  int community1 = communities[0];
  int community2 = communities[3];
  ASSERT_NE(community1, community2);
  ASSERT_EQ(community1, communities[1]);
  ASSERT_EQ(community1, communities[2]);

  ASSERT_EQ(community2, communities[4]);
  ASSERT_EQ(community2, communities[5]);
  ASSERT_EQ(community2, communities[6]);
  
}

TEST(AdjacencyListTest, LargeGrivanNewmanTest) {
  AdjacencyList g(10);
  g.set_edge(0, 0, 1);
  g.set_edge(1, 0, 2);
  g.set_edge(2, 0, 0);
  g.set_edge(1, 0, 3);
  g.set_edge(3, 0, 5);
  g.set_edge(3, 0, 4);
  g.set_edge(4, 0, 5);
  g.set_edge(5, 0, 6);
  g.set_edge(6, 0, 3);
  g.set_edge(6, 0, 7);
  g.set_edge(7, 0, 9);
  g.set_edge(7, 0, 8);
  g.set_edge(8, 0, 9);
  g.set_edge(8, 0, 1);
  AdjacencyList g1(10);
  g.girvan_newman(g1);
  // check size
  ASSERT_EQ(g1.edge_count(), 12);
  //checking the edge values
  ASSERT_EQ(g1.has_edge(0,1) or g1.has_edge(1,0), true);
  ASSERT_EQ(g1.has_edge(1,2) or g1.has_edge(2,1), true);
  ASSERT_EQ(g1.has_edge(2,0) or g1.has_edge(0,2), true);
  ASSERT_EQ(g1.has_edge(4,3) or g1.has_edge(3,4), true);
  ASSERT_EQ(g1.has_edge(5,4) or g1.has_edge(4,5), true);
  ASSERT_EQ(g1.has_edge(3,5) or g1.has_edge(5,3), true);
  ASSERT_EQ(g1.has_edge(6,5) or g1.has_edge(5,6), true);
  ASSERT_EQ(g1.has_edge(3,6) or g1.has_edge(6,3), true);
  ASSERT_EQ(g1.has_edge(7,6) or g1.has_edge(6,7), true);
  ASSERT_EQ(g1.has_edge(7,8) or g1.has_edge(8,8), true);
  ASSERT_EQ(g1.has_edge(9,7) or g1.has_edge(7,9), true);
  ASSERT_EQ(g1.has_edge(8,9) or g1.has_edge(9,8), true);
  //checking communities
  Map communities;
  g1.bfs_connected_components(communities);
  int community1 = communities[0];
  int community2 = communities[3];
  ASSERT_NE(community1, community2);
  //ccommunity 1
  ASSERT_EQ(community1, communities[1]);
  ASSERT_EQ(community1, communities[2]);
  //ccommunity 2
  ASSERT_EQ(community2, communities[4]);
  ASSERT_EQ(community2, communities[5]);
  ASSERT_EQ(community2, communities[6]);
  ASSERT_EQ(community2, communities[7]);
  ASSERT_EQ(community2, communities[8]);
  ASSERT_EQ(community2, communities[9]);
}

TEST(AdjacencyListTest, XLargeGrivanNewmanTest) {
  AdjacencyList g(14);
  g.set_edge(0, 0, 1);
  g.set_edge(1, 0, 2);
  g.set_edge(2, 0, 0);
  g.set_edge(1, 0, 3);
  g.set_edge(3, 0, 5);
  g.set_edge(3, 0, 4);
  g.set_edge(4, 0, 5);
  g.set_edge(5, 0, 6);
  g.set_edge(6, 0, 3);
  g.set_edge(6, 0, 7);
  g.set_edge(7, 0, 9);
  g.set_edge(7, 0, 8);
  g.set_edge(8, 0, 9);
  g.set_edge(8, 0, 1);
  g.set_edge(7, 0, 10);
  g.set_edge(7, 0, 11);
  g.set_edge(9, 0, 10);
  g.set_edge(10, 0, 11);
  g.set_edge(2, 0, 13);
  g.set_edge(13, 0, 12);
  g.set_edge(12, 0, 3);
  AdjacencyList g1(14);
  g.girvan_newman(g1);
  // check size
  ASSERT_EQ(g1.edge_count(), 18);
  //checking the edge values
  ASSERT_EQ(g1.has_edge(0,1) or g1.has_edge(1,0), true);
  ASSERT_EQ(g1.has_edge(1,2) or g1.has_edge(2,1), true);
  ASSERT_EQ(g1.has_edge(2,0) or g1.has_edge(0,2), true);
  ASSERT_EQ(g1.has_edge(4,3) or g1.has_edge(3,4), true);
  ASSERT_EQ(g1.has_edge(5,4) or g1.has_edge(4,5), true);
  ASSERT_EQ(g1.has_edge(3,5) or g1.has_edge(5,3), true);
  ASSERT_EQ(g1.has_edge(6,5) or g1.has_edge(5,6), true);
  ASSERT_EQ(g1.has_edge(3,6) or g1.has_edge(6,3), true);
  ASSERT_EQ(g1.has_edge(7,8) or g1.has_edge(8,8), true);
  ASSERT_EQ(g1.has_edge(9,7) or g1.has_edge(7,9), true);
  ASSERT_EQ(g1.has_edge(8,9) or g1.has_edge(9,8), true);
  //checking communities
  Map communities;
  g1.bfs_connected_components(communities);
  int community1 = communities[0];
  int community2 = communities[7];
  ASSERT_NE(community1, community2);
  //ccommunity 1
  ASSERT_EQ(community1, communities[1]);
  ASSERT_EQ(community1, communities[2]);
  ASSERT_EQ(community1, communities[3]);
  ASSERT_EQ(community1, communities[4]);
  ASSERT_EQ(community1, communities[5]);
  ASSERT_EQ(community1, communities[6]);
  ASSERT_EQ(community1, communities[12]);
  ASSERT_EQ(community1, communities[13]);
  ASSERT_EQ(community1, communities[6]);
  //ccommunity 2
  ASSERT_EQ(community2, communities[7]);
  ASSERT_EQ(community2, communities[8]);
  ASSERT_EQ(community2, communities[9]);
  ASSERT_EQ(community2, communities[10]);
  ASSERT_EQ(community2, communities[11]);
  
}

TEST(AdjacencyListTest, TriTippedGrivanNewmanTest) {
  AdjacencyList g(9);
  g.set_edge(0, 0, 1);
  g.set_edge(1, 0, 2);
  g.set_edge(2, 0, 0);
  g.set_edge(1, 0, 3);
  g.set_edge(3, 0, 5);
  g.set_edge(3, 0, 4);
  g.set_edge(4, 0, 5);
  g.set_edge(5, 0, 7);
  g.set_edge(7, 0, 8);
  g.set_edge(6, 0, 7);
  g.set_edge(6, 0, 8);
  g.set_edge(6, 0, 2);
  AdjacencyList g1(9);
  g.girvan_newman(g1);
  // check size
  ASSERT_EQ(g1.edge_count(), 9);
  //checking the edge values
  ASSERT_EQ(g1.has_edge(0,1) or g1.has_edge(1,0), true);
  ASSERT_EQ(g1.has_edge(1,2) or g1.has_edge(2,1), true);
  ASSERT_EQ(g1.has_edge(2,0) or g1.has_edge(0,2), true);
  ASSERT_EQ(g1.has_edge(4,3) or g1.has_edge(3,4), true);
  ASSERT_EQ(g1.has_edge(5,4) or g1.has_edge(4,5), true);
  ASSERT_EQ(g1.has_edge(3,5) or g1.has_edge(5,3), true);
  //checking communities
  Map communities;
  g1.bfs_connected_components(communities);
  int community1 = communities[0];
  int community2 = communities[3];
  int community3 = communities[6];
  ASSERT_NE(community1, community2);
  ASSERT_NE(community1, community3);
  ASSERT_NE(community1, community3);
  //community1
  ASSERT_EQ(community1, communities[1]);
  ASSERT_EQ(community1, communities[2]);
  //community 2
  ASSERT_EQ(community2, communities[4]);
  ASSERT_EQ(community2, communities[5]);
  //community 3
  ASSERT_EQ(community3, communities[7]);
  ASSERT_EQ(community3, communities[8]);
  
}

//I had to drastically change the adjacency list functions to work with doubles for the edge length, i did not make that change for the adjacency matrix
//so there are no tests for adjacency matrix

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

