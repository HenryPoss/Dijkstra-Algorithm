#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include "graph.hpp"
#include "my_integer.hpp"

// Example Test of shortest path algorithm from tutorial
TEST(DijkstraTest, distanceFrom0LazyInt) {
  Graph<int> G {"tinyEWD.txt"};
  std::vector<int> result = singleSourceLazyDistance<int>(G, 0);
  EXPECT_EQ(result[0], 0);
  EXPECT_EQ(result[1], 105);
  EXPECT_EQ(result[2], 26);
  EXPECT_EQ(result[3], 99);
  EXPECT_EQ(result[4], 38);
  EXPECT_EQ(result[5], 73);
  EXPECT_EQ(result[6], 151);
  EXPECT_EQ(result[7], 60);
}

// Example showing how to use MyInteger class to get performance data
TEST(DijkstraTest, distanceFrom0LazyMyInteger) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  MyInteger::clearCounts();
  std::vector<MyInteger> result = singleSourceLazyDistance<MyInteger>(G, 0);
  // see how many operations on MyInteger the algorithm made
  MyInteger::printCounts();
  EXPECT_EQ(result[0], MyInteger {0});
  EXPECT_EQ(result[1], MyInteger {105});
  EXPECT_EQ(result[2], MyInteger {26});
  EXPECT_EQ(result[3], MyInteger {99});
  EXPECT_EQ(result[4], MyInteger {38});
  EXPECT_EQ(result[5], MyInteger {73});
  EXPECT_EQ(result[6], MyInteger {151});
  EXPECT_EQ(result[7], MyInteger {60});
}

// You can generate some random graphs to help in your testing
// The graph has N vertices and p is the probability there is an
// edge between any two vertices. 
// You can vary seed to get different graphs
Graph<int> randomGraph(int N, unsigned seed, double p = 0.5) {
  std::mt19937 mt {seed};
  // set up random number generator that is 1 with probability p and
  // 0 with probability 1-p
  std::binomial_distribution<int> heads {1, p};
  // Set the minimum and maximum edge weight here
  const int minWeight {1};
  const int maxWeight {10};
  std::uniform_int_distribution<int> weight {minWeight, maxWeight};
  Graph<int> G {N};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (heads(mt)) {
        // add edge between i and j with random weight 
        // between minWeight and maxWeight
        G.addEdge(i, j, weight(mt));
      }
    }
  }
  return G;
}

Graph<double> randomGraphDouble(int N, unsigned seed, double p = 0.5) {
  std::mt19937 mt {seed};
  // set up random number generator that is 1 with probability p and
  // 0 with probability 1-p
  std::binomial_distribution<int> heads {1, p};
  // Set the minimum and maximum edge weight here
  const double minWeight {1};
  const double maxWeight {10};
  std::uniform_real_distribution<double> weight {minWeight, maxWeight};
  Graph<double> G {N};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (heads(mt)) {
        // add edge between i and j with random weight 
        // between minWeight and maxWeight
        G.addEdge(i, j, weight(mt));
      }
    }
  }
  return G;
}

Graph<MyInteger> randomGraphMyInteger(int N, unsigned seed, double p = 0.5) {
  std::mt19937 mt {seed};
  // set up random number generator that is 1 with probability p and
  // 0 with probability 1-p
  std::binomial_distribution<int> heads {1, p};
  // Set the minimum and maximum edge weight here
  const int minWeight {1};
  const int maxWeight {10};
  std::uniform_int_distribution<int> weight {minWeight, maxWeight};
  Graph<MyInteger> G {N};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (heads(mt)) {
        // add edge between i and j with random weight 
        // between minWeight and maxWeight
        G.addEdge(i, j, MyInteger {weight(mt)});
      }
    }
  }
  return G;
}

TEST(lazyDijkstra, lazyDijkstraNormalIntSubGraph) {
  Graph<int> G = randomGraph(20, 100);
  Graph<int> shortestPathTree = singleSourceLazy(G, 0);

  EXPECT_TRUE(isSubgraph(shortestPathTree, G));
}

TEST(lazyDijkstra, lazyDijkstraNormalIntTreePlusIsolated) {
  Graph<int> G = randomGraph(20, 100);
  Graph<int> shortestPathTree = singleSourceLazy(G, 0);

  EXPECT_TRUE(isTreePlusIsolated(shortestPathTree, 0));
}

TEST(lazyDijkstra, lazyDijkstraNormalIntShortestPath) {
  Graph<int> G = randomGraph(20, 100);
  Graph<int> shortestPathTree = singleSourceLazy(G, 0);

  std::vector<int> shortestPathsVector = singleSourceLazyDistance(G, 0);
  EXPECT_EQ(pathLengthsFromRoot(shortestPathTree, 0), shortestPathsVector);

  EXPECT_TRUE(allEdgesRelaxed(shortestPathsVector, G, 0));
}

TEST(lazyDijkstra, lazyDijkstraDoubleSubGraph) {
  Graph<double> G = randomGraphDouble(10, 100);
  Graph<double> shortestPathTree = singleSourceLazy(G, 0);
  
  EXPECT_TRUE(isSubgraph(shortestPathTree, G));
}

TEST(lazyDijkstra, lazyDijkstraDoubleTreePlusIsolated) {
  Graph<double> G = randomGraphDouble(10, 100);
  Graph<double> shortestPathTree = singleSourceLazy(G, 0);

  EXPECT_TRUE(isTreePlusIsolated(shortestPathTree, 0));
}

TEST(lazyDijkstra, lazyDijkstraDoubleShortestPath) {
  Graph<double> G = randomGraphDouble(10, 100);
  Graph<double> shortestPathTree = singleSourceLazy(G, 0);

  std::vector<double> shortestPathsVector = singleSourceLazyDistance(G, 0);
  EXPECT_EQ(pathLengthsFromRoot(shortestPathTree, 0), shortestPathsVector);

  EXPECT_TRUE(allEdgesRelaxed(shortestPathsVector, G, 0));
}

TEST(lazyDijkstra, lazyDijkstraMyIntegerSubGraph) {
  Graph<MyInteger> G = randomGraphMyInteger(10, 100);
  Graph<MyInteger> shortestPathTree = singleSourceLazy(G, 0);
  
  EXPECT_TRUE(isSubgraph(shortestPathTree, G));
}

TEST(lazyDijkstra, lazyDijkstraMyIntegerTreePlusIsolated) {
  Graph<MyInteger> G = randomGraphMyInteger(10, 100);
  Graph<MyInteger> shortestPathTree = singleSourceLazy(G, 0);

  EXPECT_TRUE(isTreePlusIsolated(shortestPathTree, 0));
}

TEST(lazyDijkstra, lazyDijkstraMyIntegerShortestPath) {
  MyInteger::clearCounts();

  Graph<MyInteger> G = randomGraphMyInteger(10, 100);
  Graph<MyInteger> shortestPathTree = singleSourceLazy(G, 0);

  MyInteger::printCounts();

  std::vector<MyInteger> shortestPathsVector = singleSourceLazyDistance(G, 0);
  EXPECT_EQ(pathLengthsFromRoot(shortestPathTree, 0), shortestPathsVector);

  EXPECT_TRUE(allEdgesRelaxed(shortestPathsVector, G, 0));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraNormalIntSubGraph) {
  Graph<int> G = randomGraph(20, 100);
  Graph<int> shortestPathTree = singleSourceIndex(G, 0);

  EXPECT_TRUE(isSubgraph(shortestPathTree, G));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraNormalIntTreePlusIsolated) {
  Graph<int> G = randomGraph(20, 100);
  Graph<int> shortestPathTree = singleSourceIndex(G, 0);

  EXPECT_TRUE(isTreePlusIsolated(shortestPathTree, 0));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraNormalIntShortestPath) {
  Graph<int> G = randomGraph(20, 100);
  Graph<int> shortestPathTree = singleSourceIndex(G, 0);

  std::vector<int> shortestPathsVector = singleSourceLazyDistance(G, 0);
  EXPECT_EQ(pathLengthsFromRoot(shortestPathTree, 0), shortestPathsVector);

  EXPECT_TRUE(allEdgesRelaxed(shortestPathsVector, G, 0));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraDoubleSubGraph) {
  Graph<double> G = randomGraphDouble(10, 100);
  Graph<double> shortestPathTree = singleSourceIndex(G, 0);
  
  EXPECT_TRUE(isSubgraph(shortestPathTree, G));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraDoubleTreePlusIsolated) {
  Graph<double> G = randomGraphDouble(10, 100);
  Graph<double> shortestPathTree = singleSourceIndex(G, 0);

  EXPECT_TRUE(isTreePlusIsolated(shortestPathTree, 0));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraDoubleShortestPath) {
  Graph<double> G = randomGraphDouble(10, 100);
  Graph<double> shortestPathTree = singleSourceIndex(G, 0);

  std::vector<double> shortestPathsVector = singleSourceLazyDistance(G, 0);
  EXPECT_EQ(pathLengthsFromRoot(shortestPathTree, 0), shortestPathsVector);

  EXPECT_TRUE(allEdgesRelaxed(shortestPathsVector, G, 0));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraMyIntegerSubGraph) {
  Graph<MyInteger> G = randomGraphMyInteger(10, 100);
  Graph<MyInteger> shortestPathTree = singleSourceIndex(G, 0);
  
  EXPECT_TRUE(isSubgraph(shortestPathTree, G));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraMyIntegerTreePlusIsolated) {
  Graph<MyInteger> G = randomGraphMyInteger(10, 100);
  Graph<MyInteger> shortestPathTree = singleSourceIndex(G, 0);
  
  EXPECT_TRUE(isTreePlusIsolated(shortestPathTree, 0));
}

TEST(IndexPriorityQueueDijikstra, IndexPrioirtyQueueDijkstraMyIntegerShortestPath) {
  MyInteger::clearCounts();
  
  Graph<MyInteger> G = randomGraphMyInteger(10, 100);
  Graph<MyInteger> shortestPathTree = singleSourceIndex(G, 0);

  MyInteger::printCounts();
  
  std::vector<MyInteger> shortestPathsVector = singleSourceLazyDistance(G, 0);
  EXPECT_EQ(pathLengthsFromRoot(shortestPathTree, 0), shortestPathsVector);

  EXPECT_TRUE(allEdgesRelaxed(shortestPathsVector, G, 0));
}


int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
