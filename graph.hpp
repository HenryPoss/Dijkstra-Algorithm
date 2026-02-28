#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <utility>
#include <functional>
#include <vector>
#include <string>
#include <queue>
#include <set>
#include <unordered_map>
#include <limits>
#include "my_integer.hpp"
#include "index_pq.hpp"

template <typename T>
class Graph {
 private:
  std::vector<std::unordered_map<int, T> > adjList {};
  int numVertices {};

 public:
  // empty graph with N vertices
  explicit Graph(int N);

  // construct graph from edge list in filename
  explicit Graph(const std::string& filename);

  // add an edge directed from vertex i to vertex j with given weight
  void addEdge(int i, int j, T weight);

  // removes edge from vertex i to vertex j
  void removeEdge(int i, int j);

  // is there an edge from vertex i to vertex j?
  bool isEdge(int i, int j) const;

  // return weight of edge from i to j
  // will throw an exception if there is no edge from i to j
  T getEdgeWeight(int i, int j) const;

  // returns number of vertices in the graph
  int size() const;

  // alias a const iterator to our adjacency list type to iterator
  using iterator = 
  typename std::vector<std::unordered_map<int, T> >::const_iterator;

  // cbegin returns const iterator pointing to first element of adjList
  iterator begin() const {
    return adjList.cbegin();
  }

  iterator end() const {
    return adjList.cend();
  }

  // return iterator to a particular vertex
  iterator neighbours(int a) const {
    return adjList.begin() + a;
  }
};

template <typename T>
Graph<T>::Graph(int N) : adjList(N), numVertices {N} {}

template <typename T>
Graph<T>::Graph(const std::string& inputFile) {
  std::ifstream infile {inputFile};
  if (!infile) {
    std::cerr << inputFile << " could not be opened\n";
    return;
  }
  // first line has number of vertices
  infile >> numVertices;
  adjList.resize(numVertices);
  int i {};
  int j {};
  double weight {};
  // assume each remaining line is of form
  // origin dest weight
  while (infile >> i >> j >> weight) {
    addEdge(i, j, static_cast<T>(weight));
  }
}

template <typename T>
int Graph<T>::size() const {
  return numVertices;
}

template <typename T>
void Graph<T>::addEdge(int i, int j, T weight) {
  if (i < 0 or i >= numVertices or j < 0 or j >= numVertices) {
    throw std::out_of_range("invalid vertex number");
  }
  adjList[i].insert({j, weight});
}

template <typename T>
void Graph<T>::removeEdge(int i, int j) {
  // check if i and j are valid
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    adjList[i].erase(j);
  }
}

template <typename T>
bool Graph<T>::isEdge(int i, int j) const {
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    return adjList.at(i).contains(j);
  }
  return false;
}

template <typename T>
T Graph<T>::getEdgeWeight(int i, int j) const {
  return adjList.at(i).at(j);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Graph<T>& G) {
  for (int i = 0; i < G.size(); ++i) {
    out << i << ':';
    for (const auto& [neighbour, weight] : *(G.neighbours(i))) {
      out << " (" << i << ", " << neighbour << ")[" << weight << ']';
    }
    out << '\n';
  }
  return out;
}

template <typename T>
bool isSubgraph(const Graph<T>& H, const Graph<T>& G) {
  if (H.size() > G.size()) {
    return false;
  }

  for (int i = 0; i < H.size(); i++) {
    for (auto [neighbour, weight] : *(H.neighbours(i))) {
      if (!G.isEdge(i, neighbour) || (G.getEdgeWeight(i, neighbour) != weight)) {
        return false;
      }
    }
  }
  return true;
}

template <typename T>
bool isTreePlusIsolated(const Graph<T>& G, int root) {
  std::queue<int> visitQueue;
  std::vector<bool> marked(G.size());

  visitQueue.push(root);
  
  while (visitQueue.empty() != true) {
    int x = visitQueue.front();
    visitQueue.pop();
    if (marked.at(x)) {
      return false;
    }
    marked.at(x) = true;
    for (auto [neighbour, weight] : *(G.neighbours(x))) {
      visitQueue.push(neighbour);
    }
  }
  for (int i = 0; i < G.size(); i++) {
    for (auto [neighbour, weight] : *(G.neighbours(i))) {
      if (!marked.at(neighbour)) {
        return false;
      }
    } 
  }
  return true;
}

template <typename T>
std::vector<T> pathLengthsFromRoot(const Graph<T>& tree, int root) {
  std::vector<T> bestDistanceTo(tree.size());
  std::queue<std::pair<int, T>> visitQueue {};

  bestDistanceTo.at(root) = T {};
  visitQueue.push({root, T {}});  
  
  while (visitQueue.empty() != true) {
    int currentIndex = visitQueue.front().first;
    T currentDistance = visitQueue.front().second;
    visitQueue.pop();
    for (auto [neighbour, weight] : *(tree.neighbours(currentIndex))) {
      bestDistanceTo.at(neighbour) = currentDistance + weight;
      visitQueue.push({neighbour, currentDistance + weight});   
    }
  }
  
  return bestDistanceTo;
}

template <typename T>
bool allEdgesRelaxed(const std::vector<T>& bestDistanceTo, const Graph<T>& G, 
                      int source) {
  if (bestDistanceTo.at(source) != T {}) {
     return false;
  }

  for (int i = 0; i < G.size(); i++) {
    for (auto [neighbour, weight] : *(G.neighbours(i))) {
      if (bestDistanceTo.at(neighbour) > bestDistanceTo.at(i) + weight) {
        return false;
      }
    }
  }

  return true;
}
// End of functions from Graph class

// return a number of type T to stand in for "infinity"
template <typename T>
T infinity() {
  if (std::numeric_limits<T>::has_infinity) {
    return std::numeric_limits<T>::infinity();
  } else if constexpr (std::is_same_v<T, MyInteger>) {
    return MyInteger {std::numeric_limits<int>::max()};
  } else {
    return std::numeric_limits<T>::max();
  }
}

// lazy solution as in Tutorial Week 10
template <typename T>
std::vector<T> singleSourceLazyDistance(const Graph<T>& G, int source) {
  // alias the long name for a minimum priority queue holding
  // objects of type DistAndVertex
  using DistAndVertex = std::pair<T, int>;
  using minPQ = std::priority_queue<DistAndVertex,
                                  std::vector<DistAndVertex>,
                                  std::greater<DistAndVertex> >;
  minPQ queue {};
  queue.push({T {}, source});
  // record best distance to vertex found so far
  int N = G.size();
  std::vector<T> bestDistanceTo(N, infinity<T>());
  bestDistanceTo.at(source) = T {};
  // being in visited means we have already explored a vertex's neighbours
  // the bestDistanceTo for a vertex in visited is the true distance.
  std::vector<bool> visited(N);
  while (!queue.empty()) {
    auto [dist, current] = queue.top();
    queue.pop();
    // as we use a lazy version of Dijkstra a vertex can appear multiple
    // times in the queue.  If we have already visited the vertex we
    // take out of the queue we just go on to the next one
    if (visited.at(current)) {
      continue;
    }
    visited.at(current) = true;
    // relax all outgoing edges of current
    for (const auto& [neighbour, weight] : *(G.neighbours(current))) {
      T distanceViaCurrent = bestDistanceTo.at(current) + weight;
      if (bestDistanceTo.at(neighbour) > distanceViaCurrent) {
        bestDistanceTo.at(neighbour) = distanceViaCurrent;
        // lazy dijkstra: nextPoint could already be in the queue
        // we don't update it with better distance just found.
        queue.push({distanceViaCurrent, neighbour});
      }
    }
  }
  return bestDistanceTo;
}
// Implement your solution using an index priority queue here
template <typename T>
Graph<T> singleSourceIndex(const Graph<T>& G, int source) {
  IndexPriorityQueue<T> queue(G.size());
  queue.push(T {}, source);

  std::vector<T> bestDistanceTo(G.size(), infinity<T>());
  bestDistanceTo.at(source) = T {};
  std::unordered_map<int, int> path {};
  std::unordered_map<int, T> distances {};

  while (!queue.empty()) {
    auto [distance, current] = queue.top();
    queue.pop();
    for (auto [neighbour, weight] : *(G.neighbours(current))) {
      if (bestDistanceTo.at(neighbour) > bestDistanceTo.at(current) + weight) {
        bestDistanceTo.at(neighbour) = bestDistanceTo.at(current) + weight;
        path[neighbour] = current;
        distances[neighbour] = weight;
        queue.changeKey(bestDistanceTo.at(current) + weight, neighbour);
      }
    }
  }
  
  Graph<T> shortestPathTree(G.size());

  for (auto it = path.begin(); it != path.end(); it++) {
    shortestPathTree.addEdge(it->second, it->first, distances[it->first]);
  }

  // If shortest path tree graph does not meet requirements of shortest path tree then return an empty graph.
  if (!isSubgraph(shortestPathTree, G) || !isTreePlusIsolated(shortestPathTree, source)) {
    return Graph<T> {G.size()};
  }
  std::vector<T> pathLengths = pathLengthsFromRoot(shortestPathTree, source);

  // Ensure no other shorter paths. If so, return empty graph
  if (!allEdgesRelaxed(pathLengths, G, source)) {
    return Graph<T> {G.size()};
  }
  
  return shortestPathTree;
}

// Implement your lazy solution using std::priority_queue here
template <typename T>
Graph<T> singleSourceLazy(const Graph<T>& G, int source) {
  using DistanceAndVertex = std::pair<T, int>;
  std::priority_queue<DistanceAndVertex, std::vector<DistanceAndVertex>, std::greater<DistanceAndVertex>> queue {};
  queue.push({T {}, source});

  std::vector<T> bestDistanceTo(G.size(), infinity<T>());
  bestDistanceTo.at(source) = T {};
  std::vector<bool> visited(G.size());
  std::unordered_map<int, int> path {};
  std::unordered_map<int, T> distances {};

  while(!queue.empty()) {
    auto [distance, current] = queue.top();
    queue.pop();
    if (visited.at(current)) {
      continue;
    }
    visited.at(current) = true;
    for (auto [neighbour, weight] : *(G.neighbours(current))) {
      if (bestDistanceTo.at(neighbour) > bestDistanceTo.at(current) + weight) {
        bestDistanceTo.at(neighbour) = bestDistanceTo.at(current) + weight;
        path[neighbour] = current;
        distances[neighbour] = weight;
        queue.push({bestDistanceTo.at(current) + weight, neighbour});
      }
    }
  }
  Graph<T> shortestPathTree(G.size());

  for (auto it = path.begin(); it != path.end(); it++) {
    shortestPathTree.addEdge(it->second, it->first, distances[it->first]);
  }

  // If shortest path tree graph does not meet requirements of shortest path tree then return an empty graph.
  if (!isSubgraph(shortestPathTree, G) || !isTreePlusIsolated(shortestPathTree, source)) {
    return Graph<T> {G.size()};
  }
  
  std::vector<T> pathLengths = pathLengthsFromRoot(shortestPathTree, source);

  // Ensure no other shorter paths. If so, return empty graph
  if (!allEdgesRelaxed(pathLengths, G, source)) {
    return Graph<T> {G.size()};
  }
  
  return shortestPathTree;
}

// Implement your solution using std::set here
template <typename T>
Graph<T> singleSourceSet(const Graph<T>& G, int source) {
  return Graph<T> {G.size()};
}

// put your "best" solution here
// this is the one we will use for performance testing
template <typename T>
Graph<T> singleSourceShortestPaths(const Graph<T>& G, int source) {
  IndexPriorityQueue<T> queue(G.size());
  queue.push(T {}, source);

  std::vector<T> bestDistanceTo(G.size(), infinity<T>());
  bestDistanceTo.at(source) = T {};
  std::unordered_map<int, int> path {};
  std::unordered_map<int, T> distances {};

  while (!queue.empty()) {
    auto [distance, current] = queue.top();
    queue.pop();
    for (auto [neighbour, weight] : *(G.neighbours(current))) {
      if (bestDistanceTo.at(neighbour) > bestDistanceTo.at(current) + weight) {
        bestDistanceTo.at(neighbour) = bestDistanceTo.at(current) + weight;
        path[neighbour] = current;
        distances[neighbour] = weight;
        queue.changeKey(bestDistanceTo.at(current) + weight, neighbour);
      }
    }
  }
  
  Graph<T> shortestPathTree(G.size());

  for (auto it = path.begin(); it != path.end(); it++) {
    shortestPathTree.addEdge(it->second, it->first, distances[it->first]);
  }

  // If shortest path tree graph does not meet requirements of shortest path tree then return an empty graph.
  if (!isSubgraph(shortestPathTree, G) || !isTreePlusIsolated(shortestPathTree, source)) {
    return Graph<T> {G.size()};
  }
  std::vector<T> pathLengths = pathLengthsFromRoot(shortestPathTree, source);

  // Ensure no other shorter paths. If so, return empty graph
  if (!allEdgesRelaxed(pathLengths, G, source)) {
    return Graph<T> {G.size()};
  }
  
  return shortestPathTree;
}

#endif      // GRAPH_HPP_
