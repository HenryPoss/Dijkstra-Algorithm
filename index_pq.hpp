#ifndef INDEX_PRIORITY_QUEUE_HPP_
#define INDEX_PRIORITY_QUEUE_HPP_

#include <vector>
#include <algorithm>

template <typename T>
class IndexPriorityQueue {
 private:
  // vector to hold priorities.  
  // priorities.at(i) is the priority associated to index i
  // the heap will only contain indices, we look up priorities
  // in this vector as needed
  std::vector<T> priorities {};
  // priorityQueue stores indices: priorityQueue.at(i) is an index
  // priorityQueue functions as the heap and is heap ordered: 
  // priorities.at(priorityQueue.at(i)) <= priorities.at(priorityQueue.at(2 * i)) 
  // priorities.at(priorityQueue.at(i)) <= priorities.at(priorityQueue.at(2 * i) + 1) 
  std::vector<int> priorityQueue {};
  // indexToPosition.at(i) is the position in priorityQueue of index i
  // priorityQueue.at(indexToPosition.at(i)) = i
  // indexToPosition.at(priorityQueue.at(j)) = j
  std::vector<int> indexToPosition {};
  int size_ = 0;

 public:
  explicit IndexPriorityQueue(int);
  void push(const T&, int);
  void pop();
  void erase(int);
  bool contains(int) const;
  void changeKey(const T&, int);
  std::pair<T, int> top() const;
  bool empty() const;
  int size() const;

 private:
  // TODO: you may want to add your own member functions. swim? sink?

  // swim and skin functions move the indexes around according the priority.
  void swim(int i);
  void sink(int i);
  // the following variable holds the maximum number of elements allowed. This is required to check that a user has given an index that is possible.
  int maxSize_;
};

// Useful helper functions
int leftChild(int i) {
  return 2*i;
}

int rightChild(int i) {
  return 2*i + 1;
}

int parent(int i) {
  return i/2;
}

// IndexPriorityQueue member functions
template <typename T>
IndexPriorityQueue<T>::IndexPriorityQueue(int N) : priorities {std::vector<T>(N)},
                                                   priorityQueue {std::vector<int>(N + 1)},
                                                   indexToPosition {std::vector<int>(N)},
                                                   maxSize_ {N} {
  // Initialise all indexToPosition values as -1 to indicate that the particular element doesn't exist.
  for (int i = 0; i < N; i++) {
    indexToPosition.at(i) = -1;
  }
}

template <typename T>
bool IndexPriorityQueue<T>::empty() const {
  return (size_ == 0);
}

template <typename T>
int IndexPriorityQueue<T>::size() const {
  return size_;
}

template <typename T>
void IndexPriorityQueue<T>::push(const T& priority, int index) {
  // If index already in index priority queue then return.
  if (contains(index)) {
    return;
  }
  ++size_;
  priorities.at(index) = priority;
  priorityQueue.at(size_) = index;
  indexToPosition.at(index) = size_;
  // swim so that the index gets moved to its correct position.
  swim(size_);
}

template <typename T>
void IndexPriorityQueue<T>::pop() {
  // If index priority queue is empty then return
  if (empty()) {
    return;
  }
  // swap first and last elements in priorityQueue
  std::swap(priorityQueue.at(1), priorityQueue.at(size_));
  // Update indexToPosition array to reflect new positions
  indexToPosition.at(priorityQueue.at(1)) = 1;
  indexToPosition.at(priorityQueue.at(size_)) = size_;
  --size_;
  // sink so that the elements are in the correct order
  sink(1);
  // Update indexToPosition and and priorityQueue to 'remove' the value
  indexToPosition.at(priorityQueue.at(size_ + 1)) = -1;
  priorityQueue.at(size_ + 1) = -1;
}

template <typename T>
void IndexPriorityQueue<T>::erase(int index) {
  // If the index priority queue does nto contain the index then return
  if (!contains(index)) {
    return;
  }
  // find the position of index
  int i = indexToPosition.at(index);
  // swapping the index with the last element and updating IndexToPosition
  std::swap(priorityQueue.at(i), priorityQueue.at(size_));
  indexToPosition.at(priorityQueue.at(i)) = i;
  indexToPosition.at(priorityQueue.at(size_)) = size_;
  --size_;
  // swimming and sinking to arrange the indicies
  swim(i);
  sink(i);
  // Update indexToPosition and and priorityQueue to 'remove' the value
  priorityQueue.at(size_ + 1) = -1;
  indexToPosition.at(index) = -1;
}

template <typename T>
std::pair<T, int> IndexPriorityQueue<T>::top() const {
  // Only attempt to find the top element if index priority queue is not empty
  if (!empty()) {
    return {priorities.at(priorityQueue.at(1)), priorityQueue.at(1)};
  }
  return {T {}, 0};
}

// if vertex i is not present, insert it with key
// otherwise change the associated key value of i to key
template <typename T>
void IndexPriorityQueue<T>::changeKey(const T& key, int index) {
  // index has to be a valid value
  if (index < 0 || index >= maxSize_) {
    return;
  }

  // If index is already present in index priority queue, then remove it and add it again with the new key
  if (contains(index)) {
    erase(index);
    push(key, index);
    // If index not already present then simply push it
  } else {
    push(key, index);
  }
}

template <typename T>
bool IndexPriorityQueue<T>::contains(int index) const {
  // index has to be a valid value
  if (index < 0 || index >= maxSize_) {
    return false;
  }
  // If indexToPosition of the index is equal to -1 then the index does not exist
  return (indexToPosition.at(index) != -1);
}

template <typename T>
void IndexPriorityQueue<T>::swim(int i) {
  // Find parent of index
  int p = parent(i);
  // Test of parent is greater than 0 (as 1 is the first index) and test if the parent has a higher priority value than the child
  while(p > 0 && priorities.at(priorityQueue.at(i)) < priorities.at(priorityQueue.at(p))) {
    // if priority value of parent is higher, then we swap the indicies and update their indexToPosition values
    std::swap(priorityQueue.at(p), priorityQueue.at(i));
    indexToPosition.at(priorityQueue.at(p)) = p;
    indexToPosition.at(priorityQueue.at(i)) = i;
    // update index and parent to refelct their new locations
    i = p;
    p = parent(i);
  }
}

template <typename T>
void IndexPriorityQueue<T>::sink(int i) {
  // current represents the index we start with
  int current {i};
  // check if the index priority queue has a left child
  while (leftChild(current) <= size()) {
    int left = leftChild(current);
    int right = rightChild(current);
    // Assume left child is the smaller child for now
    int swapWith = left;
    // If right child exists and its priority is smaller than left child, then make right child the index to swap with
    if (right <= size() && priorities.at(priorityQueue.at(left)) > priorities.at(priorityQueue.at(right))) {
      swapWith = right;
    }
    // If the priority of the current index is smaller (or equal to) than the priority of the smaller child then return as we do not need to go any further
    if (priorities.at(priorityQueue.at(current)) <= priorities.at(priorityQueue.at(swapWith))) {
      return;
    }
    // If the priority of the smaller child is smaller than the priority of the current index then the following code will run
    // Swap the indicies and update the indexToPosition array
    std::swap(priorityQueue.at(current), priorityQueue.at(swapWith));
    indexToPosition.at(priorityQueue.at(current)) = current;
    indexToPosition.at(priorityQueue.at(swapWith)) = swapWith;
    // update current with swapWith, as that is now the position of the index
    current = swapWith;
  }
}


#endif      // INDEX_PRIORITY_QUEUE_HPP_
