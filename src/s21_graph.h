#ifndef A2_SIMPLENAVIGATOR_V1_0_1_S21_GRAPH_H_
#define A2_SIMPLENAVIGATOR_V1_0_1_S21_GRAPH_H_

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace s21 {
class Graph {
 public:
  // Constructors:
  Graph();
  Graph(uint32_t vertices);

  // Getters:
  uint32_t GetVertices() const noexcept;
  std::vector<std::vector<uint32_t>> GetMatrix() const noexcept;

  // Public methods:
  void LoadGraphFromFile(std::string filename);
  void ExportGraphToDot(std::string filename);

  // Helpers:
  bool GraphIsDirected();
  bool GraphIsFull();
  bool GraphIsWeighted();
  void SaveGraphToFile(std::string filename);

 private:
  // Private fields:
  uint32_t vertices_;
  std::vector<std::vector<uint32_t>> adj_matrix_;
  bool directed_, weighted_;
};
}  // namespace s21

#endif  // A2_SIMPLENAVIGATOR_V1_0_1_S21_GRAPH_H_