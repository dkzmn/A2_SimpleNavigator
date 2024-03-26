#include "s21_graph.h"

namespace s21 {
Graph::Graph() : vertices_(1), directed_(false), weighted_(false) {
  adj_matrix_ = std::vector<std::vector<uint32_t>>(
      vertices_, std::vector<uint32_t>(vertices_));
}

Graph::Graph(uint32_t vertices)
    : vertices_(vertices), directed_(false), weighted_(false) {
  if (vertices < 1)
    throw std::out_of_range("There should be at least one vertix");
  adj_matrix_ = std::vector<std::vector<uint32_t>>(
      vertices_, std::vector<uint32_t>(vertices_));
}

uint32_t Graph::GetVertices() const noexcept { return vertices_; }

std::vector<std::vector<uint32_t>> Graph::GetMatrix() const noexcept {
  return adj_matrix_;
}

void Graph::LoadGraphFromFile(std::string filename) {
  std::ifstream input(filename);
  if (input) {
    uint32_t vertices, temp;
    input >> vertices;
    if (vertices < 1)
      throw std::out_of_range("There should be at least one row and/or column");
    vertices_ = vertices;
    adj_matrix_ = std::vector<std::vector<uint32_t>>(
        vertices_, std::vector<uint32_t>(vertices_));
    for (uint32_t i = 0; i != vertices_; ++i) {
      for (uint32_t j = 0; j != vertices_; ++j) {
        if (!input.eof()) {
          input >> temp;
          if (temp > 1) weighted_ = true;
          adj_matrix_[i][j] = static_cast<uint32_t>(temp);
        } else {
          throw std::length_error("Not enough data in the file");
        }
      }
    }
    input.close();
    if (GraphIsDirected()) directed_ = true;
  } else {
    throw std::invalid_argument("No such file.\n");
  }
}

void Graph::ExportGraphToDot(std::string filename) {
  std::ofstream output;
  output.open(filename);
  if (output.is_open()) {
    if (directed_)
      output << "digraph {" << std::endl;
    else
      output << "graph {" << std::endl;
    for (uint32_t i = 0; i < vertices_; ++i) {
      for (uint32_t j = 0; j < vertices_; ++j) {
        if (adj_matrix_[i][j] > 0) {
          output << "\t" << (i + 1);
          if (directed_)
            output << " -> ";
          else
            output << " -- ";
          output << (j + 1) << " [weight=" << adj_matrix_[i][j] << "];"
                 << std::endl;
        }
      }
    }
    output << "}" << std::endl;
    output.close();
  } else
    throw std::logic_error("Couldn't open file.");
}

void Graph::SaveGraphToFile(std::string filename) {
  std::ofstream output;
  output.open(filename);
  if (output.is_open()) {
    output << vertices_ << std::endl;
    for (uint32_t i = 0; i < vertices_; ++i) {
      for (uint32_t j = 0; j < vertices_ - 1; ++j) {
        output << adj_matrix_[i][j] << "  ";
      }
      output << adj_matrix_[i][vertices_ - 1] << std::endl;
    }
    output.close();
  } else
    throw std::logic_error("Couldn't open file.");
}

bool Graph::GraphIsDirected() {
  bool result = false;
  uint32_t edjes_count = 0, reverse_edjes_count = 0;
  for (uint32_t i = 0; i != vertices_; ++i) {
    for (uint32_t j = 0; j != vertices_; ++j) {
      if (adj_matrix_[i][j] > 0) {
        ++edjes_count;
        if (adj_matrix_[i][j] == adj_matrix_[j][i]) ++reverse_edjes_count;
      }
    }
  }
  if (edjes_count != reverse_edjes_count) result = true;
  return result;
}

bool Graph::GraphIsFull() {
  bool result = true;
  for (uint32_t i = 0; i != vertices_; ++i) {
    for (uint32_t j = 0; j != vertices_; ++j) {
      if (adj_matrix_[i][j] == 0 && i != j) result = false;
    }
  }
  return result;
}

bool Graph::GraphIsWeighted() {
  bool result = false;
  for (uint32_t i = 0; i != vertices_; ++i) {
    for (uint32_t j = 0; j != vertices_; ++j) {
      if (adj_matrix_[i][j] > 1) result = true;
    }
  }
  return result;
}

}  // namespace s21