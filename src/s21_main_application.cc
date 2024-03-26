#include "s21_main_application.h"

#include <chrono>
#include <iostream>

namespace s21 {

void MainApplication::PrintMenu() {
  std::cout << "=== MAIN MENU ===" << std::endl;
  std::cout << "1. Loading the original graph from a file." << std::endl;
  std::cout << "2. Graph traversal in breadth with output of the result to the "
               "console."
            << std::endl;
  std::cout
      << "3. Graph traversal in depth with output of the result to the console."
      << std::endl;
  std::cout << "4. Searching for the shortest path between any two vertices "
               "and displaying the result to the console."
            << std::endl;
  std::cout
      << "5. Searching for the shortest paths between all pairs of vertices in "
         "the graph with the output of the resulting matrix to the console."
      << std::endl;
  std::cout << "6. Searching for the minimal spanning tree in the graph with "
               "the output of the resulting adjacency matrix to the console."
            << std::endl;
  std::cout << "7. Solving the salesman problem with the output of the "
               "resulting route and its length to the console."
            << std::endl;
  std::cout << "8. Comparison of speed of the three algorithms TSP."
            << std::endl;
  std::cout << "0. Exit" << std::endl;
}

int MainApplication::GetChoice() {
  int choice;
  do {
    std::cout << "Enter your choice [0..8] ";
    std::cin >> choice;
  } while (choice < 0 || choice > 8);
  return choice;
}

bool MainApplication::CheckGraphLoaded() {
  if (g_.GetVertices() <= 1) {
    std::cout << "There are no vertices in the graph" << std::endl;
    std::cout << "Please choose 1 and enter the name of the file: "
              << std::endl;
    return false;
  }
  return true;
}
void MainApplication::EnterFileNameToLoadGraph() {
  std::cout << "Enter the name of the file (type 'v' to load "
               "./examples/valid_matrix.txt): ";
  std::string filename;
  std::cin >> filename;
  if (!filename.empty()) {
    if (filename == "v") {
      filename = "./examples/valid_matrix.txt";
    }
    std::ifstream file(filename);
    if (file.good()) {
      g_.LoadGraphFromFile(filename);
      PrintGraph(g_);
      std::cout << "Graph loaded" << std::endl;
    } else {
      std::cout << "File not found" << std::endl;
    }
    file.close();
  } else {
    std::cout << "Please enter correct filename" << std::endl;
  }
}

void MainApplication::BreadthFirstSearch() {
  if (CheckGraphLoaded()) {
    auto result = ga_.BreadthFirstSearch(g_, 1);
    for (auto i : result) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
}

void MainApplication::DepthFirstSearch() {
  if (CheckGraphLoaded()) {
    auto result = ga_.DepthFirstSearch(g_, 1);
    for (auto i : result) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
}

void MainApplication::GetShortestPathBetweenVertices() {
  if (CheckGraphLoaded()) {
    std::cout << "Enter start and end vertices: x y" << std::endl;
    size_t x, y;
    std::cin >> x >> y;
    if (x > 0 && x <= g_.GetVertices() && y > 0 && y <= g_.GetVertices()) {
      auto result = ga_.GetShortestPathBetweenVertices(g_, x, y);
      std::cout << "Route length: " << result << std::endl;
    } else {
      std::cout << "Please enter correct vertices" << std::endl;
    }
  }
}

void MainApplication::GetShortestPathsBetweenAllVertices() {
  if (CheckGraphLoaded()) {
    auto result = ga_.GetShortestPathsBetweenAllVertices(g_);
    for (auto i : result) {
      for (auto j : i) {
        std::cout << j << " ";
      }
      std::cout << std::endl;
    }
  }
}

void MainApplication::GetLeastSpanningTree() {
  if (CheckGraphLoaded()) {
    auto result = ga_.GetLeastSpanningTree(g_);
    for (auto i : result) {
      for (auto j : i) {
        std::cout << j << " ";
      }
      std::cout << std::endl;
    }
  }
}

void MainApplication::SolveTravelingSalesmanProblem() {
  if (CheckGraphLoaded()) {
    TsmResult result = ga_.SolveTravelingSalesmanProblem(g_);
    PrintTsmResult(result);
  }
}

void MainApplication::ResearchTSPworkingTime() {
  if (CheckGraphLoaded()) {
    int countRun, i = 0;
    double average = 0;
    std::cout << "Enter the number of runs: ";
    std::cin >> countRun;
    auto start = std::chrono::steady_clock::now();
    while (i < countRun) {
      TsmResult result = ga_.SolveTravelingSalesmanProblem(g_);
      average += result.distance;
      ++i;
    }
    average /= countRun;
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> firstTSMresult = end - start;
    std::cout << "1) Ant Colony Optimization result: " << firstTSMresult.count()
              << " sec. (average = " << average << ")" << std::endl;
    i = 0;
    average = 0;
    start = std::chrono::steady_clock::now();
    while (i < countRun) {
      TsmResult result = ga_.SolveTravelingSalesmanProblemAnnealing(g_);
      average += result.distance;
      ++i;
    }
    average /= countRun;
    end = std::chrono::steady_clock::now();
    std::chrono::duration<double> secondTSMresult = end - start;
    std::cout << "2) Simulated annealing result: " << secondTSMresult.count()
              << " sec. (average = " << average << ")" << std::endl;
    i = 0;
    average = 0;
    start = std::chrono::steady_clock::now();
    while (i < countRun) {
      TsmResult result = ga_.SolveTravelingSalesmanProblemGenetic(g_);
      average += result.distance;
      ++i;
    }
    average /= countRun;
    end = std::chrono::steady_clock::now();
    std::chrono::duration<double> thirdTSMresult = end - start;
    std::cout << "3) Genetic algorithm result: " << thirdTSMresult.count()
              << " sec. (average = " << average << ")" << std::endl;
  }
}

void MainApplication::PrintGraph(Graph &graph) {
  std::vector<std::vector<uint32_t>> adj_matrix = graph.GetMatrix();
  uint32_t size = graph.GetVertices();
  std::cout << " ";
  std::cout << std::endl;
  for (uint32_t i = 0; i != size; ++i) {
    for (uint32_t j = 0; j != size; ++j) {
      std::cout << adj_matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void MainApplication::PrintTsmResult(TsmResult &result) {
  std::cout << "Sequence of traversing vertices: ";
  for (uint32_t i = 0; i != result.vertices.size() - 1; ++i) {
    std::cout << result.vertices[i] + 1 << "-";
  }
  std::cout << result.vertices[0] + 1 << std::endl;
  std::cout << "Route length: " << result.distance << std::endl;
}

}  // namespace s21
