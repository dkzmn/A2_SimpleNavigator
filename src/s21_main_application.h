#ifndef A2_SIMPLENAVIGATOR_V1_0_1_S21_MAIN_APPLICATION_H_
#define A2_SIMPLENAVIGATOR_V1_0_1_S21_MAIN_APPLICATION_H_

#include "s21_graph.h"
#include "s21_graph_algorithms.h"

namespace s21 {
class MainApplication {
 public:
  void PrintMenu();
  int GetChoice();
  void EnterFileNameToLoadGraph();
  void BreadthFirstSearch();
  void DepthFirstSearch();
  void GetShortestPathBetweenVertices();
  void GetShortestPathsBetweenAllVertices();
  void GetLeastSpanningTree();
  void SolveTravelingSalesmanProblem();
  bool CheckGraphLoaded();
  void ResearchTSPworkingTime();
  // Printing methods:
  void PrintGraph(Graph &graph);
  void PrintTsmResult(TsmResult &result);

 private:
  Graph g_;
  GraphAlgorithms ga_;
};
}  // namespace s21

#endif  // A2_SIMPLENAVIGATOR_V1_0_1_S21_MAIN_APPLICATION_H_
