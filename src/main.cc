#include "s21_main_application.h"

int main() {
  s21::MainApplication m;
  int choice = 0;
  while (choice >= 0) {
    m.PrintMenu();
    choice = m.GetChoice();
    switch (choice) {
      case 0:
        std::cout << "Exit" << std::endl;
        exit(0);
        break;
      case 1:
        m.EnterFileNameToLoadGraph();
        break;
      case 2:
        m.BreadthFirstSearch();
        break;
      case 3:
        m.DepthFirstSearch();
        break;
      case 4:
        m.GetShortestPathBetweenVertices();
        break;
      case 5:
        m.GetShortestPathsBetweenAllVertices();
        break;
      case 6:
        m.GetLeastSpanningTree();
        break;
      case 7:
        m.SolveTravelingSalesmanProblem();
        break;
      case 8:
        m.ResearchTSPworkingTime();
        break;
      default:
        std::cout << "Wrong choice" << std::endl;
        break;
    }
  }
  return 0;
}
