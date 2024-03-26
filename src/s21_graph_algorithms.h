#ifndef A2_SIMPLENAVIGATOR_V1_0_1_S21_GRAPH_ALGORITHMS_H_
#define A2_SIMPLENAVIGATOR_V1_0_1_S21_GRAPH_ALGORITHMS_H_

#include <algorithm>

#include "helpers/s21_queue.h"
#include "helpers/s21_stack.h"
#include "s21_graph.h"

namespace s21 {

struct TsmResult {
  std::vector<uint32_t> vertices;
  double distance;
};

class GraphAlgorithms {
  using children = std::pair<std::vector<uint32_t>, std::vector<uint32_t>>;

 public:
  std::vector<uint32_t> DepthFirstSearch(Graph &graph, uint32_t start_vertex);
  std::vector<uint32_t> BreadthFirstSearch(Graph &graph, uint32_t start_vertex);
  uint32_t GetShortestPathBetweenVertices(Graph &graph, uint32_t vertex1,
                                          uint32_t vertex2);
  std::vector<std::vector<uint32_t>> GetShortestPathsBetweenAllVertices(
      Graph &graph);
  std::vector<std::vector<uint32_t>> GetLeastSpanningTree(Graph &graph);
  TsmResult SolveTravelingSalesmanProblem(Graph &graph);
  TsmResult SolveTravelingSalesmanProblemAnnealing(Graph &graph);
  TsmResult SolveTravelingSalesmanProblemGenetic(Graph &graph);

 private:
  // Const for SolveTravelingSalesmanProblem:
  double kDefaultPheromone = 1.0;
  int kAlpha = 1;
  int kBeta = 4;
  u_int32_t kIterationCount = 1000;
  double kEvaporationCoef = 0.4;
  double kUpdatePheromoneCoef = 240.0;

  // Const for SolveTravelingSalesmanProblemAnnealing:
  double kMaxTemperature = 10.0;
  double kMinTemperature = 0.00001;

  // Constant parameters for Genetic algorithm:
  uint32_t kPopulationSize = 1000;
  uint32_t kMutationProcent = 93;
  uint32_t kGenerations = 10000;
  uint32_t kDesiredRoute = 255;

  /// @brief Fills graph adjacency matrix with 0 (for main diagonal) and very
  /// big numbers (for absent edges)
  /// @param matrix graph adjacency matrix
  /// @param size amount og graph vertices
  /// @return adjacency matrix convenient for searching paths
  std::vector<std::vector<uint32_t>> PrepareGraphMatrix(
      std::vector<std::vector<uint32_t>> matrix, uint32_t size);

  // Helper methods for SolveTravelingSalesmanProblem:
  bool NextAntStep(TsmResult &ant, Graph &graph,
                   std::vector<std::vector<double>> &pheromone,
                   std::vector<bool> &visited);
  void UpdatePheromoneLocal(TsmResult &ant,
                            std::vector<std::vector<double>> &pheromone_delta);
  void UpdatePheromoneGlobal(std::vector<std::vector<double>> &pheromone,
                             std::vector<std::vector<double>> &pheromone_delta);
  bool CheckVisited(std::vector<bool> &visited);

  // Helper methods for SolveTravelingSalesmanProblemAnnealing:
  TsmResult GetRandomPath(Graph &graph);
  TsmResult SwapRandomElements(TsmResult path, Graph &graph);

  // Helper methods for Genetic algorithm:
  /// @brief Gets fitness parameter (in TSMP case - route length)
  /// @param solution a set of genes (a route)
  /// @return route length
  uint32_t FitnessFunction(std::vector<uint32_t> solution);
  std::vector<std::vector<uint32_t>> CreateInitialPopulation(
      uint32_t k_population_size, uint32_t graph_size);
  uint32_t GenerateRandomNumber(uint32_t size);
  std::vector<uint32_t> GenerateRandomGenotype(uint32_t size);
  std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
  SortPopulation(uint32_t size, std::vector<std::vector<uint32_t>> population,
                 std::vector<std::vector<uint32_t>> matrix);
  std::vector<uint32_t> GetPhenotype(uint32_t size,
                                     std::vector<std::vector<uint32_t>> matrix,
                                     std::vector<uint32_t> genotype);
  void GeneticAlgorithm(
      uint32_t size, std::vector<std::vector<uint32_t>> &population,
      std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
          &sorted_population,
      std::vector<std::vector<uint32_t>> matrix);

  /// @brief Creates 2 new solutions (2 new routes)
  /// @param size amount of graph vertices
  /// @param population array of routes
  /// @return 2 new sets of genes (2 new routes)
  children Crossover(uint32_t size,
                     std::vector<std::vector<uint32_t>> &population);
  std::vector<uint32_t> CreateOffspringSequence(
      uint32_t size, std::pair<uint32_t, uint32_t> parents,
      uint32_t crossover_point, std::vector<std::vector<uint32_t>> &population);
  uint32_t FillAfterCrossover(uint32_t crossover_point, uint32_t size,
                              std::vector<uint32_t> &inherited_genes,
                              std::vector<std::vector<uint32_t>> &population,
                              uint32_t parent_num,
                              std::vector<uint32_t> &offspring_seq, uint32_t &j,
                              bool parent_2);

  /// @brief Changes 2 edjes of a route randomly
  /// @param solution_genotype set of genes (route)
  /// @param size amount of graph vertices
  void Mutation(std::vector<uint32_t> &solution_genotype, uint32_t size);

  /// @brief Adds 2 new sets of genes (routes) to the population
  /// @param population array of routes
  /// @param new_solutions 2 new sets of genes (routes)
  /// @param sorted_population a structure which contains genotypes (routes),
  /// phenotypes (distances between vertices of this route) and fitness function
  /// (total route length) for each route
  /// @param matrix graph adjacency matrix
  void UpdatePopulation(
      children new_solutions,
      std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
          &sorted_population,
      std::vector<std::vector<uint32_t>> matrix);
  void AddSolutionToPopulation(
      std::vector<uint32_t> solution,
      std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
          &sorted_population,
      std::vector<std::vector<uint32_t>> matrix);

  /// @brief Removes the longest routes from the population
  /// @param sorted_population a structure which contains genotypes (routes),
  /// phenotypes (distances between vertices of this route) and fitness function
  /// (total route length) for each route
  void Selection(
      std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
          &sorted_population);
  void MinRoute(
      std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
          &sorted_population,
      TsmResult &result);

  // Other helper fields for Genetic algorithm:
  uint32_t rand_ = time(nullptr);
};

}  // namespace s21

#endif  // A2_SIMPLENAVIGATOR_V1_0_1_S21_GRAPH_ALGORITHMS_H_
