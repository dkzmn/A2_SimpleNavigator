#include "s21_graph_algorithms.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace s21 {
std::vector<uint32_t> GraphAlgorithms::DepthFirstSearch(Graph &graph,
                                                        uint32_t start_vertex) {
  int size = graph.GetVertices();
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  std::vector<uint32_t> result;
  std::vector<bool> visited(size, false);
  s21::stack<uint32_t> temp;
  temp.push(start_vertex - 1);
  visited[start_vertex - 1] = true;
  while (!temp.empty()) {
    uint32_t current_node = temp.top();
    temp.pop();
    result.push_back(current_node + 1);
    for (int i = 0; i != size; ++i) {
      if ((graph_matrix[current_node][i] > 0) && (!visited[i])) {
        visited[i] = true;
        temp.push(i);
      }
    }
  }
  return result;
}

std::vector<uint32_t> GraphAlgorithms::BreadthFirstSearch(
    Graph &graph, uint32_t start_vertex) {
  int size = graph.GetVertices();
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  std::vector<uint32_t> result;
  std::vector<bool> visited(size, false);
  s21::queue<uint32_t> temp;
  temp.push(start_vertex - 1);
  visited[start_vertex - 1] = true;
  while (!temp.empty()) {
    uint32_t current_node = temp.front();
    temp.pop();
    result.push_back(current_node + 1);
    for (int i = 0; i != size; ++i) {
      if ((graph_matrix[current_node][i] > 0) && (!visited[i])) {
        visited[i] = true;
        temp.push(i);
      }
    }
  }
  return result;
}

uint32_t GraphAlgorithms::GetShortestPathBetweenVertices(Graph &graph,
                                                         uint32_t vertex1,
                                                         uint32_t vertex2) {
  if (!graph.GraphIsWeighted())
    throw std::logic_error("The graph is not weighted!");
  int size = graph.GetVertices();
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  uint32_t max_uint = std::numeric_limits<uint32_t>::max(), min = 0,
           min_ind = 0;
  std::vector<uint32_t> dist(size, max_uint);
  std::vector<bool> visited(size, false);
  dist[vertex1 - 1] = 0;
  do {
    min = max_uint;
    min_ind = max_uint;
    for (int i = 0; i != size; ++i) {
      if ((!visited[i]) && (dist[i] < min)) {
        min = dist[i];
        min_ind = i;
      }
    }
    if (min_ind != max_uint) {
      for (int i = 0; i != size; ++i) {
        if (graph_matrix[min_ind][i] > 0) {
          uint32_t temp = min + graph_matrix[min_ind][i];
          if (temp < dist[i]) dist[i] = temp;
        }
      }
      visited[min_ind] = true;
    }
  } while (min_ind < max_uint);
  return dist[vertex2 - 1];
}

std::vector<std::vector<uint32_t>>
GraphAlgorithms::GetShortestPathsBetweenAllVertices(Graph &graph) {
  if (!graph.GraphIsWeighted())
    throw std::logic_error("The graph is not weighted!");
  uint32_t size = graph.GetVertices();
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  graph_matrix = PrepareGraphMatrix(graph_matrix, size);
  for (uint32_t k = 0; k != size; ++k) {
    for (uint32_t i = 0; i != size; ++i) {
      for (uint32_t j = 0; j != size; ++j) {
        uint32_t new_sum = 0;
        if ((graph_matrix[i][k] == std::numeric_limits<uint32_t>::max()) ||
            (graph_matrix[k][j] == std::numeric_limits<uint32_t>::max()))
          new_sum = std::numeric_limits<uint32_t>::max();
        else
          new_sum = graph_matrix[i][k] + graph_matrix[k][j];
        graph_matrix[i][j] = std::min(graph_matrix[i][j], new_sum);
      }
    }
  }
  return graph_matrix;
}

std::vector<std::vector<uint32_t>> GraphAlgorithms::GetLeastSpanningTree(
    Graph &graph) {
  if (!graph.GraphIsWeighted())
    throw std::logic_error("The graph is not weighted!");
  uint32_t size = graph.GetVertices();
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  graph_matrix = PrepareGraphMatrix(graph_matrix, size);
  std::vector<std::vector<uint32_t>> mst_matrix(size,
                                                std::vector<uint32_t>(size));
  std::vector<bool> selected(size, false);
  uint32_t max_uint = std::numeric_limits<uint32_t>::max(), min = 0,
           min_ind_row = 0, min_ind_col = 0, edges = 0;
  selected[0] = true;
  while (edges != size - 1) {
    min = max_uint;
    min_ind_row = 0;
    min_ind_col = 0;
    for (uint32_t i = 0; i != size; ++i) {
      if (selected[i]) {
        for (uint32_t j = 0; j != size; ++j) {
          if ((graph_matrix[i][j] > 0) && (!selected[j]) &&
              (graph_matrix[i][j] < min)) {
            min = graph_matrix[i][j];
            min_ind_row = i;
            min_ind_col = j;
          }
        }
      }
    }
    mst_matrix[min_ind_row][min_ind_col] = min;
    selected[min_ind_col] = true;
    ++edges;
  }
  if (!graph.GraphIsDirected()) {
    for (uint32_t i = 0; i != size; ++i) {
      for (uint32_t j = 0; j != size; ++j) {
        if (mst_matrix[i][j]) mst_matrix[j][i] = mst_matrix[i][j];
      }
    }
  }
  return mst_matrix;
}

std::vector<std::vector<uint32_t>> GraphAlgorithms::PrepareGraphMatrix(
    std::vector<std::vector<uint32_t>> matrix, uint32_t size) {
  for (uint32_t i = 0; i != size; ++i) {
    matrix[i][i] = 0;
    for (uint32_t j = 0; j != size; ++j) {
      if ((i != j) && (matrix[i][j] == 0))
        matrix[i][j] = std::numeric_limits<uint32_t>::max();
    }
  }
  return matrix;
}

TsmResult GraphAlgorithms::SolveTravelingSalesmanProblem(Graph &graph) {
  if (!graph.GraphIsFull()) throw std::logic_error("The graph is not full");
  uint32_t vert_count = graph.GetVertices();
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  TsmResult result;
  result.distance = UINT32_MAX;
  std::vector<std::vector<double>> pheromone(
      vert_count, std::vector<double>(vert_count, kDefaultPheromone));
  std::vector<std::vector<double>> pheromone_delta(
      vert_count, std::vector<double>(vert_count, 0.0));
  for (uint32_t iteration = 0; iteration < kIterationCount; ++iteration) {
    for (uint32_t i = 0; i < vert_count; ++i) {
      TsmResult ant{};
      ant.vertices.push_back(i);
      std::vector<bool> visited(vert_count, false);
      bool avaliable = true;
      while (avaliable) avaliable = NextAntStep(ant, graph, pheromone, visited);
      if (graph_matrix[ant.vertices.back()][i] > 0 && CheckVisited(visited)) {
        ant.distance += graph_matrix[ant.vertices.back()][i];
        ant.vertices.push_back(i);
        UpdatePheromoneLocal(ant, pheromone_delta);
        if (ant.distance < result.distance) {
          result = std::move(ant);
        }
      }
    }
    UpdatePheromoneGlobal(pheromone, pheromone_delta);
  }
  return result;
}

bool GraphAlgorithms::NextAntStep(TsmResult &ant, Graph &graph,
                                  std::vector<std::vector<double>> &pheromone,
                                  std::vector<bool> &visited) {
  uint32_t current_vert = ant.vertices.back();
  visited[current_vert] = true;
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  uint32_t vert_count = graph.GetVertices();
  std::vector<double> probability(vert_count, 0.0);
  double probability_sum = 0.0;
  bool found = false;
  for (uint32_t i = 0; i < vert_count; ++i) {
    if (graph_matrix[current_vert][i] != 0 && visited[i] == false) {
      probability[i] = std::pow(pheromone[current_vert][i], kAlpha) /
                       std::pow(graph_matrix[current_vert][i], kBeta);
      probability_sum += probability[i];
      found = true;
    }
  }
  if (found) {
    for (uint32_t i = 0; i < vert_count; ++i) probability[i] /= probability_sum;
    double random = (double)rand() / RAND_MAX, tmp_sum = 0.0;
    uint32_t select_idx = 0;
    while (tmp_sum + probability[select_idx] <= random)
      tmp_sum += probability[select_idx++];
    ant.vertices.push_back(select_idx);
    ant.distance += graph_matrix[current_vert][select_idx];
  }
  return found;
}

void GraphAlgorithms::UpdatePheromoneLocal(
    TsmResult &ant, std::vector<std::vector<double>> &pheromone_delta) {
  for (uint32_t i = 1; i < ant.vertices.size(); ++i) {
    pheromone_delta[ant.vertices[i - 1]][ant.vertices[i]] +=
        kUpdatePheromoneCoef / ant.distance;
  }
}

void GraphAlgorithms::UpdatePheromoneGlobal(
    std::vector<std::vector<double>> &pheromone,
    std::vector<std::vector<double>> &pheromone_delta) {
  for (uint32_t i = 0; i < pheromone.size(); ++i) {
    for (uint32_t j = 0; j < pheromone[i].size(); ++j) {
      pheromone[i][j] =
          (1 - kEvaporationCoef) * pheromone[i][j] + pheromone_delta[i][j];
      if (pheromone[i][j] < 0.1) pheromone[i][j] = 0.1;
    }
  }
}

bool GraphAlgorithms::CheckVisited(std::vector<bool> &visited) {
  return std::find(visited.begin(), visited.end(), false) == visited.end();
}

TsmResult GraphAlgorithms::SolveTravelingSalesmanProblemAnnealing(
    Graph &graph) {
  if (!graph.GraphIsFull()) throw std::logic_error("The graph is not full");
  double temperature = kMaxTemperature;
  uint32_t iteration = 1;
  TsmResult current_path = GetRandomPath(graph);
  while (temperature > kMinTemperature) {
    TsmResult new_path = SwapRandomElements(current_path, graph);
    if (new_path.distance <= current_path.distance) {
      current_path = new_path;
    } else {
      double random = (double)rand() / RAND_MAX;
      double probability =
          exp((current_path.distance - new_path.distance) / temperature);
      if (random <= probability) {
        current_path = new_path;
      }
    }
    temperature = kMaxTemperature * 0.1 / iteration;
    ++iteration;
  }
  return current_path;
}

TsmResult GraphAlgorithms::GetRandomPath(Graph &graph) {
  TsmResult res{};
  uint32_t vert_count = graph.GetVertices();
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  for (uint32_t i = 0; i < vert_count; ++i) res.vertices.push_back(i);
  std::shuffle(res.vertices.begin(), res.vertices.end(), std::random_device());
  res.vertices.push_back(res.vertices[0]);
  for (uint32_t i = 1; i < res.vertices.size(); ++i) {
    res.distance += graph_matrix[res.vertices[i - 1]][res.vertices[i]];
  }
  return res;
}

TsmResult GraphAlgorithms::SwapRandomElements(TsmResult path, Graph &graph) {
  if (path.vertices.size() < 3) throw std::logic_error("The path is too short");
  TsmResult res{};
  res.vertices = path.vertices;
  res.vertices.pop_back();
  uint32_t v1 = rand() % res.vertices.size();
  uint32_t v2 = v1;
  while (v2 == v1) {
    v2 = rand() % res.vertices.size();
  }
  std::swap(res.vertices[v1], res.vertices[v2]);
  res.vertices.push_back(res.vertices[0]);
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  for (uint32_t i = 1; i < res.vertices.size(); ++i) {
    res.distance += graph_matrix[res.vertices[i - 1]][res.vertices[i]];
  }
  return res;
}

TsmResult GraphAlgorithms::SolveTravelingSalesmanProblemGenetic(Graph &graph) {
  if (!graph.GraphIsFull()) throw std::logic_error("The graph is not full");
  uint32_t size = graph.GetVertices();
  std::vector<std::vector<uint32_t>> graph_matrix = graph.GetMatrix();
  std::vector<std::vector<uint32_t>> population =
      CreateInitialPopulation(kPopulationSize, size);
  uint32_t generation = 1;
  std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
      sorted_population = SortPopulation(size + 1, population, graph_matrix);
  TsmResult result = {std::vector<uint32_t>(sorted_population.size()),
                      std::numeric_limits<uint32_t>::max()};
  while (generation < kGenerations) {
    GeneticAlgorithm(size + 1, population, sorted_population, graph_matrix);
    ++generation;
    MinRoute(sorted_population, result);
    if (result.distance <= kDesiredRoute) break;
  }
  return result;
}

uint32_t GraphAlgorithms::FitnessFunction(std::vector<uint32_t> solution) {
  uint32_t distance = 0;
  for (size_t i = 0; i != solution.size(); ++i) {
    distance += solution[i];
  }
  return distance;
}

void GraphAlgorithms::MinRoute(
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &sorted_population,
    TsmResult &result) {
  uint32_t ind = 0;
  for (uint32_t i = 0; i != sorted_population.size(); ++i) {
    if (sorted_population[i].second[sorted_population[i].second.size() - 1] <
        result.distance) {
      result.distance =
          sorted_population[i].second[sorted_population[i].second.size() - 1];
      ind = i;
    }
  }
  result.vertices = sorted_population[ind].first;
}

std::vector<std::vector<uint32_t>> GraphAlgorithms::CreateInitialPopulation(
    uint32_t k_population_size, uint32_t graph_size) {
  std::vector<std::vector<uint32_t>> population =
      std::vector<std::vector<uint32_t>>(
          k_population_size, std::vector<uint32_t>(k_population_size));
  // Create first solution;
  for (uint32_t i = 0; i != graph_size; ++i) {
    population[0][i] = i;
  }
  population[0][graph_size] = 0;
  // Generate other solutions;
  for (uint32_t i = 1; i != k_population_size; ++i) {
    population[i] = GenerateRandomGenotype(graph_size);
  }
  return population;
}

uint32_t GraphAlgorithms::GenerateRandomNumber(uint32_t size) {
  srand((unsigned)rand_);
  rand_ = rand();
  return rand_ % size;
}

std::vector<uint32_t> GraphAlgorithms::GenerateRandomGenotype(uint32_t size) {
  std::vector<uint32_t> random_seq = std::vector<uint32_t>(size + 1, 0);
  std::vector<bool> existing = std::vector<bool>(size, false);
  existing[0] = true;
  for (uint32_t i = 1; i != size; ++i) {
    while (random_seq[i] == 0) {
      uint32_t random_num = GenerateRandomNumber(size);
      if (!existing[random_num]) random_seq[i] = random_num;
      existing[random_num] = true;
    }
  }
  random_seq[size] = 0;
  return random_seq;
}

std::vector<uint32_t> GraphAlgorithms::GetPhenotype(
    uint32_t size, std::vector<std::vector<uint32_t>> matrix,
    std::vector<uint32_t> genotype) {
  std::vector<uint32_t> phenotype = std::vector<uint32_t>(size);
  for (uint32_t i = 0; i != size - 1; ++i) {
    phenotype[i] = matrix[genotype[i + 1]][genotype[i]];
  }
  phenotype[size - 1] = matrix[genotype[0]][genotype[genotype.size() - 1]];
  return phenotype;
}

std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
GraphAlgorithms::SortPopulation(uint32_t size,
                                std::vector<std::vector<uint32_t>> population,
                                std::vector<std::vector<uint32_t>> matrix) {
  std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
      population_to_sort;
  for (size_t i = 0; i != population.size(); ++i) {
    std::vector<uint32_t> genotype = std::vector<uint32_t>(size);
    genotype = population[i];
    std::vector<uint32_t> phenotype = GetPhenotype(size, matrix, genotype);
    phenotype.push_back(FitnessFunction(phenotype));
    population_to_sort.push_back({genotype, phenotype});
  }
  std::sort(population_to_sort.begin(), population_to_sort.end(),
            [size](const auto &el1, const auto &el2) {
              return el1.second[size] < el2.second[size];
            });
  return population_to_sort;
}

void GraphAlgorithms::GeneticAlgorithm(
    uint32_t size, std::vector<std::vector<uint32_t>> &population,
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &sorted_population,
    std::vector<std::vector<uint32_t>> matrix) {
  children new_solutions = Crossover(size, population);
  Mutation(new_solutions.first, size);
  Mutation(new_solutions.second, size);
  UpdatePopulation(new_solutions, sorted_population, matrix);
  Selection(sorted_population);
}

GraphAlgorithms::children GraphAlgorithms::Crossover(
    uint32_t size, std::vector<std::vector<uint32_t>> &population) {
  std::pair<uint32_t, uint32_t> parents;
  parents.first = GenerateRandomNumber(kPopulationSize);
  parents.second = GenerateRandomNumber(kPopulationSize);
  if (parents.first == parents.second) {
    do {
      parents.second = GenerateRandomNumber(kPopulationSize);
    } while (parents.first == parents.second);
  }
  uint32_t crossover_point = GenerateRandomNumber(size - 1) + 1;
  children new_solutions;
  new_solutions.first =
      CreateOffspringSequence(size, parents, crossover_point, population);
  new_solutions.second = CreateOffspringSequence(
      size, {parents.second, parents.first}, crossover_point, population);
  return new_solutions;
}

std::vector<uint32_t> GraphAlgorithms::CreateOffspringSequence(
    uint32_t size, std::pair<uint32_t, uint32_t> parents,
    uint32_t crossover_point, std::vector<std::vector<uint32_t>> &population) {
  std::vector<uint32_t> inherited_genes;
  inherited_genes.push_back(0);
  std::vector<uint32_t> offspring_seq = std::vector<uint32_t>(size, 0);
  for (uint32_t i = 1; i != crossover_point; ++i) {
    offspring_seq[i] = population[parents.first][i];
    inherited_genes.push_back(offspring_seq[i]);
  }
  uint32_t parent_num = parents.second, j = crossover_point;
  uint32_t parent_2_genes =
      FillAfterCrossover(crossover_point, size, inherited_genes, population,
                         parent_num, offspring_seq, j, true);
  if ((parent_2_genes == 0) && (j < size)) {
    parent_num = parents.first;
    FillAfterCrossover(crossover_point, size, inherited_genes, population,
                       parent_num, offspring_seq, j, false);
  }
  return offspring_seq;
}

uint32_t GraphAlgorithms::FillAfterCrossover(
    uint32_t crossover_point, uint32_t size,
    std::vector<uint32_t> &inherited_genes,
    std::vector<std::vector<uint32_t>> &population, uint32_t parent_num,
    std::vector<uint32_t> &offspring_seq, uint32_t &j, bool parent_2) {
  uint32_t parent_2_genes = (parent_2) ? (size - crossover_point) : 0;
  for (uint32_t i = crossover_point; i != size; ++i) {
    if (std::find(inherited_genes.begin(), inherited_genes.end(),
                  population[parent_num][i]) == inherited_genes.end()) {
      offspring_seq[j] = population[parent_num][i];
      inherited_genes.push_back(offspring_seq[j]);
      ++j;
    }
    if (parent_2) --parent_2_genes;
  }
  return parent_2_genes;
}

void GraphAlgorithms::Mutation(std::vector<uint32_t> &solution_genotype,
                               uint32_t size) {
  uint32_t mutation = GenerateRandomNumber(100);
  if (mutation <= kMutationProcent) {
    uint32_t random_chromosome_1 = GenerateRandomNumber(size);
    uint32_t random_chromosome_2 = GenerateRandomNumber(size);
    if (random_chromosome_1 == random_chromosome_2) {
      do {
        random_chromosome_2 = GenerateRandomNumber(size);
      } while (random_chromosome_1 == random_chromosome_2);
    }
    uint32_t temp = solution_genotype[random_chromosome_1];
    solution_genotype[random_chromosome_1] =
        solution_genotype[random_chromosome_2];
    solution_genotype[random_chromosome_2] = temp;
  }
}

void GraphAlgorithms::UpdatePopulation(
    children new_solutions,
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &sorted_population,
    std::vector<std::vector<uint32_t>> matrix) {
  AddSolutionToPopulation(new_solutions.first, sorted_population, matrix);
  AddSolutionToPopulation(new_solutions.second, sorted_population, matrix);
}

void GraphAlgorithms::AddSolutionToPopulation(
    std::vector<uint32_t> solution_genotype,
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &sorted_population,
    std::vector<std::vector<uint32_t>> matrix) {
  std::vector<uint32_t> solution_phenotype =
      GetPhenotype(solution_genotype.size(), matrix, solution_genotype);
  solution_phenotype.push_back(FitnessFunction(solution_phenotype));
  std::pair<std::vector<uint32_t>, std::vector<uint32_t>> solution = {
      solution_genotype, solution_phenotype};
  auto pos = sorted_population.begin();
  uint32_t ind = 0;
  for (auto i = sorted_population.begin(); i != sorted_population.end(); ++i) {
    if (sorted_population[ind]
            .second[sorted_population[ind].second.size() - 1] >
        solution_phenotype[solution_phenotype.size() - 1]) {
      pos = i;
      break;
    }
    ++ind;
  }
  sorted_population.insert(pos, solution);
}

void GraphAlgorithms::Selection(
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &sorted_population) {
  while (sorted_population.size() > kPopulationSize) {
    sorted_population.erase(sorted_population.end() - 1);
  }
}

}  // namespace s21
