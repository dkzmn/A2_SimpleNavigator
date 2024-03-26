#include "tests.h"

TEST(Graph_Algorithm, DFS) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<uint32_t> result = graph_algorithm.DepthFirstSearch(test, 1);
  std::vector<uint32_t> answer;
  answer.push_back(1);
  answer.push_back(3);
  answer.push_back(6);
  answer.push_back(7);
  answer.push_back(5);
  answer.push_back(4);
  answer.push_back(2);
  bool exit = false;
  while (!exit) {
    EXPECT_EQ(result.back(), answer.back());
    result.pop_back();
    answer.pop_back();
    if (result.empty()) exit = true;
  }
}

TEST(Graph_Algorithm, BFS) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<uint32_t> result = graph_algorithm.BreadthFirstSearch(test, 1);
  std::vector<uint32_t> answer;
  answer.push_back(1);
  answer.push_back(2);
  answer.push_back(3);
  answer.push_back(4);
  answer.push_back(5);
  answer.push_back(6);
  answer.push_back(7);
  bool exit = false;
  while (!exit) {
    EXPECT_EQ(result.back(), answer.back());
    result.pop_back();
    answer.pop_back();
    if (result.empty()) exit = true;
  }
}

TEST(Graph_Algorithm, Shortest_Path) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph_weighted.txt");
  s21::GraphAlgorithms graph_algorithm;
  uint32_t result = graph_algorithm.GetShortestPathBetweenVertices(test, 1, 5);
  EXPECT_EQ(result, static_cast<uint32_t>(20));
  uint32_t result2 = graph_algorithm.GetShortestPathBetweenVertices(test, 1, 4);
  EXPECT_EQ(result2, static_cast<uint32_t>(20));
  uint32_t result3 = graph_algorithm.GetShortestPathBetweenVertices(test, 1, 6);
  EXPECT_EQ(result3, static_cast<uint32_t>(11));
}

TEST(Graph_Algorithm, Shortest_Path_2) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/min_spanning_tree_2.txt");
  s21::GraphAlgorithms graph_algorithm;
  uint32_t result = graph_algorithm.GetShortestPathBetweenVertices(test, 1, 9);
  EXPECT_EQ(result, static_cast<uint32_t>(10));
  uint32_t result2 = graph_algorithm.GetShortestPathBetweenVertices(test, 7, 9);
  EXPECT_EQ(result2, static_cast<uint32_t>(13));
  uint32_t result3 = graph_algorithm.GetShortestPathBetweenVertices(test, 3, 8);
  EXPECT_EQ(result3, static_cast<uint32_t>(17));
}

TEST(Graph_Algorithm, Shortest_Path_Err) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph.txt");
  s21::GraphAlgorithms graph_algorithm;
  EXPECT_THROW(graph_algorithm.GetShortestPathBetweenVertices(test, 1, 9), std::logic_error);
}

TEST(Graph_Algorithm, Shortest_Path_All_Verices) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph_weighted.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<std::vector<uint32_t>> result =
      graph_algorithm.GetShortestPathsBetweenAllVertices(test);
  s21::Graph answer;
  answer.LoadGraphFromFile("examples/shortest_path_between_all.txt");
  std::vector<std::vector<uint32_t>> answer_matrix = answer.GetMatrix();
  for (uint32_t i = 0; i < test.GetVertices(); ++i) {
    for (uint32_t j = 0; j < test.GetVertices(); ++j) {
      EXPECT_EQ(answer_matrix[i][j], result[i][j]);
    }
  }
}

TEST(Graph_Algorithm, Shortest_Path_All_Verices_2) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/min_spanning_tree_2.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<std::vector<uint32_t>> result =
      graph_algorithm.GetShortestPathsBetweenAllVertices(test);
  s21::Graph answer;
  answer.LoadGraphFromFile("examples/mst_2_shortest_all.txt");
  std::vector<std::vector<uint32_t>> answer_matrix = answer.GetMatrix();
  for (uint32_t i = 0; i < test.GetVertices(); ++i) {
    for (uint32_t j = 0; j < test.GetVertices(); ++j) {
      EXPECT_EQ(answer_matrix[i][j], result[i][j]);
    }
  }
}

TEST(Graph_Algorithm, Shortest_Path_All_Verices_3) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/min_spanning_tree.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<std::vector<uint32_t>> result =
      graph_algorithm.GetShortestPathsBetweenAllVertices(test);
  s21::Graph answer;
  answer.LoadGraphFromFile("examples/mst_shortest_all.txt");
  std::vector<std::vector<uint32_t>> answer_matrix = answer.GetMatrix();
  for (uint32_t i = 0; i < test.GetVertices(); ++i) {
    for (uint32_t j = 0; j < test.GetVertices(); ++j) {
      EXPECT_EQ(answer_matrix[i][j], result[i][j]);
    }
  }
}

TEST(Graph_Algorithm, Shortest_Path_All_Verices_Err) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph.txt");
  s21::GraphAlgorithms graph_algorithm;
  EXPECT_THROW(graph_algorithm.GetShortestPathsBetweenAllVertices(test), std::logic_error);
}

TEST(Graph_Algorithm, Min_Spanning_Tree) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/min_spanning_tree.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<std::vector<uint32_t>> result =
      graph_algorithm.GetLeastSpanningTree(test);
  s21::Graph answer;
  answer.LoadGraphFromFile("examples/mst_result.txt");
  std::vector<std::vector<uint32_t>> answer_matrix = answer.GetMatrix();
  for (uint32_t i = 0; i < test.GetVertices(); ++i) {
    for (uint32_t j = 0; j < test.GetVertices(); ++j) {
      EXPECT_EQ(answer_matrix[i][j], result[i][j]);
    }
  }
}

TEST(Graph_Algorithm, Min_Spanning_Tree_Directed) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph_wd.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<std::vector<uint32_t>> result =
      graph_algorithm.GetLeastSpanningTree(test);
  s21::Graph answer;
  answer.LoadGraphFromFile("examples/small_graph_wd_result.txt");
  std::vector<std::vector<uint32_t>> answer_matrix = answer.GetMatrix();
  for (uint32_t i = 0; i < test.GetVertices(); ++i) {
    for (uint32_t j = 0; j < test.GetVertices(); ++j) {
      EXPECT_EQ(answer_matrix[i][j], result[i][j]);
    }
  }
}

TEST(Graph_Algorithm, Min_Spanning_Tree_2) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph_weighted.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<std::vector<uint32_t>> result =
      graph_algorithm.GetLeastSpanningTree(test);
  s21::Graph answer;
  answer.LoadGraphFromFile("examples/small_graph_weighted_mst.txt");
  std::vector<std::vector<uint32_t>> answer_matrix = answer.GetMatrix();
  for (uint32_t i = 0; i < test.GetVertices(); ++i) {
    for (uint32_t j = 0; j < test.GetVertices(); ++j) {
      EXPECT_EQ(answer_matrix[i][j], result[i][j]);
    }
  }
}

TEST(Graph_Algorithm, Min_Spanning_Tree_3) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/min_spanning_tree_2.txt");
  s21::GraphAlgorithms graph_algorithm;
  std::vector<std::vector<uint32_t>> result =
      graph_algorithm.GetLeastSpanningTree(test);
  s21::Graph answer;
  answer.LoadGraphFromFile("examples/mst_2_result.txt");
  std::vector<std::vector<uint32_t>> answer_matrix = answer.GetMatrix();
  for (uint32_t i = 0; i < test.GetVertices(); ++i) {
    for (uint32_t j = 0; j < test.GetVertices(); ++j) {
      EXPECT_EQ(answer_matrix[i][j], result[i][j]);
    }
  }
}

TEST(Graph_Algorithm, Min_Spanning_Tree_Err) {
  s21::Graph test;
  test.LoadGraphFromFile("examples/small_graph.txt");
  s21::GraphAlgorithms graph_algorithm;
  EXPECT_THROW(graph_algorithm.GetLeastSpanningTree(test), std::logic_error);
}

TEST(Graph_Algorithm, TSP_Ant) {
  s21::Graph g;
  g.LoadGraphFromFile("examples/valid_matrix.txt");
  s21::GraphAlgorithms ga;
  s21::TsmResult r = ga.SolveTravelingSalesmanProblem(g);
  EXPECT_EQ(r.vertices.back(), r.vertices.front());
  EXPECT_LT(r.distance, 300.0);
}

TEST(Graph_Algorithm, TSP_Annealing) {
  s21::Graph g;
  g.LoadGraphFromFile("examples/valid_matrix.txt");
  s21::GraphAlgorithms ga;
  s21::TsmResult r = ga.SolveTravelingSalesmanProblemAnnealing(g);
  EXPECT_EQ(r.vertices.back(), r.vertices.front());
  EXPECT_LT(r.distance, 300.0);
}

TEST(Graph_Algorithm, TSP_Genetic) {
  s21::Graph g;
  g.LoadGraphFromFile("examples/valid_matrix.txt");
  s21::GraphAlgorithms ga;
  s21::TsmResult r = ga.SolveTravelingSalesmanProblemGenetic(g);
  EXPECT_EQ(r.vertices.back(), r.vertices.front());
  EXPECT_LT(r.distance, 300.0);
}