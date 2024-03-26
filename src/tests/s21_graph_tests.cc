#include "tests.h"

TEST(Graph, Constructor_Default) {
  s21::Graph graph;
  EXPECT_EQ(graph.GetVertices(), static_cast<uint32_t>(1));
}

TEST(Graph, Constructor_Parametrized) {
  s21::Graph graph(5);
  EXPECT_EQ(graph.GetVertices(), static_cast<uint32_t>(5));
}

TEST(Graph, Constructor_Parametrized_Wrong) {
  EXPECT_THROW(s21::Graph graph(0), std::out_of_range);
}

TEST(Graph, Load_From_File) {
  s21::Graph graph;
  graph.LoadGraphFromFile("examples/valid_matrix.txt");
  EXPECT_EQ(graph.GetVertices(), static_cast<uint32_t>(11));
  std::vector<std::vector<uint32_t>> answer_matrix(
      graph.GetVertices(), std::vector<uint32_t>(graph.GetVertices()));
  answer_matrix[0][0] = 0;
  answer_matrix[0][1] = 29;
  answer_matrix[0][2] = 20;
  answer_matrix[0][3] = 21;
  answer_matrix[0][4] = 16;
  answer_matrix[0][5] = 31;
  answer_matrix[0][6] = 100;
  answer_matrix[0][7] = 12;
  answer_matrix[0][8] = 4;
  answer_matrix[0][9] = 31;
  answer_matrix[0][10] = 18;
  answer_matrix[1][0] = 29;
  answer_matrix[1][1] = 0;
  answer_matrix[1][2] = 15;
  answer_matrix[1][3] = 29;
  answer_matrix[1][4] = 28;
  answer_matrix[1][5] = 40;
  answer_matrix[1][6] = 72;
  answer_matrix[1][7] = 21;
  answer_matrix[1][8] = 29;
  answer_matrix[1][9] = 41;
  answer_matrix[1][10] = 12;
  answer_matrix[2][0] = 20;
  answer_matrix[2][1] = 15;
  answer_matrix[2][2] = 0;
  answer_matrix[2][3] = 15;
  answer_matrix[2][4] = 14;
  answer_matrix[2][5] = 25;
  answer_matrix[2][6] = 81;
  answer_matrix[2][7] = 9;
  answer_matrix[2][8] = 23;
  answer_matrix[2][9] = 27;
  answer_matrix[2][10] = 13;
  answer_matrix[3][0] = 21;
  answer_matrix[3][1] = 29;
  answer_matrix[3][2] = 15;
  answer_matrix[3][3] = 0;
  answer_matrix[3][4] = 4;
  answer_matrix[3][5] = 12;
  answer_matrix[3][6] = 92;
  answer_matrix[3][7] = 12;
  answer_matrix[3][8] = 25;
  answer_matrix[3][9] = 13;
  answer_matrix[3][10] = 25;
  answer_matrix[4][0] = 16;
  answer_matrix[4][1] = 28;
  answer_matrix[4][2] = 14;
  answer_matrix[4][3] = 4;
  answer_matrix[4][4] = 0;
  answer_matrix[4][5] = 16;
  answer_matrix[4][6] = 94;
  answer_matrix[4][7] = 9;
  answer_matrix[4][8] = 20;
  answer_matrix[4][9] = 16;
  answer_matrix[4][10] = 22;
  answer_matrix[5][0] = 31;
  answer_matrix[5][1] = 40;
  answer_matrix[5][2] = 25;
  answer_matrix[5][3] = 12;
  answer_matrix[5][4] = 16;
  answer_matrix[5][5] = 0;
  answer_matrix[5][6] = 95;
  answer_matrix[5][7] = 24;
  answer_matrix[5][8] = 36;
  answer_matrix[5][9] = 3;
  answer_matrix[5][10] = 37;
  answer_matrix[6][0] = 100;
  answer_matrix[6][1] = 72;
  answer_matrix[6][2] = 81;
  answer_matrix[6][3] = 92;
  answer_matrix[6][4] = 94;
  answer_matrix[6][5] = 95;
  answer_matrix[6][6] = 0;
  answer_matrix[6][7] = 90;
  answer_matrix[6][8] = 101;
  answer_matrix[6][9] = 99;
  answer_matrix[6][10] = 84;
  answer_matrix[7][0] = 12;
  answer_matrix[7][1] = 21;
  answer_matrix[7][2] = 9;
  answer_matrix[7][3] = 12;
  answer_matrix[7][4] = 9;
  answer_matrix[7][5] = 24;
  answer_matrix[7][6] = 90;
  answer_matrix[7][7] = 0;
  answer_matrix[7][8] = 15;
  answer_matrix[7][9] = 25;
  answer_matrix[7][10] = 13;
  answer_matrix[8][0] = 4;
  answer_matrix[8][1] = 29;
  answer_matrix[8][2] = 23;
  answer_matrix[8][3] = 25;
  answer_matrix[8][4] = 20;
  answer_matrix[8][5] = 36;
  answer_matrix[8][6] = 101;
  answer_matrix[8][7] = 15;
  answer_matrix[8][8] = 0;
  answer_matrix[8][9] = 35;
  answer_matrix[8][10] = 18;
  answer_matrix[9][0] = 31;
  answer_matrix[9][1] = 41;
  answer_matrix[9][2] = 27;
  answer_matrix[9][3] = 13;
  answer_matrix[9][4] = 16;
  answer_matrix[9][5] = 3;
  answer_matrix[9][6] = 99;
  answer_matrix[9][7] = 25;
  answer_matrix[9][8] = 35;
  answer_matrix[9][9] = 0;
  answer_matrix[9][10] = 38;
  answer_matrix[10][0] = 18;
  answer_matrix[10][1] = 12;
  answer_matrix[10][2] = 13;
  answer_matrix[10][3] = 25;
  answer_matrix[10][4] = 22;
  answer_matrix[10][5] = 37;
  answer_matrix[10][6] = 84;
  answer_matrix[10][7] = 13;
  answer_matrix[10][8] = 18;
  answer_matrix[10][9] = 38;
  answer_matrix[10][10] = 0;
  std::vector<std::vector<uint32_t>> result_matrix = graph.GetMatrix();
  for (uint32_t i = 0; i < graph.GetVertices(); ++i) {
    for (uint32_t j = 0; j < graph.GetVertices(); ++j) {
      EXPECT_EQ(result_matrix[i][j], answer_matrix[i][j]);
    }
  }
}

TEST(Graph, Load_From_File_Wrong) {
  s21::Graph graph;
  EXPECT_THROW(graph.LoadGraphFromFile("examples/valid_matrixxx.txt"),
               std::invalid_argument);
}

TEST(Graph, Load_From_File_Wrong_Data) {
  s21::Graph graph;
  EXPECT_THROW(graph.LoadGraphFromFile("examples/invalid_matrix.txt"),
               std::length_error);
}

TEST(Graph, Load_From_File_Wrong_Vertices_Amount) {
  s21::Graph graph;
  EXPECT_THROW(graph.LoadGraphFromFile("examples/invalid_matrix2.txt"),
               std::out_of_range);
}

TEST(Graph, Save_Graph_To_File) {
  s21::Graph graph;
  graph.LoadGraphFromFile("examples/valid_matrix.txt");
  graph.SaveGraphToFile("examples/valid_matrix_copy.txt");
  s21::Graph exported;
  exported.LoadGraphFromFile("examples/valid_matrix_copy.txt");
  std::vector<std::vector<uint32_t>> initial_matrix = graph.GetMatrix();
  std::vector<std::vector<uint32_t>> exported_matrix = exported.GetMatrix();
  for (uint32_t i = 0; i < graph.GetVertices(); ++i) {
    for (uint32_t j = 0; j < graph.GetVertices(); ++j) {
      EXPECT_EQ(initial_matrix[i][j], exported_matrix[i][j]);
    }
  }
}

TEST(Graph, Export_Graph_To_Dot) {
  s21::Graph graph;
  graph.LoadGraphFromFile("examples/valid_matrix.txt");
  graph.ExportGraphToDot("examples/valid_matrix.dot");
  graph.LoadGraphFromFile("examples/small_graph.txt");
  graph.ExportGraphToDot("examples/small_graph.dot");
}