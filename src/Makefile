CXX=g++
CXX_FLAGS=-std=c++17 -Wall -Wextra -Werror -O2
SAN_FLAGS=-fsanitize=address -pedantic 
COV_FLAGS=--coverage 
LD_FlAGS=-lgtest

ALL_SOURCES = ./*.cc ./helpers/*.tpp
ALL_HEADERS = ./*.h ./helpers/*.h

SOURCE = s21_graph_algorithms.cc s21_graph.cc
LIB_FILES = $(SOURCE:.cc=.a)

TEST_DIR = ./tests
TEST = ./tests/test
TESTS = $(wildcard $(TEST_DIR)/*.cc)

all: clean test

s21_graph.a: s21_graph.h s21_graph.cc
	$(CXX) $(CXX_FLAGS) -c s21_graph.cc
	@ar rcs s21_graph.a s21_graph.o

s21_graph_algorithms.a: s21_graph_algorithms.h s21_graph_algorithms.cc
	$(CXX) $(CXX_FLAGS) -c s21_graph_algorithms.cc
	@ar rcs s21_graph_algorithms.a s21_graph_algorithms.o

mainapp: $(LIB_FILES)
	$(CXX) $(CXX_FLAGS) $(SAN_FLAGS) main.cc s21_main_application.cc $(LIB_FILES) -o mainapp.out
	./mainapp.out

test: clean $(LIB_FILES)
	$(CXX) $(CXX_FLAGS) $(SAN_FLAGS) $(TESTS) $(LIB_FILES) $(LD_FlAGS) -o $(TEST)
	$(TEST)

gcov_report: clean
	$(CXX) $(CXX_FLAGS) $(COV_FLAGS) $(TESTS) $(SOURCE) $(LD_FlAGS) -o $(TEST)
	$(TEST)
	mkdir ./tests/report/
	gcovr . -s --html --html-details -o ./tests/report/coverage_report.html
	rm -rf *.gcno *.gcda ./tests/*.gcno ./tests/*.gcda *.o
	open ./tests/report/coverage_report.html

style:
	cppcheck $(ALL_SOURCES) $(ALL_HEADERS) --language=c++
	cp ../materials/linters/.clang-format .clang-format
	clang-format -i $(ALL_SOURCES) $(ALL_HEADERS)
	clang-format -n $(ALL_SOURCES) $(ALL_HEADERS)
	rm -rf .clang-format

leaks_check: clean test
ifeq ($(OS), Darwin)
	leaks --atExit -- $(TEST)
else
	valgrind --leak-check=full --log-file=RESULT_VALGRIND.txt $(TEST)
endif

clean:
	@rm -rf *.o *.a *.out $(TEST) RESULT_VALGRIND.txt ./examples/*copy.txt ./examples/*.dot coverage*
	@rm -rf *.gcno *.gcda ./tests/*.gcno ./tests/*.gcda ./tests/report
