CXX=gcc -lstdc++
FLAGS=-Wall -Werror -Wextra -std=c++17
GTEST=-lgtest
GCOVFLAGS=--coverage
DEBUG=-DDEBUG -g
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	OPEN_CMD = xdg-open
endif

ifeq ($(UNAME_S),Darwin)
	OPEN_CMD = open
endif

all: s21_matrix_oop.a test gcov

s21_matrix_oop.a: s21_matrix_oop.o
	@ar rcs s21_matrix_oop.a s21_matrix_oop.o

s21_matrix_oop.o: s21_matrix_oop.cc
	@$(CXX) $(FALGS) $(GTEST) $(DEBUG) $(GCOVFLAGS) -c s21_matrix_oop.cc

build_test: s21_matrix_oop.a
	@$(CXX) $(FLAGS) $(DEBUG) $(GCOVFLAGS) tests.cc $(GTEST) s21_matrix_oop.a -o test

main: s21_matrix_oop.a
	@$(CXX) $(FLAGS) $(GTEST) $(DEBUG) $(GCOVFLAGS) s21_matrix_oop.a main.cpp -o main

test: build_test
	./test

gcov: build_test
	./test
	@rm -rf tests.gcda tests.gcno
	@geninfo . -b . -o ./coverage.info
	@genhtml coverage.info -o gcov_report
	@$(OPEN_CMD) gcov_report/index.html

style:
	@clang-format -style=Google -n *.cc *.h 

clean:
	@-rm -rf *.o *.gcno *.gcda *.info *.gcov *.a *.dSYM  report test
