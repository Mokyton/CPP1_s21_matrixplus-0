CFLAGS = -Wall -Werror -Wextra

all: s21_matrix_oop.a

 build: S21Matrix.cc
	gcc -c $(CFLAGS) -lstdc++ -std=c++17 S21Matrix.cc -o S21Matrix.o

s21_matrix_oop.a: build
	ar rcs s21_matrix_oop.a S21Matrix.o
	ranlib s21_matrix_oop.a


clean:
	rm -rf *.0 *.out s21_matrix_oop.a
	rm -rf *.gcda *.gcno report gcov_report.* gcov_report

rb:
	make clean
	make build
	./a.out

leaks: test
	leaks -atExit -- ./test

test: s21_matrix_oop.a test.cc
	gcc -lgtest -lstdc++ -std=c++17 s21_matrix_oop.a test.cc -o test
	./test

gcov_report:
	gcc  --coverage test.cc S21Matrix.cc -o gcov_report -lgtest -lstdc++
	./gcov_report
	lcov -t "./gcov_report" -o gcov_report.info -c -d .
	genhtml -o report gcov_report.info
	open ./report/index.html
