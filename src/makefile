CPPflags= c++ -std=c++11

TESTPROG= testcode.out
#MAINPROG= mainprog.out
LIB = -fopenmp
all:	${TESTPROG} ${MAINPROG}

${TESTPROG} :	   		tests-main.o  test-functions.o
					${CPPflags} tests-main.o test-functions.o ${LIB} -o ${TESTPROG}

tests-main.o :			tests-main.cpp 
		        	${CPPflags} -c tests-main.cpp


test-functions.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp








