CPPflags= c++ -std=c++11

TESTPROG= testcode.out
#MAINPROG= mainprog.out
LIB = -fopenmp
all:	${TESTPROG} ${MAINPROG}

${TESTPROG} :	   		tests-main.o  test-functions.o integration.o
					${CPPflags} tests-main.o test-functions.o integration.o ${LIB} -o ${TESTPROG}

tests-main.o :			tests-main.cpp 
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp

integration.o :		integration.cpp 
		        	${CPPflags} -c integration.cpp






