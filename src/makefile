CPPflags = g++-9 -std=c++11
MAINPROG = mainprog.out
TESTPROG = testcode.out
#MAINPROG= mainprog.out
LIB = -fopenmp
all:	${TESTPROG} ${MAINPROG}

${TESTPROG} :	   		tests-main.o  test-functions.o
					${CPPflags} tests-main.o test-functions.o ${LIB} -o ${TESTPROG}

tests-main.o :			tests-main.cpp 
		        	${CPPflags} -c tests-main.cpp


test-functions.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp


${MAINPROG} :	   	main.cpp
					${CPPflags} main.cpp ${LIB} -o ${MAINPROG}





