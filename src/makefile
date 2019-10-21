UNAME := $(shell uname)

ifeq ($(UNAME), Linux) # Linux
$(info "Compiling for Linux")
CPPflags= g++ -std=c++11 -Wall -O3
LIB = -fopenmp
endif
ifeq ($(UNAME), Darwin) # Mac
$(info "Compiling for Mac for openMP support")
CPPflags= g++ -std=c++11 -Wall -O3
LIB = -lomp
endif




TESTPROG= testcode.out
MAINPROG= mainprog.out

all:	${TESTPROG} ${MAINPROG}

${TESTPROG} :	   		tests-main.o  test-functions.o
					${CPPflags} tests-main.o test-functions.o ${LIB} -o ${TESTPROG}

tests-main.o :			tests-main.cpp
		        	${CPPflags} -c tests-main.cpp ${LIB}

test-functions.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp ${LIB}


${MAINPROG} :	   	main.cpp
					${CPPflags} main.cpp ${LIB} -o ${MAINPROG}
