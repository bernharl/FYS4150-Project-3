UNAME := $(shell uname)

ifeq ($(UNAME), Linux) # Linux
$(info "Comiling for Linux")
CPPflags= g++ -std=c++11 -Wall
endif
ifeq ($(UNAME), Darwin) # Mac
$(info "Compiling for Mac using g++-9 for openMP support")
CPPflags= /usr/local/bin/g++-9 -std=c++11 -L/usr/local/Cellar/armadillo/9.600.5/lib/ -Wall
endif




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


${MAINPROG} :	   	main.cpp
					${CPPflags} main.cpp ${LIB} -o ${MAINPROG}





