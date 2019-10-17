# FYS4150-Project-3
Third project in FYS4150 Computational Physics

## Dependencies
* MacOS: The newest version of homebrew GCC (You need to be able to run g++-9).
* Does not work on Windows
* OpenMP is used as the parallelisation backend.

## Build Instructions:

* Run main_script.sh, this will build the source code, run tests and the main executable, generate plots using the Python library Matplotlib, and (optional) build the latex report for the project.
* Should you want to run compile the c++ source code manually, run "make" in src/, this should generate two executables, one for unit tests using Catch2 (https://github.com/catchorg/Catch2) and one for running the main calculations.
&nbsp;

## Run Instructions:
* src/mainprog.out Is compiled when running make. Run this to do all calculations and generate necessary data for plotting
* src/testcode.out Run this for unit tests
&nbsp;

## Structure: 

* src/integration.cpp Contains all integration methods used.
* src/main.cpp Runs the methods in src/integration.cpp. Results are saved in src/montecarlo_paro(1-3).txt, src/montecarlo.txt, src/montecarlo_improved.txt and Exercise_a_b.txt
* src/test-functions.cpp Contains all unit tests run by Catch2.
* src/plot1a.py Generates all relevant plots for gauss-legendre and gauss-laguerre. Plots are saved as vector graphics (.pdf) in doc/Figures
* src/plotmontecarlo.py Plots all relevant plots for monte carlo integration. Plots are saved as vector graphics (.pdf) in doc/Figures
&nbsp;
 
