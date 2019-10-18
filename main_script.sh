#!/bin/bash
cd src

make # Compile everything
echo "Run tests? (y/n)"
read yntest
if [ "$yntest" == "y" ] # If y, run tests
then
  ./testcode.out
fi


echo "Generate new data? (y/n)"
read yn
if [ "$yn" == "y" ] # If y, run c++ code
then

  ./mainprog.out

fi

echo "Generate compiler flag data? (y/n)"
read ynflag
if [ "$ynflag" == "y" ]
then
  os=$(uname)
  if [ "$os" == "Darwin" ]
  then
    echo "Compiling for MacOS"
    g++-9 monte_carlo_compiler_flags.cpp -std=c++11 -fopenmp -O1 -o monte_carlo_compiler_flags_O1.out
    ./monte_carlo_compiler_flags_O1.out 4 1

    g++-9 monte_carlo_compiler_flags.cpp -std=c++11 -fopenmp -O2 -o monte_carlo_compiler_flags_O2.out
    ./monte_carlo_compiler_flags_O2.out 4 2

    g++-9 monte_carlo_compiler_flags.cpp -std=c++11 -fopenmp -O3 -o monte_carlo_compiler_flags_O3.out
    ./monte_carlo_compiler_flags_O3.out 4 3
  elif [ "$os" == "Linux"]
  then
    echo "Compiling for Linux"
    g++ monte_carlo_compiler_flags.cpp -std=c++11 -fopenmp -O1 -o monte_carlo_compiler_flags_O1.out
    ./monte_carlo_compiler_flags_O1.out 4 1

    g++ monte_carlo_compiler_flags.cpp -std=c++11 -fopenmp -O2 -o monte_carlo_compiler_flags_O2.out
    ./monte_carlo_compiler_flags_O2.out 4 2

    g++ monte_carlo_compiler_flags.cpp -std=c++11 -fopenmp -O3 -o monte_carlo_compiler_flags_O3.out
    ./monte_carlo_compiler_flags_O3.out 4 3
  else
    echo "Your OS is not currently supported"
  fi

fi

echo "Generate plots? (y/n)"
read ynplot
if [ "$ynplot" == "y" ] # If y, run python script for plots
then
  python3 plot1a.py
  python3 plotmontecarlo.py
fi

echo "Build report? (y/n)"
read ynreport
# If y, compile TeX document. The compilation is run many times because
# bibtex is usually non-cooperative...
if [ "$ynreport" == "y" ]
then
  cd ../doc/
  pdflatex -synctex=1 -interaction=nonstopmode ComphysProj3.tex
  bibtex Comphysproj3.aux
  pdflatex -synctex=1 -interaction=nonstopmode ComphysProj3.tex
  bibtex ComphysProj3.aux
  pdflatex -synctex=1 -interaction=nonstopmode ComphysProj3.tex
fi
