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

  ./mainprog.out 1
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
  pdflatex -synctex=1 -interaction=nonstopmode CompPhysProj3.tex
  bibtex CompPhysProj3.aux
  pdflatex -synctex=1 -interaction=nonstopmode CompPhysProj3.tex
  bibtex CompPhysProj3.aux
  pdflatex -synctex=1 -interaction=nonstopmode CompPhysProj3.tex
fi
