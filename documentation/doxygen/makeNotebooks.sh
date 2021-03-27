#!/bin/bash
# Run commands to generate notebooks, possibly in parallel.
# List dependencies below as <tutorial> <tutorial it depends on>


while [ $# -gt 0 ]; do
  case $1 in
    -j*)
      nJobs=${1#-j}
      shift
      ;;
    *)
      inputFile=$1
      shift
      ;;
  esac
done

if [ ! -f "$inputFile" ]; then
  echo "Usage: makeNotebooks <fileWithNotebooks to generate> [-j jobs]"
  exit 1
fi

# Prepare jupyter runs
mkdir -p ~/.jupyter

# If notebooks depend on the output of others, run those first.
# List dependencies below in the form
# <Notebook with dependencies> <dependency> [<dependency> ...]
while read notebook dependencies; do
	if grep -q ${notebook} $inputFile; then
    for dependency in $dependencies; do
      cmd=$(grep $dependency $inputFile)
      if [ -n "$cmd" ]; then
        echo "Running $cmd as depedency of $notebook"
        ${PYTHON_EXECUTABLE:-python3} $cmd  && sed -i'.back' "\#${dependency}#d" $inputFile
      fi
    done
  fi
done <<EOF
math/testUnfold5d.C   math/testUnfold5a.C math/testUnfold5b.C math/testUnfold5c.C
xml/xmlreadfile.C   xml/xmlnewfile.C
roofit/rf503_wspaceread.C   roofit/rf502_wspacewrite.C roofit/rf502_wspacewrite.py
io/readCode.C   io/importCode.C
fit/fit1.C   hist/fillrandom.C
fit/myfit.C   fit/fitslicesy.C
foam/foam_demopers.C   foam/foam_demo.C
tree/staff.C   tree/cernbuild.C
tree/cernstaff.C   tree/cernbuild.C
hist/hbars.C   tree/cernbuild.C
pyroot/ntuple1.py   pyroot/hsimple.py
pyroot/h1draw.py   pyroot/hsimple.py
pyroot/fit1.py   pyroot/fillrandom.py
tmva/TMVAClassificationApplication.C   tmva/TMVAClassification.C
tmva/TMVAClassificationCategory.C   tmva/TMVAClassification.C
tmva/TMVAClassificationCategoryApplication.C   tmva/TMVAClassificationCategory.C
tmva/TMVAMulticlass.C   tmva/TMVAMultipleBackgroundExample.C
tmva/TMVAMulticlassApplication.C   tmva/TMVAMulticlass.C
tmva/TMVARegressionApplication.C   tmva/TMVARegression.C
EOF

# Run rest in parallel
xargs -L 1 -P ${nJobs:-1} ${PYTHON_EXECUTABLE:-python3} < $inputFile

rm ${inputFile}.back

