# JupyROOT
A software layer to integrate Jupyter notebooks and ROOT.

## Installation
1. [Install ROOT6](https://root.cern.ch/building-root) (> 6.05)
2. Install dependencies: pip install jupyter metakernel

## Start using ROOTbooks
Set up the ROOT environment (`. $ROOTSYS/bin/thisroot.[c]sh`) and type in your
shell:
```
root --notebook
```
This will start a ROOT-flavoured notebook server in your computer.

Alternatively, if you would like to use the Jupyter command directly, you 
can do on Linux:
```
cp -r $ROOTSYS/etc/notebook/kernels/root ~/.local/share/jupyter/kernels
jupyter notebook
```
and on OSx:
```
cp -r $ROOTSYS/etc/notebook/kernels/root /Users/$USER/Library/Jupyter/kernels/
jupyter notebook
```

Once the server is up, you can use ROOT with two kernels:

1. ROOT C++: new kernel provided by ROOT
2. Python 2: already provided by Jupyter

##  C++ ROOTbook
ROOT offers a C++ kernel that transforms the notebook in a ROOT prompt.
Embedded graphics, syntax highlighting and tab completion are among
the features provided by this kernel.

An example of how you would plot a histogram in a C++ ROOTbook is:
```cpp
TCanvas c;
TH1F h("h","ROOT Histo;X;Y",64,-4,4);
h.FillRandom("gaus");
h.Draw();
c.Draw();
``` 

## Python ROOTbook
If you prefer to use Python, you can create a new Python 2 kernel and
import the ROOT libraries:
```python
import ROOT
```
And then you could write something like:
```python
c = ROOT.TCanvas("c")
h = ROOT.TH1F("h","ROOT Histo;X;Y",64,-4,4)
```
Additionally, you can mix Python and C++ in the same notebook
by using the **%%cpp** magic:
```cpp
%%cpp
h->FillRandom("gaus");
h->Draw();
c->Draw();
```
