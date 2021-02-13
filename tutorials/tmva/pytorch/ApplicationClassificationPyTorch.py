#!/usr/bin/env python
## \file
## \ingroup tutorial_tmva_pytorch
## \notebook -nodraw
## This tutorial shows how to apply a trained model to new data.
##
## \macro_code
##
## \date 2020
## \author Anirudh Dagar <anirudhdagar6@gmail.com> - IIT, Roorkee


from ROOT import TMVA, TFile, TString
from array import array
from subprocess import call
from os.path import isfile


# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()
reader = TMVA.Reader("Color:!Silent")


# Load data
if not isfile('tmva_class_example.root'):
    call(['curl', '-O', 'http://root.cern.ch/files/tmva_class_example.root'])

data = TFile.Open('tmva_class_example.root')
signal = data.Get('TreeS')
background = data.Get('TreeB')

branches = {}
for branch in signal.GetListOfBranches():
    branchName = branch.GetName()
    branches[branchName] = array('f', [-999])
    reader.AddVariable(branchName, branches[branchName])
    signal.SetBranchAddress(branchName, branches[branchName])
    background.SetBranchAddress(branchName, branches[branchName])


# Define predict function
def predict(model, test_X, batch_size=32):
    # Set to eval mode
    model.eval()
   
    test_dataset = torch.utils.data.TensorDataset(torch.Tensor(test_X))
    test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    predictions = []
    with torch.no_grad():
        for i, data in enumerate(test_loader):
            X = data[0]
            outputs = model(X)
            predictions.append(outputs)
        preds = torch.cat(predictions)
   
    return preds.numpy()


load_model_custom_objects = {"optimizer": None, "criterion": None, "train_func": None, "predict_func": predict}


# Book methods
reader.BookMVA('PyTorch', TString('dataset/weights/TMVAClassification_PyTorch.weights.xml'))


# Print some example classifications
print('Some signal example classifications:')
for i in range(20):
    signal.GetEntry(i)
    print(reader.EvaluateMVA('PyTorch'))
print('')

print('Some background example classifications:')
for i in range(20):
    background.GetEntry(i)
    print(reader.EvaluateMVA('PyTorch'))
