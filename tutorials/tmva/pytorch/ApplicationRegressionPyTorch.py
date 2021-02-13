#!/usr/bin/env python
## \file
## \ingroup tutorial_tmva_pytorch
## \notebook -nodraw
## This tutorial shows how to apply a trained model to new data (regression).
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
if not isfile('tmva_reg_example.root'):
    call(['curl', '-O', 'http://root.cern.ch/files/tmva_reg_example.root'])

data = TFile.Open('tmva_reg_example.root')
tree = data.Get('TreeR')

branches = {}
for branch in tree.GetListOfBranches():
    branchName = branch.GetName()
    branches[branchName] = array('f', [-999])
    tree.SetBranchAddress(branchName, branches[branchName])
    if branchName != 'fvalue':
        reader.AddVariable(branchName, branches[branchName])


# Book methods
reader.BookMVA('PyTorch', TString('dataset/weights/TMVARegression_PyTorch.weights.xml'))


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


# Print some example regressions
print('Some example regressions:')
for i in range(20):
    tree.GetEntry(i)
    print('True/MVA value: {}/{}'.format(branches['fvalue'][0],reader.EvaluateMVA('PyTorch')))
