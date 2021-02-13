#!/usr/bin/env python
## \file
## \ingroup tutorial_tmva_pytorch
## \notebook -nodraw
## This tutorial shows how to do multiclass classification in TMVA with neural
## networks trained with PyTorch.
##
## \macro_code
##
## \date 2020
## \author Anirudh Dagar <anirudhdagar6@gmail.com> - IIT, Roorkee


from ROOT import TMVA, TFile, TTree, TCut, gROOT
from os.path import isfile

import torch
from torch import nn


# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

output = TFile.Open('TMVA.root', 'RECREATE')
factory = TMVA.Factory('TMVAClassification', output,
    '!V:!Silent:Color:DrawProgressBar:Transformations=D,G:AnalysisType=multiclass')


# Load data
if not isfile('tmva_example_multiple_background.root'):
    createDataMacro = str(gROOT.GetTutorialDir()) + '/tmva/createData.C'
    print(createDataMacro)
    gROOT.ProcessLine('.L {}'.format(createDataMacro))
    gROOT.ProcessLine('create_MultipleBackground(4000)')

data = TFile.Open('tmva_example_multiple_background.root')
signal = data.Get('TreeS')
background0 = data.Get('TreeB0')
background1 = data.Get('TreeB1')
background2 = data.Get('TreeB2')

dataloader = TMVA.DataLoader('dataset')
for branch in signal.GetListOfBranches():
    dataloader.AddVariable(branch.GetName())

dataloader.AddTree(signal, 'Signal')
dataloader.AddTree(background0, 'Background_0')
dataloader.AddTree(background1, 'Background_1')
dataloader.AddTree(background2, 'Background_2')
dataloader.PrepareTrainingAndTestTree(TCut(''),
        'SplitMode=Random:NormMode=NumEvents:!V')


# Generate model
# Define model
model = nn.Sequential()
model.add_module('linear_1', nn.Linear(in_features=4, out_features=32))
model.add_module('relu', nn.ReLU())
model.add_module('linear_2', nn.Linear(in_features=32, out_features=4))
model.add_module('softmax', nn.Softmax(dim=1))


# Set loss and optimizer
loss = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD


# Define train function
def train(model, train_loader, val_loader, num_epochs, batch_size, optimizer, criterion, save_best, scheduler):
    trainer = optimizer(model.parameters(), lr=0.01)
    schedule, schedulerSteps = scheduler
    best_val = None

    for epoch in range(num_epochs):
        # Training Loop
        # Set to train mode
        model.train()
        running_train_loss = 0.0
        running_val_loss = 0.0
        for i, (X, y) in enumerate(train_loader):
            trainer.zero_grad()
            output = model(X)
            target = torch.max(y, 1)[1]
            train_loss = criterion(output, target)
            train_loss.backward()
            trainer.step()

            # print train statistics
            running_train_loss += train_loss.item()
            if i % 32 == 31:    # print every 32 mini-batches
                print("[{}, {}] train loss: {:.3f}".format(epoch+1, i+1, running_train_loss / 32))
                running_train_loss = 0.0

        if schedule:
            schedule(optimizer, epoch, schedulerSteps)

        # Validation Loop
        # Set to eval mode
        model.eval()
        with torch.no_grad():
            for i, (X, y) in enumerate(val_loader):
                output = model(X)
                target = torch.max(y, 1)[1]
                val_loss = criterion(output, target)
                running_val_loss += val_loss.item()

            curr_val = running_val_loss / len(val_loader)
            if save_best:
               if best_val==None:
                   best_val = curr_val
               best_val = save_best(model, curr_val, best_val)

            # print val statistics per epoch
            print("[{}] val loss: {:.3f}".format(epoch+1, curr_val))
            running_val_loss = 0.0

    print("Finished Training on {} Epochs!".format(epoch+1))

    return model


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


load_model_custom_objects = {"optimizer": optimizer, "criterion": loss, "train_func": train, "predict_func": predict}


# Store model to file
# Convert the model to torchscript before saving
m = torch.jit.script(model)
torch.jit.save(m, "model.pt")
print(m)


# Book methods
factory.BookMethod(dataloader, TMVA.Types.kFisher, 'Fisher',
        '!H:!V:Fisher:VarTransform=D,G')
factory.BookMethod(dataloader, TMVA.Types.kPyTorch, "PyTorch",
        'H:!V:VarTransform=D,G:FilenameModel=model.pt:NumEpochs=20:BatchSize=32')


# Run TMVA
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

# Plot ROC Curves
roc = factory.GetROCCurve(dataloader)
roc.SaveAs('ROC_MulticlassPyTorch.png')