#!/usr/bin/env python
## \file
## \ingroup tutorial_tmva_keras
## \notebook -nodraw
## This tutorial shows how to do regression in TMVA with neural networks
## trained with keras.
##
## \macro_code
##
## \date 2017
## \author TMVA Team

from ROOT import TMVA, TFile, TTree, TCut
from subprocess import call
from os.path import isfile

from keras.models import Sequential
from keras.layers.core import Dense, Activation
from keras.regularizers import l2
from keras.optimizers import SGD

# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

output = TFile.Open('TMVA.root', 'RECREATE')
factory = TMVA.Factory('TMVARegression', output,
        '!V:!Silent:Color:DrawProgressBar:Transformations=D,G:AnalysisType=Regression')

# Load data
if not isfile('tmva_reg_example.root'):
    call(['curl', '-O', 'http://root.cern/files/tmva_reg_example.root'])

data = TFile.Open('tmva_reg_example.root')
tree = data.Get('TreeR')

dataloader = TMVA.DataLoader('dataset')
for branch in tree.GetListOfBranches():
    name = branch.GetName()
    if name != 'fvalue':
        dataloader.AddVariable(name)
dataloader.AddTarget('fvalue')

dataloader.AddRegressionTree(tree, 1.0)
dataloader.PrepareTrainingAndTestTree(TCut(''),
        'nTrain_Regression=4000:SplitMode=Random:NormMode=NumEvents:!V')

# Generate model

# Define model
model = Sequential()
model.add(Dense(64, activation='tanh', W_regularizer=l2(1e-5), input_dim=2))
model.add(Dense(1, activation='linear'))

# Set loss and optimizer
model.compile(loss='mean_squared_error', optimizer=SGD(lr=0.01))

# Store model to file
model.save('model.h5')
model.summary()

# Book methods
factory.BookMethod(dataloader, TMVA.Types.kPyKeras, 'PyKeras',
        'H:!V:VarTransform=D,G:FilenameModel=model.h5:NumEpochs=20:BatchSize=32')
factory.BookMethod(dataloader, TMVA.Types.kBDT, 'BDTG',
        '!H:!V:VarTransform=D,G:NTrees=1000:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=4')

# Run TMVA
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
