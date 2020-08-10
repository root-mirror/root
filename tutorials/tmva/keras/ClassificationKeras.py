#!/usr/bin/env python
## \file
## \ingroup tutorial_tmva_keras
## \notebook -nodraw
## This tutorial shows how to do classification in TMVA with neural networks
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
from keras.layers import Dense, Activation
from keras.regularizers import l2
from keras.optimizers import SGD

# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

output = TFile.Open('TMVA.root', 'RECREATE')
factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:Transformations=D,G:AnalysisType=Classification')

# Load data
if not isfile('tmva_class_example.root'):
    call(['curl', '-O', 'http://root.cern/files/tmva_class_example.root'])

data = TFile.Open('tmva_class_example.root')
signal = data.Get('TreeS')
background = data.Get('TreeB')

dataloader = TMVA.DataLoader('dataset')
for branch in signal.GetListOfBranches():
    dataloader.AddVariable(branch.GetName())

dataloader.AddSignalTree(signal, 1.0)
dataloader.AddBackgroundTree(background, 1.0)
dataloader.PrepareTrainingAndTestTree(TCut(''),
                                      'nTrain_Signal=4000:nTrain_Background=4000:SplitMode=Random:NormMode=NumEvents:!V')

# Generate model

# Define model
model = Sequential()
model.add(Dense(64, activation='relu', W_regularizer=l2(1e-5), input_dim=4))
model.add(Dense(2, activation='softmax'))

# Set loss and optimizer
model.compile(loss='categorical_crossentropy',
              optimizer=SGD(lr=0.01), metrics=['accuracy', ])

# Store model to file
model.save('model.h5')
model.summary()

# Book methods
factory.BookMethod(dataloader, TMVA.Types.kFisher, 'Fisher',
                   '!H:!V:Fisher:VarTransform=D,G')
factory.BookMethod(dataloader, TMVA.Types.kPyKeras, 'PyKeras',
                   'H:!V:VarTransform=D,G:FilenameModel=model.h5:NumEpochs=20:BatchSize=32')

# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
