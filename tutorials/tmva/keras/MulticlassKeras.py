#!/usr/bin/env python

from ROOT import TMVA, TFile, TTree, TCut, gROOT
from os.path import isfile

from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.regularizers import l2
from keras.optimizers import SGD

# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

output = TFile.Open('TMVAMulticlass.root', 'RECREATE')
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

dataloader = TMVA.DataLoader('TMVAMulticlass')
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
model = Sequential()
model.add(Dense(32, activation='relu', W_regularizer=l2(1e-5), input_dim=4))
#model.add(Dense(32, activation='relu', W_regularizer=l2(1e-5)))
model.add(Dense(4, activation='softmax'))

# Set loss and optimizer
model.compile(loss='categorical_crossentropy', optimizer=SGD(lr=0.01), metrics=['accuracy',])

# Store model to file
model.save('model.h5')
model.summary()

# Book methods
factory.BookMethod(dataloader, TMVA.Types.kFisher, 'Fisher',
        '!H:!V:Fisher:VarTransform=D,G')
factory.BookMethod(dataloader, TMVA.Types.kPyKeras, "PyKeras",
        'H:!V:VarTransform=D,G:FilenameModel=model.h5:NumEpochs=20:BatchSize=32')

# Run TMVA
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
