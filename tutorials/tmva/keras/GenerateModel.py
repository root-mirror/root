#!/usr/bin/env python
## \file
## \ingroup tutorial_tmva_keras
## \notebook -nodraw
## This tutorial shows how to define and generate a keras model for use with
## TMVA.
##
## \macro_code
##
## \date 2017
## \author TMVA Team

from keras.models import Sequential
from keras.layers.core import Dense, Activation
from keras.regularizers import l2
from keras.optimizers import SGD
from keras.utils import plot_model

# Setup the model here
num_input_nodes = 4
num_output_nodes = 2
num_hidden_layers = 1
nodes_hidden_layer = 64
l2_val = 1e-5

model = Sequential()

# Hidden layer 1
# NOTE: Number of input nodes need to be defined in this layer
model.add(Dense(nodes_hidden_layer, activation='relu', W_regularizer=l2(l2_val), input_dim=num_input_nodes))

# Hidden layer 2 to num_hidden_layers
# NOTE: Here, you can do what you want
for k in range(num_hidden_layers-1):
    model.add(Dense(nodes_hidden_layer, activation='relu', W_regularizer=l2(l2_val)))

# Ouput layer
# NOTE: Use following output types for the different tasks
# Binary classification: 2 output nodes with 'softmax' activation
# Regression: 1 output with any activation ('linear' recommended)
# Multiclass classification: (number of classes) output nodes with 'softmax' activation
model.add(Dense(num_output_nodes, activation='softmax'))

# Compile model
# NOTE: Use following settings for the different tasks
# Any classification: 'categorical_crossentropy' is recommended loss function
# Regression: 'mean_squared_error' is recommended loss function
model.compile(loss='categorical_crossentropy', optimizer=SGD(lr=0.01), metrics=['accuracy',])

# Save model
model.save('model.h5')

# Additional information about the model
# NOTE: This is not needed to run the model

# Print summary
model.summary()

# Visualize model as graph
try:
    plot_model(model, to_file='model.png', show_shapes=True)
except:
    print('[INFO] Failed to make model plot')
