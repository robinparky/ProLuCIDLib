from __future__ import division
import sys
import pickle
import time
import math
import gc



import tensorflow as tf
import mesh-tensorflow as mtf
from tensorflow import keras
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 5 :
    print(len(sys.argv))
    print("Error with command line inputs")
    sys.exit(0)
else:
    save = sys.argv[1]
    outputPath = sys.argv[2]
    batchSize = int(sys.argv[3])
    numEpochs = int(sys.argv[4])

start = time.time()

print("Running program")
print("Training Data with ", numEpochs, " epochs.")
print("\n")

'''--------------------------------------------------------------------------'''
print("Gathering Data from Files")
with tf.device('/cpu:0'):
    with open (save + 'spectrumsList', 'rb') as sp:
        spectrums = pickle.load(sp)
    sp.close()
    with open (save + 'labelList', 'rb') as lp:
        labelList = pickle.load(lp)
    sp.close()
    with open (save + 'indexList', 'rb') as lp:
        indexList = pickle.load(lp)
    sp.close()
    with open (save + 'idList', 'rb') as lp:
        idList = pickle.load(lp)
    sp.close()
    with open (save + 'binArray', 'rb') as lp:
        binArray = pickle.load(lp)
    sp.close()
    with open (save + 'outputLabels', 'rb') as lp:
        labels = np.array(pickle.load(lp))
    sp.close()
'''--------------------------------------------------------------------------'''


#Metrics for printing(mostly)
inputs = len(spectrums)
totalBins = len(binArray[0])
inputLayers = len(binArray)
outputLayers = len(labels)


print("Inputs (spectrums): ",inputs)
print("\tShapes: ", binArray.shape);
print("\tType: ", binArray.dtype);
print("\tInput Layers: ", inputLayers, "\n")

print("Outputs:", )
print("\tNum: ", len(labels));
print("\tType: ", labels.dtype);
print("\tOutput Layers: ", outputLayers)


graph = mtf.Graph()
mesh = mtf.Mesh(graph, "my_mesh")
batch_dim = mtf.Dimension("batch", inputs)
bin_dim = mtf.Dimension("specSize", totalBins)
hidden_dim = mtf.Dimension("hidden", outputLayers * 3)
classes_dim = mtf.Dimension("classes", len(labels))

spectrumTensor = mtf.import_tf_tensor(mesh, spectrums, shape=[batch_dim, bin_dim])
labels = mtf.import_tf_tensor(mesh, labels, [batch_dim]

w1 = mtf.get_variable(mesh, "w1", [bin_dim, hidden_dim])
w2 = mtf.get_variable(mesh, "w2", [hidden_dim, hidden_dim])
w2 = mtf.get_variable(mesh, "w2", [hidden_dim, classes_dim])

h1 = mtf.relu(mtf.einsum(spectrums, w1, output_shape=[batch_dim, hidden_dim]))
h2 = mtf.relu(mtf.einsum(h1, w2, output_shape=[batch_dim, hidden_dim]))
logits = mtf.einsum(h2, w3, output_shape=[batch_dim, classes_dim])

loss = mtf.reduce_mean(mtf.layers.softmax_cross_entropy_with_logits(logits, mtf.one_hot(labels, classes_dim), classes_dim))

w1_grad, w2_grad, w3_grad = mtf.gradients([loss], [w1, w2, w3])

update_w1_op = mtf.assign(w1, w1 - w1_grad * 0.001)
update_w2_op = mtf.assign(w2, w2 - w2_grad * 0.001)
update_w3_op = mtf.assign(w3, w3 - w3_grad * 0.001)


devices = ["gpu:0", "gpu:1", "gpu:2", "gpu:3"]
mesh_shape = [("all_processors", 4)]
layout_rules = [("hidden", "all_processors")]
mesh_impl = mtf.placement_mesh_impl.PlacementMeshImpl(mesh_shape, layout_rules, devices)
lowering = mtf.Lowering(graph, {mesh:mesh_impl})
tf_update_ops = [lowering.lowered_operation(update_w1_op), lowering.lowered_operation(update_w2_op)]

'''--------------------------------------------------------------------------'''

#Train the data

#model.fit_generator(generator(batchSize), steps_per_epoch = np.ceil(len(binArray)/batchSize), epochs = numEpochs)

model.fit(binArray, indexList, batch_size = batchSize, epochs = numEpochs)
#model.save_weights(outputPath)
model.save(outputPath)


print ("Finished Training")
train = time.time()
print ("Elapsed Time: " + str(round(train - start)) + "\n")
