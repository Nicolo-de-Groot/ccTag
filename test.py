import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import tensorflow as tf
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
import json
from keras.models import model_from_json
import os
import sys

def read_data(fname,tname):
    input_data = pd.read_csv("./" + fname + ".csv")
    test_data = pd.read_csv("./" + tname + ".csv")

#    print(input_data.head(5))
    # make numpy arrays out of the pandas dataframes
    input_data = input_data.values
    test_data = test_data.values
    print(input_data.shape)
    # shuffle pions and muons
    np.random.seed(1234)
    np.random.shuffle(input_data)
    nevt = input_data.shape[0]

    train_x = input_data[:,:-1]
    train_y = input_data[:,-1]
    test_x = test_data[:,:-1]
    test_y = test_data[:,-1]

    return train_x, train_y, test_x, test_y

def normalize_data(dataset):
    mu = np.mean(dataset, axis = 0)
    sigma = np.std(dataset, axis = 0)
    return (dataset - mu)/sigma, mu, sigma

def normalize_test(dataset, mu, sigma):
    return (dataset - mu)/sigma

try:
    fname = sys.argv[1]
except IndexError:
    fname = "train3"
    print("No filename specified, using barrel")
try:
    tname = sys.argv[2]
except IndexError:
    tname = fname

json_file = open(fname+'.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
model = model_from_json(loaded_model_json, custom_objects={"GlorotUniform": tf.keras.initializers.glorot_uniform})

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# load model
model.load_weights(fname+".h5")
print("Model Loaded");

train_x, train_y, test_x, test_y = read_data(fname,tname)
train_x, mu, sigma = normalize_data(train_x)
test_x_n = normalize_test(test_x, mu, sigma)

#print results
test_loss, test_acc = model.evaluate(test_x_n, test_y)
y_pred = model.predict(test_x_n).ravel()
fpr, tpr, thresholds = roc_curve(test_y, y_pred)
inveff = 1/fpr
auc_cid = auc(fpr, tpr)
print('Test accuracy:', test_acc)
print('Test loss:', test_loss)
print('AUC', auc_cid)
plt.figure(1)
plt.subplot(221)
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr, tpr, label='ccTag (area = {:.3f})'.format(auc_cid))
#plt.plot(fpr_rf, tpr_rf, label='RF (area = {:.3f})'.format(auc_rf))
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc='best')
#plt.show()

# Zoom in view of the upper left corner.
#plt.figure(2)
plt.subplot(222)
plt.xlim(0, 0.5)
plt.ylim(0.5, 1)
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr, tpr, label='ccTag(area = {:.3f})'.format(auc_cid))
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve (zoomed in at top left)')
plt.legend(loc='best')
#plt.show()
plt.subplot(223)
plt.xlim(0.05, 0.9)
plt.plot(tpr, inveff)
plt.xlabel('Signal Efficiency')
plt.ylabel('1/ BG Efficiency')
plt.yscale('log')

plt.subplot(224)
y_sig = y_pred[test_y ==1]
y_bg = y_pred[test_y ==0]
n, bins, patches = plt.hist(x=y_pred, bins=50, color='#0504aa',
                            alpha=0.7, rwidth=1.00)
maxfreq = n.max()
n, bins, patches = plt.hist(x=y_bg, bins=50, color='#ff0400',
                            alpha=0.7, rwidth=1.00)
plt.xlabel('NNout')
plt.ylabel('#Events')
plt.ylim(top=np.ceil(maxfreq / 10) * 12 if maxfreq % 10 else maxfreq + 10)

plt.subplots_adjust(wspace=0.35, hspace=0.45)
plt.savefig(fname + ".png")
plt.show()
