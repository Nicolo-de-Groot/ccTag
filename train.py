import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import tensorflow as tf
#from keras.models import save_model
import json
import os
import sys
import datetime

def read_data(fname, tname):
    input_data = pd.read_csv("./" + fname + ".csv")
    test_data = pd.read_csv("./" + tname + ".csv")

    print(input_data.head(5))
    # make numpy arrays out of the pandas dataframes
    input_data = input_data.values
    test_data = test_data.values
    print(input_data.shape)
    # shuffle pions and muons
    rand=1234
    np.random.seed(rand)
    np.random.shuffle(input_data)
    np.random.shuffle(test_data)
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
    tname = sys.argv[2]
except IndexError:
    fname = "train3"
    fname = "test3"
    print("No filename specified, using train3")

#use hyperparameter optimation to fix these numbers to a good fit (see bookmarks)
#originally 3 relu, 1 sigmoid; 10 - 0.1, 9 - 0.1, 6 - 0.2, 1; epochs = 15
dense1 = [20, 0.2]
dense2 = [16, 0.2]
dense3 = [12, 0.2]
dense4= [10, 0.2]
dense5 = [1, 0]
activ = [tf.nn.relu, tf.nn.relu, tf.nn.relu, tf.nn.relu, tf.nn.sigmoid]

#see log_keras function

model = tf.keras.models.Sequential([
    tf.keras.layers.Dense(dense1[0], activation=activ[0]),
    tf.keras.layers.Dropout(dense1[1]),
    tf.keras.layers.Dense(dense2[0], activation=activ[1]),
    tf.keras.layers.Dropout(dense2[1]),
#    tf.keras.layers.Dense(dense3[0], activation=activ[2]),
#    tf.keras.layers.Dropout(dense3[1]),
    tf.keras.layers.Dense(dense4[0], activation=activ[3]),
    tf.keras.layers.Dropout(dense4[1]),
    tf.keras.layers.Dense(dense5[0], activation=activ[4])
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
train_x, train_y, test_x, test_y = read_data(fname,tname)
train_x, mu, sigma = normalize_data(train_x)
test_x_n = normalize_test(test_x, mu, sigma)
epochsNum = 16
history = model.fit(train_x, train_y, batch_size=64, epochs=epochsNum, validation_data = (test_x_n,test_y))

model.summary()
#print results
test_loss, test_acc = model.evaluate(test_x_n, test_y)
print('Test accuracy:', test_acc)
print('Test loss:', test_loss)
plt.plot(history.history['acc'])
plt.plot(history.history['val_acc'])
plt.title("Model accuracy")
plt.ylabel("Accuracy")
plt.xlabel("Epoch")
plt.savefig('acc_plot_%s.png' % (fname))
#plt.show()

# save model

with open(fname + '.json', "w") as json_file:
    json_file.write(model.to_json())
model.save_weights(fname + '.h5', save_format='h5')
print("model saved")

def log_keras(log_name='keras_log'):
   with open('{0:s}.out'.format(log_name),'a') as log:
      log.write(str(datetime.datetime.now()) + '\n')
      iterArray = np.append(dense1, str(activ[0]))
      iterArray = np.append(iterArray, np.append(dense2, str(activ[1])))
      iterArray = np.append(iterArray, np.append(dense3, str(activ[2])))
      iterArray = np.append(iterArray, np.append(dense4, str(activ[3])))
      iterArray = np.append(iterArray, np.append(dense5, str(activ[4])))
      iterArray = iterArray.reshape(5,3)
      for i in iterArray:
         for j in i:
            log.write(str(j) + '\t')
         log.write('\n')
      log.write(str(test_acc) + '\t' + str(epochsNum))
      log.write('\n')

log_keras()
