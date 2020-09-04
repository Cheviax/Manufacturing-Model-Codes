from __future__ import absolute_import, division, print_function


import pandas as pd

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np
import seaborn as sns

#Path of the dataset to update
dataset_path = (r'/users/tlombardo/CGMD/PSO/PSO_ML1/ML_data_sets.csv', sep=',')

dataset = pd.read_csv(dataset_path, sep=',')
train_dataset = dataset.sample(frac=0.8,random_state=2)
test_dataset = dataset.drop(train_dataset.index)
X_train = train_dataset.iloc[:, -3:]
y_train = train_dataset.iloc[:, :7]

X_train_mean, X_train_std = X_train.mean(), X_train.std()
y_train_mean, y_train_std = y_train.mean(), y_train.std()


sns.pairplot(train_dataset[["Epsi_AM", "Sigma_AM", "Epsi_CBD", "Rho_CBD", "Kn",
                            "Gamman", "Xu", "slurry_density", "visco1(18)", "visco2(125)"]], diag_kind="kde")

def norm(x, mean, std):
  return (x - mean) / std
X_train_norm = norm(X_train, X_train_mean, X_train_std)
y_train_norm = norm(y_train, y_train_mean, y_train_std)
def norm1(x, x_min, x_max):
  return ( x - x_min ) / ( x_max - x_min )

def build_model():
  model = keras.Sequential([
    layers.Dense(64, activation=tf.nn.relu, input_shape=[len(X_train.keys())]),
    layers.Dense(64, activation=tf.nn.relu),
    layers.Dense(7)
  ])

  optimizer = tf.keras.optimizers.RMSprop(0.001)

  model.compile(loss='mean_squared_error',
                optimizer=optimizer,
                metrics=['mean_absolute_error', 'mean_squared_error'])
  return model
model = build_model()
model.summary()

class PrintDot(keras.callbacks.Callback):
  def on_epoch_end(self, epoch, logs):
    if epoch % 100 == 0: print('')
    print('.', end='')

EPOCHS = 600
full_batch = int(X_train.size/3)
print("Batch size = ", full_batch)

history = model.fit(X_train, y_train, epochs=EPOCHS, validation_split = 0.2, verbose=0, batch_size=full_batch)
hist = pd.DataFrame(history.history)
hist['epoch'] = history.epoch
hist.tail()

def plot_history(history):
    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch

X_test_norm = norm1(test_dataset.iloc[:, -3:],test_dataset.iloc[:, -3:].min(), test_dataset.iloc[:, -3:].max())
y_test_norm = norm1(test_dataset.iloc[:, :7],test_dataset.iloc[:, :7].min(), test_dataset.iloc[:, :7].max())
y_test = test_dataset.iloc[:, :7]
test_predictions = model.predict(test_dataset.iloc[:, -3:])
X_Wanted = np.array([[1, 1, 1]])
prediction = model.predict(X_Wanted)
print("Test set prediction :")
print(test_predictions)
test_labels = y_test_norm

print("Desired prediction : ")
print(prediction)
with open (r"/users/tlombardo/CGMD/PSO/PSO_ML1/predicted_out.txt", 'w') as f:
    for i in prediction:
        for j in i:
            f.write(str(j) + ",")
