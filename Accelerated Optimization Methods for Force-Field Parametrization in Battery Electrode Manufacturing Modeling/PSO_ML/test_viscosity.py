import numpy as np
import pandas as pd
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split

#Path of the dataset to update
train_data = pd.read_csv(r'/users/tlombardo/CGMD/PSO/PSO_ML1/ML_data_sets.csv', sep=',')
train_data.head()
print (train_data.shape)
X = train_data.iloc[:, -3:]
y = train_data.iloc[:, :7]
print (X.shape)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=15)
print(X_train.shape, X_test.shape)
nn = MLPRegressor(activation='relu', solver='sgd',hidden_layer_sizes=(5000), random_state=7, max_iter=500000, learning_rate='adaptive', tol=1e-3, alpha=0.0001)
print(nn.fit(X_train, y_train))
print(X_test.shape)
print(X_test)
with open(r'C:\Users\Teo Lombardo\Desktop\ML\Viscosity_Curve\results.txt', 'a') as f:
    liste = str(nn.predict(X_test)).split(']')

    new_liste = []
    for item in liste:
        item1 = ''.join(k for k in item if (k != '\n' and k != '[' ))
        new_liste.append(item1)
    liste = new_liste
    for item in liste:
        try:
            if item[0] == ' ':
                item1 = item[1:]
            else:
                item1 = item
        except:
            item1 = item
        f.write(str(item1) + '\n')
    f.close()
print(nn.predict(X_test))
X_wanted = pd.DataFrame(np.array([[1, 1, 1]]))
print(X_wanted.shape)
o = nn.predict(X_wanted)
print (o)
with open (r"C:\Users\Teo Lombardo\Desktop\test_DNN\predicted_out.txt", 'w') as f:
    for i in o:
        for j in i:
            f.write(str(j) + ",")

