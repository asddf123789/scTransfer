import scipy
import scipy.io
import csv
import sklearn.preprocessing
import numpy as np

meta1 = []
with open('./data/Kidney_pdata.csv', newline='') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
     for row in spamreader:
         meta1.append(row)
meta2 = []
with open('./data/Kidney_rna_microwell_pdata.csv', newline='') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
     for row in spamreader:
         meta2.append(row)

meta1.pop(0)
meta2.pop(0)
label0 = []
label1 = []
label2 = []
for i in meta1:
    label1.append(i[-3].replace(" ", "_").replace("\"", ""))
    label0.append(i[-1].replace(" ", "_").replace("\"", ""))
for i in meta2:
    label2.append(i[8].split('_')[0].replace(" ", "_").replace("\"", ""))

enc = sklearn.preprocessing.LabelEncoder()
enc.fit(label2)
new_label2 = enc.transform(label2)
ohe = sklearn.preprocessing.OneHotEncoder()
ohe.fit(new_label2.reshape(-1,1))
ohe_label2 = ohe.transform(new_label2.reshape(-1,1)).toarray()


d1 = []
with open('./data/Kidney_atac_pc.csv', newline='') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
     for row in spamreader:
        row.pop(0)
        d1.append(row)

d2 = []
with open('./data/Kidney_rna_microwell_pc.csv', newline='') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
     for row in spamreader:
        row.pop(0)
        d2.append(row)

d1.pop(0)
d2.pop(0)
d1 = np.array(d1).T.astype(float)
d2 = np.array(d2).T.astype(float)

from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
# split data into train and test sets
seed = 8
test_size = 0.25
X_train, X_test, y_train, y_test = train_test_split(d2, new_label2, test_size=test_size, random_state=seed)
for i in range(5,20):
    print("K=",i)
    knn = KNeighborsClassifier(n_neighbors = i).fit(X_train, y_train)
    y_pred = knn.predict(X_test)
    predictions = [round(value) for value in y_pred]
    # evaluate predictions
    accuracy = accuracy_score(y_test, predictions)
    print("Validation Accuracy: %.2f%%" % (accuracy * 100.0))
    y_pred = knn.predict(X_train)
    predictions = [round(value) for value in y_pred]
    # evaluate predictions
    accuracy = accuracy_score(y_train, predictions)
    print("Training Accuracy: %.2f%%" % (accuracy * 100.0))

knn = KNeighborsClassifier(n_neighbors = 10).fit(d2, new_label2)
knn_preds = knn.predict(d1)

preds = list(knn_preds)
preds = enc.inverse_transform(preds)
preds = list(preds)
f = open('/Users/fanny/Documents/Shendure/mouse/mlproject/knnclassifier_rna_as_ref.csv', 'w')
f.write("\n".join(preds))
f.close()

