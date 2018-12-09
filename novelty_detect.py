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

meta1.pop(0)
label0 = []
label1 = []
for i in meta1:
    label1.append(i[-3].replace(" ", "_").replace("\"", ""))
    label0.append(i[-1].replace(" ", "_").replace("\"", ""))


d1 = []
with open('./data/Kidney_atac_pc.csv', newline='') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
     for row in spamreader:
        row.pop(0)
        d1.append(row)

d1.pop(0)
d1 = np.array(d1).T.astype(float)

enc = sklearn.preprocessing.LabelEncoder()
enc.fit(label1)
new_label1 = enc.transform(label1)
test_size = 0.2
seed = 8
for l in range(len(enc.classes_)):
    print(enc.classes_[l])
    X_outlier = d1[new_label1 == l]
    print(X_outlier.shape)
    X = d1[new_label1 != l]
    y = new_label1[new_label1 != l]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=seed)
    nu = 0.05
    gamma = 0.005

    clf = svm.OneClassSVM(nu=nu, kernel="rbf", gamma=gamma)
    clf.fit(X_train)
    y_pred_tr = clf.predict(X_train)
    y_test_pred = clf.predict(X_test)
    y_outlier_pred = clf.predict(X_outlier)
    n_error_tr = y_pred_tr[y_pred_tr == -1].size
    n_error_ts = y_test_pred[y_test_pred == -1].size
    n_error_outlier = y_outlier_pred[y_outlier_pred == 1].size
    print(n_error_tr, '/', y_pred_tr.size, '=', n_error_tr / y_pred_tr.size)
    print(n_error_ts, '/', y_test_pred.size, '=', n_error_ts / y_test_pred.size)
    print(n_error_outlier, '/', y_outlier_pred.size, '=', n_error_outlier / y_outlier_pred.size)


