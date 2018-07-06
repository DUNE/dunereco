from sklearn.metrics import classification_report, confusion_matrix
import numpy as np

test_info1 = np.load('/scratch/cvn/seresnet34/output/seresnet34_0.6lr.np').item()
test_info2 = np.load('/scratch/cvn/seresnet34/output/seresnet34_v3.np').item()

Y_pred1 = test_info1['Y_pred']
Y_pred2 = test_info2['Y_pred']

y_test1 = test_info1['y_test']
y_test2 = test_info2['y_test']

if np.array_equal(y_test1, y_test2):
    y_test = y_test1

else:
    print("different test sets!")
    exit(0)

Y_pred = (Y_pred1 + Y_pred2) / 2.0
#Y_pred = np.maximum(Y_pred1, Y_pred2)
#Y_pred = np.average(Y_pred1, Y_pred2)

y_pred  = np.argmax(Y_pred,  axis=1).reshape((Y_pred.shape[0],  1))  
y_pred1 = np.argmax(Y_pred1, axis=1).reshape((Y_pred1.shape[0], 1))  
y_pred2 = np.argmax(Y_pred2, axis=1).reshape((Y_pred2.shape[0], 1))  

inter_target_names = ['interac. 0', 'interac. 1', 'interac. 2', 'interac. 3', 'interac. 4', 'interac. 5', 
                         'interac. 6', 'interac. 7', 'interac. 8', 'interac. 9', 'interac. 10', 'interac. 11', 'interac. 13']

print('Classification report (interaction types):\n')

print(classification_report(y_test, y_pred, target_names=inter_target_names))
print(classification_report(y_test, y_pred1, target_names=inter_target_names))
print(classification_report(y_test, y_pred2, target_names=inter_target_names))

