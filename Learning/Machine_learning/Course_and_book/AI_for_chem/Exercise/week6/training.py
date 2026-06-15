from sklearn.model_selection import cross_validate  # 导入交叉验证
from sklearn.svm import SVC
from sklearn.metrics import classification_report, accuracy_score
import numpy as np

def train_the_model(model,X,y):
    """X,y is your training data; this function will train with cross validation and returns the best estimator"""
    cv_results = cross_validate(model, X, y, cv=5, 
                            return_train_score=True, 
                            return_estimator=True)
    
    # print mean score
    ts = cv_results['test_score']
    print("In train dataset, model accuracy is {:.4f} (+/- {:.4f}) ".format(np.mean(ts), 2* np.std(ts)))
    estimators = cv_results['estimator'] 

    best_estimator = max(estimators, key=lambda e: cv_results['test_score'][estimators.index(e)])
    return best_estimator

def test_the_model(model,X,y):
    """X,y is your test data; this function will train with do the prediction and print accuracy, returning the predicted value."""   
    # 预测测试集
    y_pred = model.predict(X)
    # 计算准确率
    accuracy = accuracy_score(y, y_pred)
    print(f"Accuracy: {accuracy:.3f}")
    # shou the report
    # classification_report_result = classification_report(y_test, y_pred, digits=3, zero_division=0)
    # print(classification_report_result)
    return y_pred