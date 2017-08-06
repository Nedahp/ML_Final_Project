import SequenceSegmentsImport as ss
import numpy as np
from sklearn import metrics
from sklearn.cross_validation import StratifiedKFold
import collections
import operator
from sklearn import svm

# Description: 
            # this function trains a Support Vector Machine(SVM) model based on the input data.
            # the function use cross validation to calculate the accuracy of the tranied model.
# input:
        # name of the input training data file
# output:
        # accuracy of the built SVM model
def SVM(fileName):
	
	# setup the data
	dataSet = ss.importSeqSegements(fileName)
	sequenceSegmentList, featuresNames, seqSegFeaturesValues, classLabel = ss.getAverageMatix(dataSet)
	classLabel = np.array(classLabel)
	
	# get indices for stratified-k-fold cross validation
	skf = StratifiedKFold(classLabel, n_folds = 5)
	
    # calculate accuracy for each of the folds
	roc_auc = []
    
	for train_index, test_index in skf:
	
		# divide the dataSet into training and testing sets according to the fold index
		X_train, X_test = seqSegFeaturesValues[train_index], seqSegFeaturesValues[test_index]
		y_train, y_test = classLabel[train_index], classLabel[test_index]
		
		# train SVM model
		SVM_model = svm.SVC()
		SVM_model.fit(X_train,y_train)
		
		# use the trained model to predict lables for the test data
		predicted = SVM_model.predict(X_test)
        
        # use the area under the ROCcurve as a measurement for accuracy 
		false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(y_test, predicted)
		roc_auc.append(metrics.auc(false_positive_rate, true_positive_rate))
	
    # average resulted acuracy in k folds as the final accuracy
	rocAucAvg = np.mean(roc_auc)
	return rocAucAvg

	
def main():
    # call SVM on all seven datasets
	rocAucAvgList = []
	for i in range(1,7):
		fileName = "Set"+ str(i) + ".txt"
		rocAucAvgList.append(SVM(fileName))
	print rocAucAvgList
	
	
	
	
if __name__ == '__main__':
	main()	