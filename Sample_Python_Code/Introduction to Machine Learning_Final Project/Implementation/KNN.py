import SequenceSegmentsImport as ss
import numpy as np
import math
from sklearn import metrics
from sklearn.cross_validation import StratifiedKFold
import collections
import operator
from sklearn.neighbors import KNeighborsClassifier

# Description: 
            # this function trains a K Nearest Neighbors(KNN) model based on the input data.
            # the function use cross validation to calculate the accuracy of the tranied model.
# input:
        # name of the input training data file
# output:
        # accuracy of the built KNN model
def KNN(fileName):
	
	# setup the data
	dataSet = ss.importSeqSegements(fileName)
	sequenceSegmentList, featuresNames, seqSegFeaturesValues, classLabel = ss.getAverageMatix(dataSet)
	classLabel = np.array(classLabel)
	
	# get indices for stratified-k-fold cross validation
	skf = StratifiedKFold(classLabel, n_folds = 5)
    
    # set the number of centroids
	k = 13
    
    # calculate accuracy for each of the folds
	roc_auc = []
	for train_index, test_index in skf:
	
		# divide the dataSet into training and testing sets according to the fold index
		X_train, X_test = seqSegFeaturesValues[train_index], seqSegFeaturesValues[test_index]
		y_train, y_test = classLabel[train_index], classLabel[test_index]
		
		# train KNN model
		kNN_model = KNeighborsClassifier(n_neighbors = k)
		kNN_model.fit(X_train,y_train)
		
        # use the trained model to predict lables for the test data
		predicted = kNN_model.predict(X_test)

        # use the area under the ROCcurve as a measurement for accuracy 
		false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(y_test, predicted)
		roc_auc.append(metrics.auc(false_positive_rate, true_positive_rate))
	
    # average resulted acuracy in k folds as the final accuracy
	rocAucAvg = np.mean(roc_auc)
	return rocAucAvg

	
def main():
    # call KNN on all seven datasets
	rocAucAvgList = []
	for i in range(1,7):
		fileName = "Set"+ str(i) + ".txt"
		rocAucAvgList.append(KNN(fileName))
	print rocAucAvgList
	
if __name__ == '__main__':
	main()	
		
	