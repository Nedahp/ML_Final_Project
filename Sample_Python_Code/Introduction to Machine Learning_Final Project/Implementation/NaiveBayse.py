import SequenceSegmentsImport as ss
import numpy as np
from sklearn import metrics
from sklearn.cross_validation import StratifiedKFold
from sklearn import datasets
from sklearn import metrics
from sklearn.naive_bayes import GaussianNB
import collections
import operator

# Description: 
            # this function trains a Naive Bayse(NB) model based on the input data.
            # the function use cross validation to calculate the accuracy of the tranied model.
# input:
        # name of the input training data file
# output:
        # accuracy of the built NB model
def naieveBayes(fileName):
	# setup the data
	dataSet = ss.importSeqSegements(fileName)
	sequenceSegmentList, featuresNames, seqSegFeaturesValues, classLabel = ss.getAverageMatix(dataSet)
	seqSegFeaturesDiscreteVals = ss.descretizeAvgMatrix(seqSegFeaturesValues)
	classLabel = np.array(classLabel)
	
	# get indices for stratified-k-fold cross validation
	skf = StratifiedKFold(classLabel, n_folds=5)
    # calculate accuracy for each of the folds
	roc_auc = []
	for train_index, test_index in skf:
		
		# train on 4 folds and test on 1 remaining fold
		X_train, X_test = seqSegFeaturesDiscreteVals[train_index], seqSegFeaturesDiscreteVals[test_index]
		y_train, y_test = classLabel[train_index], classLabel[test_index]
			
		# fit a Naive Bayes model to the data
		model = GaussianNB()
		model.fit(X_train, y_train)
		
		# make predictions
		expected = y_test
		predicted = model.predict(X_test)
        
		# generate the ROCcurve
		false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(y_test, predicted)
        
        # calculate the area under ROCcurve
		roc_auc.append(metrics.auc(false_positive_rate, true_positive_rate))
        
	# final accuracy is the avrage accuracy of five folds	
	rocAucAvg = np.mean(roc_auc)	
	return rocAucAvg	



# Description: 
            # this function uses wrapper forward selection for feature selection
# input:
        # training data, class lables
# output:
        # selected features
def featureSelection(seqSegFeaturesDiscreteVals,classLabel):
    # wrapper forward selection
    (row, col) = seqSegFeaturesDiscreteVals.shape
    selectedFeatures = []
    acc = 0
    continuing = True
    # termination condition
    while (continuing):
        # initialize the data matrix with the data associated with features already selected
        data = np.zeros(shape=(row, col))
        for i in range(len(selectedFeatures)):
            for j in range(row):
                data[j,i] = seqSegFeaturesDiscreteVals[j,selectedFeatures[i]] 
        # delete all empty columns, save one for adding a new feature data values in future        
        emptyCols = range(len(selectedFeatures)+1,col)
        data = np.delete(data, emptyCols, 1)
        # select the next feature  
        accList = {}
        for i in range(col):
            if (i in selectedFeatures):
                # do nothing, the feature is already selected   
                continue
            else:
                # add a new column to data
                newColID = len(selectedFeatures)
                tempData = data
                for rowID in range(row):
                    tempData[rowID,newColID] = seqSegFeaturesDiscreteVals[rowID,i]
                # train a NB model on the current data
                model = GaussianNB()
                model.fit(tempData, classLabel)
                predicted = model.predict(tempData)
                confusionMatrix = metrics.confusion_matrix(classLabel, predicted)
                # compute the accuracy 
                currentAccuracy = float(confusionMatrix[0,0] + confusionMatrix[1,1])/float(np.sum(confusionMatrix))
                # store the accuracy resulted by trying each feature
                accList.update({i:currentAccuracy})   
                # identify what feature will generate the highest accuracy
        sorted_accList = sorted(accList.items(), key=operator.itemgetter(1))
        l = len(sorted_accList)
        t = sorted_accList.pop(l-1)
        maxAccuracy = t[1] 
        # check the termination condition
        if (maxAccuracy <= acc):
            continuing = False
        else:
            acc = maxAccuracy
            selectedFeatures.append(t[0])
            continuing = True
    # return the selected features        
    return selectedFeatures   


      
# Description: 
            # this function trains a Feature Selection Naive Bayse(FSNB) model based on the input data.
            # the function uses feature selection to select a subset of features which lead to the highest accuracy.
            # the function use cross validation to calculate the accuracy of the tranied model.
# input:
        # name of the input training data file
# output:
        # accuracy of the built FSNB model        
def featureSelection_naiveBayse(fileName):
    # setup the data
    dataSet = ss.importSeqSegements(fileName)
    sequenceSegmentList, featuresNames, seqSegFeaturesValues, classLabel = ss.getAverageMatix(dataSet)
    seqSegFeaturesDiscreteVals = ss.descretizeAvgMatrix(seqSegFeaturesValues)
    classLabel = np.array(classLabel)
    
    # feature selection
    selectedFeatures = featureSelection(seqSegFeaturesDiscreteVals,classLabel)
    
    # extract the data associated with the selected features
    (row, col) = seqSegFeaturesDiscreteVals.shape
    selectedFeatures_data = np.zeros(shape=(row, len(selectedFeatures)))
    for i in range(len(selectedFeatures)):
        for j in range(row):
            selectedFeatures_data[j,i] = seqSegFeaturesDiscreteVals[j,selectedFeatures[i]] 
            
    # get indices for stratified-k-fold cross validation
    skf = StratifiedKFold(classLabel, n_folds = 5)
    
    # calculate accuracy for each of the folds
    roc_auc = []
    for train_index, test_index in skf:
        
        # train on 4 folds and test on 1 remaining fold
        X_train, X_test = selectedFeatures_data[train_index], selectedFeatures_data[test_index]
        y_train, y_test = classLabel[train_index], classLabel[test_index]
        
        # fit a Naive Bayes model to the data
        model = GaussianNB()
        model.fit(X_train, y_train)
        
        # make predictions
        expected = y_test
        predicted = model.predict(X_test)
 
        # generate the ROCcurve
        false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(y_test, predicted)
        # calculate the area under ROCcurve
        roc_auc.append(metrics.auc(false_positive_rate, true_positive_rate))
        
    # final accuracy is the avrage accuracy of five folds
    rocAucAvg = np.mean(roc_auc)
    return rocAucAvg
    
    
def main():
    # call NB on all seven datasets
    NB_rocAucAvgList = []
    for i in range(1,7):
        fileName = "Set"+ str(i) + ".txt"
        NB_rocAucAvgList.append(naieveBayes(fileName))
    print NB_rocAucAvgList
    
    # call FSNB on all seven datasets
    FSNB_rocAucAvgList = []
    for i in range(1,7):
        fileName = "Set"+ str(i) + ".txt"
        FSNB_rocAucAvgList.append(featureSelection_naiveBayse(fileName))
    print FSNB_rocAucAvgList
    
if __name__ == '__main__':
    main()
