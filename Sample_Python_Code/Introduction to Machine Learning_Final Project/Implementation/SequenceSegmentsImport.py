import os
import AAFeaturesExtraction as aa
import numpy as np

# This function imports a list of sequence segments and their details from a txfile
# input: 
# 	fileName: the name of the text file containing the sequence dataFile
# output:
#	seqSegData: a list containing a sequence segment, protein ID, K residue position and class label	
def importSeqSegements(fileName):
	
	# get the file path of the data file
	filePath = "Sequence Segments Data" + os.sep + fileName
	
	seqSegFile = open(filePath,'r')
	seqSegFileLines = seqSegFile.read().splitlines()
	
	# get the index from which the sequence segments data starts
	startIdx = seqSegFileLines.index('Sequence Segment	Protein	Position	Central K Site (1: Ubiquitination; 0: Non_Ubiquitination)')
	
	seqSegData = []
	for i in range(startIdx+1, len(seqSegFileLines)):
		sequenceItems = seqSegFileLines[i].split()
		seqSegData.append(sequenceItems)	
	
	return seqSegData


# This function imports sequenece segment matrix, and generates segment PCP matrix
# input: 
# 	sequenece segment matrix
# output:
#	segment PCP matrix
def getAverageMatix(seqSegData):
	
	featuresNames, featuresValues = aa.importFeaturesData("AAIndex_FeaturesData.txt")
	n_features = len(featuresNames)
	featuresValues = np.array(featuresValues)
	
	aaIndexDict = {"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, "H":8 , "I":9, "L":10, "K":11, "M":12, "F":13, "P":14, "S":15, "T":16, "W":17, "Y":18, "V":19};	
	aminoAcids = aaIndexDict.keys()
	
	seqSegFraturesValues = np.empty([len(seqSegData),n_features])
	
	seqIdx = 0
	sequenceSegmentList = []
	classLabel = []
	
	for item in seqSegData:
		
		sequenceSegment = item[0]
		sequenceSegmentList.append(item[0])
		classLabel.append(int(item[-1]))
		
		aaFeaturesValues = np.zeros([len(sequenceSegment), n_features])
		# print aaFeaturesValues.shape
		index = 0
		
		for char in sequenceSegment:
			if char in aminoAcids:
				# print "current amino acid is %r and had index: %r" %(char, aaIndexDict[char])
				aaFeaturesValues[index] = np.transpose(featuresValues[:,aaIndexDict[char]])
				# print aaFeaturesValues
				# raw_input()
				index+=1
			else:
				aaFeaturesValues = np.delete(aaFeaturesValues, -1,0)
		seqSegFraturesValues[seqIdx] = np.mean(aaFeaturesValues, axis = 0)
		# print seqSegFraturesValues[seqIdx]
		seqIdx+=1
		# raw_input()	
		
		
	return sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel
	
    
# This function imports feature values from segment PCP matrix, and descretize them
# input: 
# 	feature values from segment PCP matrix
# output:
#	descretized feature values from segment PCP matrix
def descretizeAvgMatrix(seqSegFraturesValues):
	
	minValue = np.amin(seqSegFraturesValues, axis = 0)
	maxValue = np.amax(seqSegFraturesValues, axis = 0)
	valuesRange = np.ptp(seqSegFraturesValues, axis = 0)
	rangeIntervalCut = valuesRange/3
	
	(rows, cols) = seqSegFraturesValues.shape
	seqSegFraturesDiscreteVals = np.empty([rows, cols])
	
	for i in range(0,cols):
		seqIdx = 0
		
		for value in seqSegFraturesValues[:,i]:
			if minValue[i] <= value <= minValue[i]+rangeIntervalCut[i]:
				seqSegFraturesDiscreteVals[seqIdx,i] = 0
			elif minValue[i]+rangeIntervalCut[i] < value <= minValue[i]+(rangeIntervalCut[i]*2):
				seqSegFraturesDiscreteVals[seqIdx,i] = 1
			else:
				seqSegFraturesDiscreteVals[seqIdx,i] = 2
			seqIdx+=1	
		
	return seqSegFraturesDiscreteVals	
	

# This function writes segment PCP matrix to a file
# input: 
# 	file name and segment PCP matrix
# output:
#	no output
def writeAvgMatrixToFile(fileName, sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel):
	print len(sequenceSegmentList)
	print len(seqSegFraturesValues)
	
	dataFile = open(fileName,'w')
	dataFile.truncate()

	# writing a header line	
	line = "sequence Segment"
	for feature in featuresNames:
		line = line + "," + feature
	line = line + ", Class Label\n"
	dataFile.write(line)
	
	# writing the average matrix and the class label
	for seqIdx, sequenceSegment in enumerate(sequenceSegmentList):
		line = sequenceSegment
		
		for value in seqSegFraturesValues[seqIdx,:]:
			line = line + "," + str(value)
		
		line = line + "," + classLabel[seqIdx] + "\n"
		dataFile.write(line)
	
	dataFile.close()

	
def main():
	# combinedData = combineData()	
	# sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel = getAverageMatix(combinedData)
	# writeAvgMatrixToFile(sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel)
	
	dataSet1 = importSeqSegements("Set1.txt")
	print len(dataSet1)
	sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel = getAverageMatix(dataSet1)
	fileName = "AAIndex_AverageMatrix_DataSet1.csv"
	writeAvgMatrixToFile(fileName,sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel)
	
	dataSet2 = importSeqSegements("Set2.txt")
	print len(dataSet2)
	sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel = getAverageMatix(dataSet2)
	fileName = "AAIndex_AverageMatrix_DataSet2.csv"
	writeAvgMatrixToFile(fileName,sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel)
	
	dataSet3 = importSeqSegements("Set3.txt")
	print len(dataSet3)
	sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel = getAverageMatix(dataSet3)
	fileName = "AAIndex_AverageMatrix_DataSet3.csv"
	writeAvgMatrixToFile(fileName,sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel)
	
	dataSet4 = importSeqSegements("Set4.txt")
	print len(dataSet4)
	sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel = getAverageMatix(dataSet4)
	fileName = "AAIndex_AverageMatrix_DataSet4.csv"
	writeAvgMatrixToFile(fileName,sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel)
	
	dataSet5 = importSeqSegements("Set5.txt")
	print len(dataSet5)
	sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel = getAverageMatix(dataSet5)
	fileName = "AAIndex_AverageMatrix_DataSet5.csv"
	writeAvgMatrixToFile(fileName,sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel)
	
	dataSet6 = importSeqSegements("Set6.txt")
	print len(dataSet6)
	sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel = getAverageMatix(dataSet6)
	fileName = "AAIndex_AverageMatrix_DataSet6.csv"
	writeAvgMatrixToFile(fileName,sequenceSegmentList, featuresNames, seqSegFraturesValues, classLabel)
	
	
if __name__ == '__main__':
	main()	