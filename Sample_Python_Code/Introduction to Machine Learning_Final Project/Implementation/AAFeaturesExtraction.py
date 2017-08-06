import urllib2

# This is a function that reads all the text on AAIndexURL page
# output:
	# featuresComponents: a list that contains an item for each feature along with 
	# 	all the related information for it
def readAAindexURL():

	AAIndexURL = "ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1"
	urlResponse = urllib2.urlopen(AAIndexURL)
	featuresComponents = urlResponse.read().split('//')

	return featuresComponents

	
# This function extracts the features Names, and the features values for each Amino Acid column
# output:
	# featuresNames: a list of the names of the features
	# featuresValues: a 2D list that containts values of each Amino Acid for every feature
	# 				 the #rows = to #features Names, and #columns = #Amino Acids
def extractFeaturesValues(featuresComponents):
	
	featuresNames = []
	featuresValues = []
	
	# each feature was saved as a full item in a list, they were originally split on (//)
	# the following split each item with respect to new lines and extracts features names and values
	featureItems = featuresComponents[0].splitlines()
	nameItems = featureItems[0].split()
	featuresNames.append(nameItems[1])
		
	values = featureItems[-2].split() + featureItems[-1].split()
	for i in range(0, len(values)):
		values[i] = float(values[i])
		
	featuresValues.append(values)
	featuresComponents.pop(0)
	del featuresComponents[-1]
	print "feature original count: ", (len(featuresComponents)+1)
	
	NACounter = 0
	for feature in featuresComponents:
		featureItems = feature.splitlines()
		
		nameItems = featureItems[1].split()
		values = featureItems[-2].split() + featureItems[-1].split()
	
		# any feature that has at least one NA as part of 
		# the Amino Acids values is ignored from the features list
		NAFlag = False
		for i in range(0, len(values)):
			if values[i] == 'NA':
				NAFlag = True
				NACounter +=1
				break
			else:
				values[i] = float(values[i])
			
		if not NAFlag:
			featuresNames.append(nameItems[1])
			featuresValues.append(values)
	
	# print featuresNames
	# print "The NA counter: ", NACounter
	# print "the final feature count: ", (len(featuresComponents) - NACounter +1)
	
	return featuresNames, featuresValues

# This function saves the featuresNames and featuresValues to 
# a txt file taking the form of a matrix
def exportDataToFile(featuresNames, featuresValues):
	dataFile = open("AAIndex_FeaturesData.txt", 'w')
	dataFile.truncate()
	
	for i, feature in enumerate(featuresNames):
		line = feature + " " 
		
		for value in featuresValues[i]:
			line = line + str(value) + " "
		
		line = line + "\n"	
		dataFile.write(line)
	dataFile.close()	

# This function exports the features data from a txt file and saves the data in two lists	
# output:
	# featuresNames: a list of the names of the features
	# featuresValues: a 2D list that containts values of each Amino Acid for every feature
	# 				 the #rows = to #features Names, and #columns = #Amino Acids
def importFeaturesData(fileName):
	dataFile = open(fileName, 'r')
	readDataLines = dataFile.read().splitlines()
	
	featuresNames = []
	featuresValues = []
	
	for line in readDataLines:
		lineItems = line.split()
		featuresNames.append(lineItems[0])
		
		lineItems.pop(0)
		values = []
		for i, item in enumerate(lineItems):
			values.append(float(item))
		
		
		featuresValues.append(values)	
		
	return featuresNames, featuresValues

	
def main():
	featuresComponents = readAAindexURL()
	featuresNames, featuresValues = extractFeaturesValues(featuresComponents)
	exportDataToFile(featuresNames, featuresValues)
	featuresNames, featuresValues = importFeaturesData("AAIndex_FeaturesData.txt")
	
if __name__ == "__main__":
	main()	
	