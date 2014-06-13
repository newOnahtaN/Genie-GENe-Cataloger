import time
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from xlrd import open_workbook
from xlwt import Workbook

#GENe RNA Cataloger
#Nathan Owen, ncowen@email.wm.edu, 757-752-4220
#Written for Dr. Saha of The College of William and Mary's Biology Department
#Github: https://github.com/newOnahtaN/Genie-GENe-Cataloger



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def findAcession(string):
	'''Given the standard string returned by the .title() method of the descriptions class from each blast record class, this function locates the acession number found between the second set of | | bars.'''
	barCount = 0
	returnString = ""
	for char in string:
		if char == "|":
			barCount += 1 
			continue
		if barCount == 4:
			break
		if barCount == 3:
			returnString = returnString + char

	return returnString

def cleanTitle(string):
	'''Given the standard string returned by the .title() method of the descriptions class from each blast record class, this function removes everything before the start of the actual title of the gene, indicated by the fourth | bar.'''
	barCount = 0
	returnString = ""
	index = 0
	for char in string:
		if char == "|":
			barCount += 1 
		index += 1
		if barCount == 4:
			break

	returnString = string[index:]
	return returnString

def writeData(row,column):
	sheetOne.write(row,column + 0, cleanTitle(description.title))
	sheetOne.write(row,column + 1, sequence)
	sheetOne.write(row,column + 2, description.e)
	sheetOne.write(row,column + 3, findAcession(description.title))
	sheetOne.write(row,column + 4, description.score)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#User input, file to open and the column that has the sequences, file to save to

#Make sure is excel file for both
fileToOpen = 'C:\Users\Nate\Desktop\GENe\SequencesVerySmall.xlsx'
saveAs = 'GENe Return Very Small.xls'

#The database that the user would like to use according to this webpage: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#db
#Use 'nr' as default
database = 'nr'

#make sure column actually has sequences
col = 0

#Open up an excel workbook, identify sequences and load them into a list to be queried

workBook = open_workbook(fileToOpen)
sequences = []

for sheet in workBook.sheets():
	for row in range(sheet.nrows):
		sequences.append(sheet.cell(row,col).value)


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create an excel document to write each hit into. Create headers.

newBook = Workbook()
sheetOne = newBook.add_sheet('Sheet One')
sheetOne.write(0,2,'First Hit')
sheetOne.write(0,8,'Second Hit')
sheetOne.write(0,14,'Third Hit')

x = 0
for y in range(3):
	sheetOne.write(1,x + 0,'Title')
	sheetOne.write(1,x + 1,'Sequence')
	sheetOne.write(1,x + 2,'e-value')
	sheetOne.write(1,x + 3,'Acession')
	sheetOne.write(1,x + 4,'Score')
	x += 6

sheetOne.write(1,18,'Load Time')
sheetOne.write(1,19,'Time of Day')


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using each sequence gathered from the intial excel sheet, query each one by one to the NCBI Blast server.

for sequence in sequences:

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Send of the sequence for query, and read it's return
	print "Working...."

	timeOne = time.time()

	serverWasDown = False
	while True: 

		try:
			#The function below is where all of the time is consumed. 
			resultHandle = NCBIWWW.qblast("blastn",database, sequence)
			if serverWasDown:
				print "Server is up and running again."
			break

		except:
			print "Server connection lost, waiting 60 seconds to try agiain.  Please make sure the computer has a working network connection."
			serverWasDown = True
			time.sleep(60)

	blastRecord = NCBIXML.read(resultHandle)

	print "\n\n\n"

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Dr. Saha's simple algorithm for sorting through each sequence's hits. Collects three hits for each query.
	first = True
	second = True
	count = 0
	column = 0 
	row = 2 
	for description in blastRecord.descriptions:

		if "laevis" in description.title.lower() and first:

			print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
			writeData(row,column)
			column += 6
			first = False
			count += 1

		if "tropicalis" in description.title.lower() and second:

			print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
			writeData(row,column)
			column += 6
			second = False
			count += 1

		if "laevis" not in description.title.lower() and "tropicalis" not in description.title.lower():

			print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
			writeData(row,column)
			column += 6
			count += 1

		if count == 3:
			break


	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Stop the timer, print it out and write it to the excel sheet. Save when done. 
	timeTwo = time.time()
	print 'Query Time: ', round(timeTwo - timeOne,1)
	sheetOne.write(row,18, round(timeTwo - timeOne,1))

	#Include time of day so analysis of when the best and worst times to use the server can be done.
	sheetOne.write(row,19, time.ctime())
	row += 1

	#Save the results into an excel spreadsheet, named by the input from the beginning of the program.
	while True:
		try: 
			newBook.save(saveAs)
			break
		except IOError:
			cont = raw_input("Please make sure the Excel spreadsheet is closed so that the program may continue. Press any key to continue after you have done so.")

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






