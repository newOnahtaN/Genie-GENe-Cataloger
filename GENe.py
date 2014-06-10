import time
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from xlrd import open_workbook
from xlwt import Workbook

#GENe RNA Cataloger
#Nathan Owen, ncowen@email.wm.edu, 757-752-4220
#Written for Dr. Saha of The College of William and Mary's Biology Department
#Github: https://github.com/newOnahtaN/Genie-GENe-Cataloger



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




#Open up an excel workbook, identify sequences and load them into a list to be queried

workBook = open_workbook('C:\Users\Nate\Desktop\GENe\SequencesVerySmall.xlsx')

sequences = []

for sheet in workBook.sheets():
	for row in range(sheet.nrows):
		for col in range(sheet.ncols):
			sequences.append(sheet.cell(row,col).value)


#Create an excel document to write into with all end data

newBook = Workbook()
sheetOne = newBook.add_sheet('Sheet One')
sheetOne.write(0,2,'First Hit')
sheetOne.write(0,8,'Second Hit')
sheetOne.write(0,14,'Third Hit')

sheetOne.write(1,0,'Title')
sheetOne.write(1,1,'Sequence')
sheetOne.write(1,2,'e-value')
sheetOne.write(1,3,'Acession')
sheetOne.write(1,4,'Score')

sheetOne.write(1,6,'Title')
sheetOne.write(1,7,'Sequence')
sheetOne.write(1,8,'e-value')
sheetOne.write(1,9,'Acession')
sheetOne.write(1,10,'Score')

sheetOne.write(1,12,'Title')
sheetOne.write(1,13,'Sequence')
sheetOne.write(1,14,'e-value')
sheetOne.write(1,15,'Acession')
sheetOne.write(1,16,'Score')
row = 2 


for sequence in sequences:

	print "Working...."

	timeOne = time.time()
	resultHandle = NCBIWWW.qblast("blastn","nr", sequence)

	blastRecord = NCBIXML.read(resultHandle)

	print "\n\n\n"
	first = True
	second = True
	count = 0
	column = 0


	#Dr. Saha's simple algorithm, collect three data sets always. 
	for description in blastRecord.descriptions:
		if "laevis" in description.title.lower() and first:

			print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"

			sheetOne.write(row,column + 0, cleanTitle(description.title))
			sheetOne.write(row,column + 1, sequence)
			sheetOne.write(row,column + 2, description.e)
			sheetOne.write(row,column + 3, findAcession(description.title))
			sheetOne.write(row,column + 4, description.score)

			column += 6
			first = False
			count += 1

		if "tropicalis" in description.title.lower() and second:
			print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"

			sheetOne.write(row,column + 0, cleanTitle(description.title))
			sheetOne.write(row,column + 1, sequence)
			sheetOne.write(row,column + 2, description.e)
			sheetOne.write(row,column + 3, findAcession(description.title))
			sheetOne.write(row,column + 4, description.score)

			column += 6
			second = False
			count += 1

		if "laevis" not in description.title.lower() and "tropicalis" not in description.title.lower():
			print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"

			sheetOne.write(row,column + 0, cleanTitle(description.title))
			sheetOne.write(row,column + 1, sequence)
			sheetOne.write(row,column + 2, description.e)
			sheetOne.write(row,column + 3, findAcession(description.title))
			sheetOne.write(row,column + 4, description.score)

			column += 6
			count += 1

		if count == 3:
			break

	timeTwo = time.time()
	print 'Query Time: ', round(timeTwo - timeOne,1)
	row += 1

	newBook.save('GENe Return.xls')






