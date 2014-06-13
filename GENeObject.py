import time
import __future__
import pickle
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from xlrd import open_workbook
from xlwt import Workbook

#GENe RNA Cataloger
#Nathan Owen, ncowen@email.wm.edu, 757-752-4220
#Written for Dr. Saha of The College of William and Mary's Biology Department
#Github: https://github.com/newOnahtaN/Genie-GENe-Cataloger

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

class GENe(object):


	def __init__(self):
		#User input, file to open and the column that has the sequences, file to save to.
		#Will make into input soon

		#Make sure is excel file for both
		self.saveAs = 'GENe Return Very Small.xls'

		#The database that the user would like to use according to this webpage: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#db
		#Use 'nr' as default because that is what this algorithm should use
		self.database = 'nr'
		self.newBook = None
		self.sheet = None
		self.openBook()



	def readBook(self):
		'''Open up an excel workbook, identify sequences and load them into a list to be queried'''
		#make sure column actually has sequences
		#make sure this is actually an excel file
		fileToOpen = 'C:\Users\Nate\Desktop\GENe\SequencesVerySmall.xlsx'
		col = 0
		workBook = open_workbook(fileToOpen)
		sequences = []

		for sheet in workBook.sheets():
			for row in range(sheet.nrows):
				sequences.append(sheet.cell(row,col).value)

		return sequences



	def openBook(self):
		'''Create an excel document to write each hit into. Create headers.'''
		self.newBook = Workbook()
		self.sheet = self.newBook.add_sheet('Sheet One')
		self.sheet.write(1,0,'Sequence')
		self.sheet.write(0,3,'First Hit')
		self.sheet.write(0,9,'Second Hit')
		self.sheet.write(0,15,'Third Hit')

		x = 0
		for y in range(3):
			self.sheet.write(1,x + 1,'Title')
			self.sheet.write(1,x + 2,'e-value')
			self.sheet.write(1,x + 3,'Acession')
			self.sheet.write(1,x + 4,'Score')
			self.sheet.write(1,x + 5,'Short Gene Name (May be incorrect)')
			x += 6

		self.sheet.write(1,20,'Load Time')
		self.sheet.write(1,21,'Time of Day')



	def save(self):
		"Save the Excel Sheet"
		while True:
			try: 
				self.newBook.save(self.saveAs)
				break
			except IOError:
				cont = raw_input("Please make sure the Excel spreadsheet is closed so that the program may continue. Press any key to continue after you have done so.")



	def findAcession(self,string):
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



	def cleanTitle(self,string):
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

	def findShortName(self,string):
		'''Given the standard string returned by the .title() method of the descriptions class from each blast record class, this function will attempt to find the shortened name of the gene in the name.  It will not always be succesful. 
		The current logic is that most shortened names exist within parantheses, and do not contain spaces. This is open to improvement.'''
		record = False
		returnString = ''
		for char in string:
			if char == ')':
				record = False
				if ' ' in returnString or 'ilurana' in returnString:
					returnString = ''
					continue
				else:
					break
			if record == True:
				returnString = returnString + char
			if char == '(':
				record = True

		return returnString



	def writeData(self,row,column,description):
		'''Only to be used in the filter function, writes metadata from each BlastRecord instance into the excel sheet.'''
		self.sheet.write(row,column + 1, self.cleanTitle(description.title))
		self.sheet.write(row,column + 2, description.e)
		self.sheet.write(row,column + 3, self.findAcession(description.title))
		self.sheet.write(row,column + 4, description.bits/description.score)
		self.sheet.write(row,column + 5, self.findShortName(description.title))



	def stopTime(self,row,timeOne):
		'''Stop the timer, print it out and write it to the excel sheet in the appropriate row.'''
		timeTwo = time.time()
		print 'Query Time: ', round(timeTwo - timeOne,1)
		self.sheet.write(row,20, round(timeTwo - timeOne,1))

		#Include time of day so analysis of when the best and worst times to use the server can be done.
		self.sheet.write(row,21, time.ctime())



	def queryLocal(self,sequences):

		sequenceFile = open("Sequences.fasta", "w")

		#Split up into multiple files, instead of just one that will make the stop time function make sense here. In order to do that, split up the 'sequences' python list, and when you do so, name each queried segment 'sequences' so 
		that the sequences can be easily written into the excel file in the for loop below.

		for sequence in sequences:
			sequenceFile.write('>')
			sequenceFile.write(str(num))
			sequenceFile.write('\n')
			sequenceFile.write(sheet.cell(row,col).value)
			sequenceFile.write('\n')
			num += 1
		sequenceFile.close()

		timeOne = time.time()
		
		#need to give user the option to change the database 
		blastnCommandLine = NcbiblastnCommandline(cmd= 'blastn', query= "Sequences.fasta", db= 'nt' , evalue=.05, outfmt=5, out="sequences.xml")
		print blastnCommandLine
		stdout, stderr = blastnCommandLine()
		xmlFile = open('sequences.xml')

		blastRecords = NCBIXML.parse(xmlFile)
		#for testing: blastRecords = NCBIXML.parse(open('Sequences.xml'))

		#Row starts at 2 because rows 0 and 1 are filled by the header
		row = 2
		for blastRecord in blastRecords:
			self.filterNames(blastRecord, row)
			self.sheet.write(row,0,sequences[row-2])
			print '\n\n\n'
			row += 1
		
		self.stopTime(row,timeOne)
		self.save()

		
		


	def queryServer(self, sequences):
		'''Using each sequence gathered from the intial excel sheet, query each one by one to the NCBI Blast server. Method takes a python list of sequences in regular text, unicode, or fasta format. 
		In order to call this function, the readBook() method needs to be called first in order to obtain the python list of sequences.'''
		#Row starts at 2 because rows 0 and 1 are filled by the header
		row = 2
		for sequence in sequences:
			self.sheet.write(row,0,sequence)

			print "Working...."
			timeOne = time.time()

			serverWasDown = False
			while True: 

				try:
					#The function below is where all of the time is consumed. 
					resultHandle = NCBIWWW.qblast("blastn", self.database, sequence)
					if serverWasDown:
						print "Server is up and running again."
					break

				except:
					print "Server connection lost, waiting 60 seconds to try agiain.  Please make sure the computer has a working network connection."
					serverWasDown = True
					time.sleep(60)

			#for testing: blastRecord = pickle.load(open('recordDump.txt', 'rb'))
			blastRecord = NCBIXML.read(resultHandle)
			#for testing: pickle.dump(blastRecord, open('recordDump.txt', 'wb'))

			print "\n\n\n"

			self.filterNames(blastRecord,row)
			self.stopTime(row,timeOne)
			self.save()
			row += 1



	def filterNames(self,blastRecord, row):
		'''Dr. Saha's simple algorithm for sorting through each sequence's hits. Collects three hits for each query. If the program is to be improved to include criteria for other organisms, if more data is to be collected, or
		if the decision making process for what should be recorded by this cataloger changes, this is the method that should be altered.  As it is, it filters through the descriptions for each hit to see if they contain the
		keywords 'laevis' or 'tropicalis'.'''
		first = True
		second = True
		count = 0
		column = 0  
		index = 0
		for description in blastRecord.descriptions:

			if "laevis" in description.title.lower() and first:

				print "Title: ", self.cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", self.findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments
				print 'Short Name: ', self.findShortName(description.title), '\n'
				self.writeData(row,column,description)
				column += 6
				first = False
				count += 1

			if "tropicalis" in description.title.lower() and second:

				print "Title: ", self.cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", self.findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments
				self.writeData(row,column,description)
				print 'Short Name: ', self.findShortName(description.title), '\n'
				column += 6
				second = False
				count += 1

			if "laevis" not in description.title.lower() and "tropicalis" not in description.title.lower():

				print "Title: ", self.cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", self.findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments
				self.writeData(row,column,description)
				print 'Short Name: ', self.findShortName(description.title), '\n'
				column += 6
				count += 1

			if count == 3:
				break

			index += 1



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
	newCatalog = GENe()
	sequences = newCatalog.readBook()
	newCatalog.queryLocal(sequences)






