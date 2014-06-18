import sys
import time
import __future__
import pickle
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import *
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
		'''Instructions for use of this program: Most information about the purpose of this program can be found in it's readme.  If one wishes to use this program without the provided GUI, 
		then they must call three instructions from this class, an example of which is at the bottom of this program.  Firstly, create an instance of the GENe class, then call the readBook method on it
		(provided that it has the correct fileToOpen and saveAs instance variables), which returns a list of sequences. After that, call either localQuery or serverQuery on the class instance with 
		the list of sequences passed in as a parameter.  The program will then run the according BLAST search and write into the appropriate Excel file.'''


		#The excel file that this program saves to, in addition to a backup excel saved in this program's directory.
		self.saveAs = 'GENe Return.xls'

		#The excel file that this program reads from to gather sequences.
		self.fileToOpen = 'JustSequences.xlsx'#"SequencesVerySmall.xlsx"     

		#The database that the user would like to use according to this webpage: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#db
		#Use 'nr' as default because that is what Saha's algorithm should use
		self.serverDatabase = 'nr'

		#The database that the user would like to use according to this webpage: ftp://ftp.ncbi.nlm.nih.gov/blast/db/README
		#Use 'nt' as default because that is what Saha's algorithm should use
		self.localDatabase = 'nt'

		#The search type that the user would like to use. Available choices are blastn, blastp, blastx, tblastn, tblastx. This variable works for both Local and Server Blasts.
		self.searchType = 'blastn'

		#The column in the excel file provided just above that contains the sequences.  Remember to remind the user that the leftmost column starts with 0, not 1.
		self.seqCol = 0

		#The e-value cap/limit. Only applies to local blast.
		self.eValCap = 1

		#This Boolean decides whether to sort each hit by Dr. Saha's algorithm or to just return the top three results.
		self.topThree = False
		
		#Reference variables for Python Excel to create a new workbook.
		self.newBook = None
		self.sheet = None

		#Total number of queries. This variable isn't set until the readBook function is called.
		self.numberOfQueries = None


		self.openBook()


	def openBook(self):
		'''Opens an excel workbook to write each hit into. Creates headers.'''
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



	def readBook(self):
		'''Reads the provided excel workbook, identifies sequences and loads them into a list to be used by a query function.'''
		#make sure column actually has sequences

		try:
			workBook = open_workbook(self.fileToOpen)
		except IOError:
			raise IOError("File provided was not an Excel file. Please provide an excel file to be read from.")


		sequences = []
		for sheet in workBook.sheets():
			for row in range(sheet.nrows):
				sequences.append(sheet.cell(row,self.seqCol).value)


		if sequences == []:
			raise IOError('No sequences were loaded: The column was empty')

		self.numberOfQueries = len(sequences)

		return sequences


	def queryLocal(self,sequencesList):
		'''Given a python list of sequences, arbitrarily large, this method will call on BLAST+, a locally defined software program that emulates NCBI's BLAST alogrithm on a computer. This is a very consumptuous process, and 
		will more than likely lock up a computer, albeit a very fast option for query in comparison to the queryServer method. In order to use this option, the user must download BLAST+ from NCBI and the correct databases they
		wish to BLAST over.  The default for this program is the 'nt' database.'''

		'''#This bit of code splits up the list of sequences into portions of 50 as to not overload the computer when querying locally.  It also aids in displaying progress, rather than doing all at once.
		that = sequencesList[0:50]
		print that
		listOfSequenceLists = []
		start = 0
		end = 0
		while True:
			try:
				listOfSequenceLists.append(sequencesList[start:end])
				start += 50
				end += 50
			except:
				listOfSequenceLists.append(sequencesList[start:])
				break

		print listOfSequenceLists
		for sequences in sequencesList:'''



		'''sequenceFile = open("Sequences.fasta", "w")
		num = 0
		for sequence in sequences:
			sequenceFile.write('>')
			sequenceFile.write(str(num))
			sequenceFile.write('\n')
			sequenceFile.write(sequence)
			sequenceFile.write('\n')
			num += 1
		sequenceFile.close()

		timeOne = time.time()
		
		#blastn, blastp, blastx, tblastn, tblastx

		if self.searchType == 'blastn':
			blastCommandLine = NcbiblastnCommandline(cmd= self.searchType, query= "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out="localblastresults.xml")

		elif self.searchType == 'blastp':
			blastCommandLine = NcbiblastpCommandline(cmd= self.searchType, query= "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out="localblastresults.xml")

		elif self.searchType == 'blastx':
			blastCommandLine = NcbiblastxCommandline(cmd= self.searchType, query= "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out="localblastresults.xml")

		elif self.searchType == 'tblastn':
			blastCommandLine = NcbitblastnCommandline(cmd= self.searchType, query= "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out="localblastresults.xml")

		elif self.searchType == 'tblastx':
			blastCommandLine = NcbitblastxCommandline(cmd= self.searchType, query= "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out="localblastresults.xml")

		print blastCommandLine
		sys.stdout.flush()

		stdout, stderr = blastCommandLine()'''
		xmlFile = open('testxml.xml')

		blastRecords = NCBIXML.parse(xmlFile)


		#Row starts at 2 because rows 0 and 1 are filled by the header
		row = 2
		for blastRecord in blastRecords:
			self.filterNames(blastRecord, row)
			self.sheet.write(row,0,sequencesList[row-2])
			self.save()
			print '\n\n\n'
			sys.stdout.flush()
			row += 1
		
		self.stopTime(row,timeOne)
		self.save()

	def queryServer(self, sequences):
		'''Using each sequence gathered from the intial excel sheet, query each one by one to the NCBI Blast server. This method is very slow, but does not require that the user install BLAST+ or setup any
		corresponding databases locally on their computer.'''
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
					resultHandle = NCBIWWW.qblast(self.searchType, self.serverDatabase, sequence)
					if serverWasDown:
						print "Server is up and running again."
						sys.stdout.flush()
					break

				except:
					print "Server connection lost, waiting 60 seconds to try agiain.  Please make sure the computer has a working network connection."
					sys.stdout.flush()
					serverWasDown = True
					time.sleep(60)

			#for testing: blastRecord = pickle.load(open('recordDump.txt', 'rb'))
			blastRecord = NCBIXML.read(resultHandle)
			#for testing: pickle.dump(blastRecord, open('recordDump.txt', 'wb'))

			print "\n\n\n"
			sys.stdout.flush()

			self.filterNames(blastRecord,row)
			self.stopTime(row,timeOne)
			self.save()
			row += 1



	def filterNames(self,blastRecord, row):
		'''This contains Dr. Saha's simple algorithm for sorting through each sequence's hits. Collects three hits for each query. If the program is to be improved to include criteria for other organisms, if more data is to be collected, or
		if the decision making process for what should be recorded by this cataloger changes, this is the method that should be altered.  As it is, it filters through the descriptions for each hit to see if they contain the
		keywords 'laevis' or 'tropicalis'.  The user can also specify to ignore this alorithm and just return the top three results by setting self.topThree to True.'''

		if self.topThree == False:

			first = True
			second = True
			count = 0
			column = 0  
			print row-1, "of", self.numberOfQueries

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
					print 'Short Name: ', self.findShortName(description.title), '\n'

					self.writeData(row,column,description)
					column += 6
					second = False
					count += 1

				if "laevis" not in description.title.lower() and "tropicalis" not in description.title.lower():

					print "Title: ", self.cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", self.findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments
					print 'Short Name: ', self.findShortName(description.title), '\n'

					self.writeData(row,column,description)
					column += 6
					count += 1

				if count == 3:
					break

			if count == 0:
				self.sheet.write(row,1,"No Results")



		else:
			#Simply return top three results

			count = 0
			print row-1, "of", self.numberOfQueries

			for description in blastRecord.descriptions:
				
				print "Title: ", self.cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", self.findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments
				print 'Short Name: ', self.findShortName(description.title), '\n'

				self.writeData(row,column,description)
				column += 6
				first = False
				count += 1

				if count == 3:
					break

			if count == 0:
				self.sheet.write(row,1,"No Results")



	def save(self):
		'''Save the Excel Sheet that is being written into.'''
		while True:
			try: 
				self.newBook.save(self.saveAs)
				self.newBook.save('Backup.xls')
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
		'''Only to be used in a filter function, writes metadata from each BlastRecord instance into the excel sheet.'''
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



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
	newCatalog = GENe()
	sequences = newCatalog.readBook()
	newCatalog.queryLocal(sequences)






