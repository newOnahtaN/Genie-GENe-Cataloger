import os
import wx
import sys
import time
import shlex
import subprocess
import __future__
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import *
from Bio import SeqIO
from xlrd import open_workbook
from xlwt import Workbook
from threading import Thread

#GENe RNA Cataloger - Mac Version
#Nathan Owen, ncowen@email.wm.edu, 757-752-4220
#Written for Dr. Saha of The College of William and Mary's Biology Department
#Github: https://github.com/newOnahtaN/Genie-GENe-Cataloger

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

class GENe(Thread):


	def __init__(self, GUIWindow = wx.Window):
		'''Instructions for use of this program: Most information about the purpose of this program can be found in it's readme.  If one wishes to use this program without the provided GUI, 
		then they must call three instructions from this class, an example of which is at the bottom of this program.  Firstly, create an instance of the GENe class (wxPython must be installed), 
		then call the readBook method on it (provided that it has the correct fileToOpen and saveAs instance variables). After that, call either localQuery or serverQuery on the class instance.  
		The program will then run the according BLAST search and write into the appropriate Excel file.'''

		#The excel file that this program saves to, in addition to a backup excel saved in this program's directory.
		self.saveAs = 'C:\Users\Nate\Desktop\Overnight Local.xls'

		#The directory that all of the working files for this program saves to.
		self.saveDirectory = self.saveDirectory()

		#The excel file that this program reads from to gather sequences.
		self.fileToOpen = 'SequencesSmall.xlsx' 

		#The variable to be used by the .run() method when using the GUI that determines whether a to perform a Local or Server query. Default is server because it requires no set up.
		#Must be set to either 'queryServer' or 'queryLocal'. Unles you are using threads, this variable does not to be set correctly if you are going to directly call either the
		#queryServer or queryLocal methods.
		self.queryType = 'queryServer'

		#The database that the user would like to use according to this webpage: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#db
		#Use 'nr' as default because that is what Saha's algorithm should use
		self.serverDatabase = 'nr'

		#The database that the user would like to use according to this webpage: ftp://ftp.ncbi.nlm.nih.gov/blast/db/README
		#Use 'nt' as default because that is what Saha's algorithm should use
		self.localDatabase = 'nt'

		#The type of Blast search that the user would like to use. Available choices are blastn, blastp, blastx, tblastn, tblastx. This variable works for both Local and Server Blasts.
		self.searchType = 'blastn'

		#The column in the excel file provided just above that contains the sequences.  Remember to remind the user that the leftmost column starts with 0, not 1.
		self.seqCol = 0

		#The e-value cap/limit. Only applies to local blast.
		self.eValCap = 3

		#This Boolean decides whether to sort each hit by Dr. Saha's algorithm or to just return the top three results.
		self.xenopusAlgo = True
		
		#Reference variables for Python Excel to create a new workbook.
		self.newBook = None
		self.sheet = None

		#The list of sequences to be queried.  This variable is filled by the .readBook() method and must be called before anything can be queried.
		self.sequences = None

		#Total number of queries. This variable isn't set until the readBook method is called.
		self.numberOfQueries = None

		#An instance variable for use in the GUI - lets the GUI know if an instance of this program is already running.
		self.running = False

		#Sets up this class as a thread so that the GUI can update itself while this program is running.  self.GUI is what gets updated.
		Thread.__init__(self)
		self.GUI = GUIWindow



		self.openBook()


	def run(self):
		"""To be used with the GUI so that it may run this program in a seperate thread using the threading module. Accepts either the string 'queryLocal' or 'queryServer' 
		as well as a list of sequences.  Middleman to the queryLocal and queryServer methods."""
		if self.queryType == 'queryLocal':
			self.queryLocal()

		elif self.queryType == 'queryServer':
			self.queryServer()

		else:
			print "Please set either 'queryLocal' or 'queryServer' for self.queryType."


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
			self.sheet.write(1,x + 4,'Number of Alignments')
			self.sheet.write(1,x + 5,'Short Gene Name (May be incorrect)')
			x += 6

		self.sheet.write(1,20,'Load Time')
		self.sheet.write(1,21,'Time of Day')



	def readBook(self):
		'''Reads the provided excel workbook, identifies sequences and loads them into a list to be used by a query function.'''
		#make sure column actually has sequences
		try:
			workBook = open_workbook(self.fileToOpen)
		except:
			self.guiError("Could not open file. Please provide an excel file to be read from. If the file was a CSV file or some other variant, open it with Excel and save it as an Excel Worksheet.")
			raise IOError("Could not open file. Please provide an excel file to be read from. If the file was a CSV file or some other variant, open it with Excel and save it as an Excel Worksheet.")
			return


		sequences = []
		try:
			for sheet in workBook.sheets():
				for row in range(sheet.nrows):
					sequences.append(sheet.cell(row,self.seqCol).value)

		except:
			print 'Please choose a column that contains sequences. The leftmost column is zero.'
			self.guiError('Please choose a column that contains sequences. The leftmost column is zero.')
			return


		if sequences == []:
			self.guiError('No sequences were loaded: The column was empty')
			raise IOError('No sequences were loaded: The column was empty')
			return


		self.running = True
		self.numberOfQueries = len(sequences)
		self.sequences = sequences


	def queryLocal(self):
		'''Given a python list of sequences, arbitrarily large, this method will call on BLAST+, a locally defined software program that emulates NCBI's BLAST alogrithm on a computer. This is a very consumptuous process, and 
		will more than likely lock up a computer, albeit a very fast option for query in comparison to the queryServer method. In order to use this option, the user must download BLAST+ from NCBI and the correct databases they
		wish to BLAST over.  The default for this program is the 'nt' database.'''

		if self.sequences == None:
			raise IOError('You must first call the .readbook() method on the GENe instance before using this method so that sequences can be collected.')


		sequenceFile = open(self.saveDirectory + "Sequences.fasta", "w")
		num = 0
		for sequence in self.sequences:
			sequenceFile.write('>')
			sequenceFile.write(str(num))
			sequenceFile.write('\n')
			sequenceFile.write(sequence)
			sequenceFile.write('\n')
			num += 1
		sequenceFile.close()

		timeOne = time.time()
		
		#Make GUI loadbar pulse for a while because this next operation's progress can not be easily gauged. 
		self.loadBarUpdate(-1)


		try:

			#blastn, blastp, blastx, tblastn, tblastx
			self.guiError("Reminder: You may want to check the Error Log if this is your first time running this kind of a local search. If the local blast fails for some reason, GENe will not show it, it will just load forever. The error log will show if the local blast has failed.\n\n(This is not an error. If the error log shows no details of an error, then everything is working properly.)")
			#cmd= '/Users/TARDIS/BLAST/ncbi-blast-2.2.29+/bin/' + self.searchType,
			if self.searchType == 'blastn':
				blastCommandLine = NcbiblastnCommandline(query=self.saveDirectory + "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out=self.saveDirectory + "localblastresults.xml")

			elif self.searchType == 'blastp':
				blastCommandLine = NcbiblastpCommandline(query=self.saveDirectory + "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out=self.saveDirectory + "localblastresults.xml")

			elif self.searchType == 'blastx':
				blastCommandLine = NcbiblastxCommandline(query=self.saveDirectory + "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out=self.saveDirectory + "localblastresults.xml")

			elif self.searchType == 'tblastn':
				blastCommandLine = NcbitblastnCommandline(cmd= self.searchType, query=self.saveDirectory + "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out=self.saveDirectory + "localblastresults.xml")

			elif self.searchType == 'tblastx':
				blastCommandLine = NcbitblastxCommandline(query=self.saveDirectory + "Sequences.fasta", db= self.localDatabase , evalue=self.eValCap, outfmt=5, out=self.saveDirectory + "localblastresults.xml")

		except:

			self.loadBarUpdate(-3)
			self.guiError("Local Blast failed.")
			return



		print blastCommandLine		
		
		stdout, stderr = blastCommandLine()
		

		xmlFile = open(self.saveDirectory + 'localblastresults.xml')
		
	
		blastRecords = NCBIXML.parse(xmlFile)
		

		#Blast records is a list of blast records.  Each blast record corresponds to one sequence, and holds all the hits for each.
		#The following operations parse through those Blast Records.


		#Row starts at 2 because rows 0 and 1 are filled by the header
		row = 2
		for blastRecord in blastRecords:

			self.filterNames(blastRecord, row)

			#Write the sequence itself to the excel sheet
			self.sheet.write(row,0,self.sequences[row-2])

			if row % 500 == 0:
				self.sheet.flush_row_data()
				self.save()

			print '\n\n\n'
			

			row += 1
		

		self.stopTime(row,timeOne)
		self.save()
		self.running = False

		#Tell GUI loadbar that the program is done.
		self.loadBarUpdate(-2)



	def queryServer(self):
		'''Using each sequence gathered from the intial excel sheet, query each one by one to the NCBI Blast server. This method is very slow, but does not require that the user install BLAST+ or setup any
		corresponding databases locally on their computer.'''

		if self.sequences == None:
			raise IOError('You must first call the .readbook() method on the GENe instance before using this method so that sequences can be collected.')

		#Row starts at 2 because rows 0 and 1 are filled by the header
		row = 2
		for sequence in self.sequences:
			self.sheet.write(row,0,sequence)

			timeOne = time.time()

			print 'Working...'

			serverWasDown = False
			while True: 

				try:
					#The function below is where all of the time is consumed. 

					resultHandle = NCBIWWW.qblast(self.searchType, self.serverDatabase, sequence)
					if serverWasDown:
						print "Server is up and running again."
						self.guiError("Server is up and running again.")
					break

				except ValueError:
					print "The data provided could not be BLASTED. Please make sure the data that was provided was a column of sequences, that has no rows that are either empty or have anything other than sequences in them."
					self.guiError("The data provided could not be BLASTED. Please make sure the data that was provided was a column of sequences, that has no rows that are either empty or have anything other than sequences in them.")
					self.loadBarUpdate(-3)
					return

				except:
					print "Server connection lost, waiting 60 seconds to try agiain.  Please make sure the computer has a working network connection."
					self.guiError("Server connection lost, waiting 60 seconds to try agiain.  Please make sure the computer has a working network connection.")
					
					serverWasDown = True
					time.sleep(60)

			blastRecord = NCBIXML.read(resultHandle)

			print "\n\n\n"
			

			self.filterNames(blastRecord,row)
			self.stopTime(row,timeOne)
			self.save()
			row += 1

		self.running = False
		
		#Tell GUI loadbar that the program is done.
		self.loadBarUpdate(-2)


	def filterNames(self,blastRecord, row):
		'''This contains Dr. Saha's simple algorithm for sorting through each sequence's hits. Collects three hits for each query. If the program is to be improved to include criteria for other organisms, if more data is to be collected, or
		if the decision making process for what should be recorded by this cataloger changes, this is the method that should be altered.  As it is, it filters through the descriptions for each hit to see if they contain the
		keywords 'laevis' or 'tropicalis'.  The user can also specify to ignore this alorithm and just return the top three results by setting self.topThree to True.'''

		if self.xenopusAlgo == True:

			first = True
			second = True
			count = 0
			column = 0  
			print row-1, "of", self.numberOfQueries

			for description in blastRecord.descriptions:

				if "laevis" in description.title.lower() and first and float(description.e) < self.eValCap:

					print "Title: ", self.cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", self.findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments
					print 'Short Name: ', self.findShortName(description.title), '\n'

					self.writeData(row,column,description)
					column += 6
					first = False
					count += 1

				if "tropicalis" in description.title.lower() and second and float(description.e) < self.eValCap:

					print "Title: ", self.cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", self.findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments
					print 'Short Name: ', self.findShortName(description.title), '\n'

					self.writeData(row,column,description)
					column += 6
					second = False
					count += 1

				if "laevis" not in description.title.lower() and "tropicalis" not in description.title.lower() and float(description.e) < self.eValCap:

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
			column = 0
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


		self.loadBarUpdate(row - 1)
		


	def save(self):
		'''Save the Excel Sheet that is being written into.'''

		while True:

			try: 
				self.newBook.save(self.saveAs)
				self.newBook.save(self.saveDirectory + 'Backup.xls')
				break

			except IOError:

				print "Please make sure the Excel spreadsheet is closed so that the program may continue. Will wait 10 seconds and try again."
				self.guiError("Please make sure the Excel spreadsheet is closed so that the program may continue. Will wait 10 seconds and try again.")
				time.sleep(10)



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
		self.sheet.write(row,column + 4, description.num_alignments)
		self.sheet.write(row,column + 5, self.findShortName(description.title))



	def stopTime(self,row,timeOne):
		'''Stop the timer, print it out and write it to the excel sheet in the appropriate row.'''

		timeTwo = time.time()
		print 'Query Time: ', round(timeTwo - timeOne,1)
		self.sheet.write(row,20, round(timeTwo - timeOne,1))

		#Include time of day so analysis of when the best and worst times to use the server can be done.
		self.sheet.write(row,21, time.ctime())


	def saveDirectory(self):
		'''This method identifies or creates a directory to save the working files for this program to.'''
		if 'win32' in sys.platform:
			homeFolder = os.path.expanduser('~')
			GENefolder = homeFolder + '\GENe\\'
			if not os.path.exists(GENefolder):
				os.makedirs(GENefolder)
	
			return GENefolder
		else:
			return ''


	def exit(self):
		sys.exit()




	#GUI Methods

	def guiError(self, Errorstring):
		"Pass an error message to the GUI if there is one."
		try: 
			self.GUI.errorPop(Errorstring)
			
		except:
			pass

	def loadBarUpdate(self, progress):
		"Update the load bar of the GUI if there is one."
		try:
			self.GUI.progressUpdate(progress)
			
		except:
			pass
		






#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
	newCatalog = GENe()
	newCatalog.readBook()
	newCatalog.queryServer()






