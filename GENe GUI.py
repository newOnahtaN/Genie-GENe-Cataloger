import wx
import sys
import webbrowser
from GENe import *
from wx.lib.pubsub import setuparg1

class GENeGUI(wx.Frame):

	def __init__(self,parent,id):
		wx.Frame.__init__(self,parent,id,'GENe BLAST Automation and Gene Cataloger', size=(600,300), style=wx.MINIMIZE_BOX|wx.SYSTEM_MENU|wx.CAPTION|wx.CLOSE_BOX|wx.CLIP_CHILDREN)
		self.Center()
		self.panel = wx.Panel(self)

		#Check to see if there are other instances of this same program. If so, close out.
		self.name = "GENeGUI"
		self.instance = wx.SingleInstanceChecker(self.name)
		if self.instance.IsAnotherRunning():
			wx.MessageBox("Another instance is running already.", "Only one GENe instance at a time.")
			time.sleep(5)
			sys.exit()

		#IO Info
		self.openFileName = ''
		self.saveFileName = ''


		#Static Texts
		self.openFileText = wx.StaticText(self.panel, label = '', pos = (120, 16))
		self.saveFileText = wx.StaticText(self.panel, label = '', pos = (155, 57))
		self.colScroll = wx.StaticText(self.panel, label = 'Column Containing Sequences', pos = (421, 5))
		self.progressText = wx.StaticText(self.panel, label = "Waiting to run." , pos = (10, 225))
		self.localDBText = wx.StaticText(self.panel, label = "Local Database:" , pos = (60, 180))
		self.serverDBText = wx.StaticText(self.panel, label = "Server Database:" , pos = (290, 180))
		self.searchText = wx.StaticText(self.panel, label = "BLAST Search Type:" , pos = (290, 143))
		self.evalText = wx.StaticText(self.panel, label = "e-value maximum:" , pos = (60, 143))



		#Load Bar
		self.progressBar = wx.Gauge(self.panel, -1, 1, pos= (10,245), size= (400,20))



		#Buttons
		catalogButton = wx.Button(self.panel, label = "Run GENe", pos = (415,217), size = (136, 50))
		openFileButton = wx.Button(self.panel, label = "Open Excel File", pos = (10,10), size = (100, 30))
		saveFileButton = wx.Button(self.panel, label = "Choose Save Destination", pos = (10,50), size = (140,30))
		helpButton = wx.Button(self.panel, label = "Help", pos = (550,217), size = (40, 50))

		self.Bind(wx.EVT_BUTTON, self.saveExcelFile, saveFileButton)
		self.Bind(wx.EVT_BUTTON, self.openExcelFile, openFileButton)
		self.Bind(wx.EVT_BUTTON, self.runCataloger, catalogButton)
		self.Bind(wx.EVT_BUTTON, self.openREADME, helpButton)
		self.Bind(wx.EVT_CLOSE, self.closewindow)




		#Check Boxes
		self.localCheck = wx.CheckBox(self.panel, -1, "Local Blast Search (Requires BLAST+)", pos = (35,90))
		self.serverCheck = wx.CheckBox(self.panel, -1, "NCBI Server Blast Search (No installation required)", pos = (275,90))
		self.serverCheck.SetValue(True)
		self.queryType = 'queryServer'

		self.xenoCheck = wx.CheckBox(self.panel, -1, "Use Xenopus Laevis algorithm - When disabled, GENe simply records top three BLAST hits.", pos = (35, 115))
		self.xenoCheck.SetValue(True)
		self.xenopusAlgo = True

		self.Bind(wx.EVT_CHECKBOX, self.selectXeno, self.xenoCheck)
		self.Bind(wx.EVT_CHECKBOX, self.selectLocal, self.localCheck)
		self.Bind(wx.EVT_CHECKBOX, self.selectServer, self.serverCheck)




		#Choices
		self.dbChoice = wx.Choice(self.panel, -1, pos = (410, 178), choices = ['nr', 'refseq', 'swissprot', 'pat', 'month', 'pdb', 'refseq_mrna', 'refseq_genomic', 'est', 'est_human', \
		'est_mouse', 'est_others', 'gss', 'htgs', 'pat', 'dbsts', 'chromosome', 'genome', 'HTGS', 'RefSeq RNA', 'RefSeq Protein', 'Build RNA', 'Build Protein', 'ESTs'])
		self.dbChoice.SetStringSelection('nr')
		self.serverDatabase = 'nr'

		self.searchChoice = wx.Choice(self.panel, -1, pos = (410, 140), choices = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'])
		self.searchChoice.SetStringSelection('blastn')
		self.searchType = 'blastn'




		#Text Entry
		self.dbTextBox = wx.TextCtrl (self.panel, -1, 'nt', size = (30, -1), pos = (180, 178))
		self.localDatabase = 'nt'

		self.evalTextBox = wx.TextCtrl(self.panel, -1, '3', size = (30,-1), pos = (180, 140))
		self.eValCap = None



		#Scroll Counter Boxes
		self.seqColumn = wx.SpinCtrl(self.panel, value='0', pos = (480,24), size = (60, -1))
		self.seqColumn.SetRange(0,1000)
		self.seqCol = 0





		#Initialize GENe instance
		self.newCatalog = GENe(self)


	def runCataloger(self,event):
		"Main method. Passess in most of the GUI's instnace variables to GENe and runs GENe in a thread."

		#Make sure there is only one instance of GENe running as a result of this program.
		if self.newCatalog.running == False:

			#Collect contents from text boxes and GUI choices.
			self.localDatabase = self.dbTextBox.GetValue()
			self.eValCap = float(self.evalTextBox.GetValue())
			self.serverDatabase = self.dbChoice.GetStringSelection()
			self.searchType = self.searchChoice.GetStringSelection()
			self.seqCol = int(self.seqColumn.GetValue())


			#Set instance variables in GENe
			self.newCatalog.seqCol = self.seqCol
			self.newCatalog.eValCap = self.eValCap
			self.newCatalog.queryType = self.queryType
			self.newCatalog.searchType = self.searchType
			self.newCatalog.xenopusAlgo = self.xenopusAlgo
			self.newCatalog.localDatabase = self.localDatabase
			self.newCatalog.serverDatabase = self.serverDatabase


			#Force user to pick a file if they haven't already.
			if self.openFileName == '':
				self.openExcelFile(None)
				return

			if self.saveFileName == '':
				self.saveExcelFile(None)
				return


			#Run the GENe program
			self.newCatalog.readBook()
			self.progressBar.SetRange(self.newCatalog.numberOfQueries) 

			if self.queryType == 'queryServer':
				self.progressText.SetLabel("Initializing Server BLAST Query...")

			if self.queryType == 'queryLocal':
				self.progressText.SetLabel("Running Local BLAST Query. This may take a very long time.")

			self.newCatalog.start()



	def progressUpdate(self, progress):
		'''Updates progress bar'''

		if progress == -1:
			self.progressBar.Pulse()

		elif progress == -2:
			self.progressText.SetLabel("Complete!")
		
		else:
			self.progressBar.SetValue(progress)
			progressString = str(progress) + ' of ' + str(self.newCatalog.numberOfQueries) + " catalogged."
			self.progressText.SetLabel(progressString) 



	def openExcelFile(self, event):

		dialog = wx.FileDialog(self, message= "Open an Excel File", style = wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			self.openFileName = dialog.GetPath()

		if len(self.openFileName) > 45:
			shortName = self.openFileName[0:45] + '...'
		else:
			shortName = self.openFileName

		self.openFileText.SetLabel(shortName)
		self.newCatalog.fileToOpen = self.openFileName

		dialog.Destroy()



	def saveExcelFile(self,event):

		dialog = wx.FileDialog(self, message= "Choose save Destination", style = wx.SAVE)
		dialog.SetFilename('GENe Results.xls')

		if dialog.ShowModal() == wx.ID_OK:
			self.saveFileName = dialog.GetPath()

		#Make sure that the file being saved is in excel format.
		if '.xls' not in self.saveFileName:
			self.saveFileName = self.saveFileName + '.xls'

		#Shorten up the name to display on the Window
		if len(self.saveFileName) > 50:
			shortName = self.saveFileName[0:50] + '...'
		else:
			shortName = self.saveFileName

		self.saveFileText.SetLabel(shortName)
		self.newCatalog.saveAs = self.saveFileName

		dialog.Destroy()


	def selectLocal(self, event):
		'''Checkbox method for local selection'''
		if not self.localCheck.IsChecked():
			self.localCheck.SetValue(True)

		else:
			self.serverCheck.SetValue(False)
			self.queryType = 'queryLocal'


	def selectServer(self, event):
		'''Checkbox method for server selection'''
		if not self.serverCheck.IsChecked():
			self.serverCheck.SetValue(True)
		else:
			self.localCheck.SetValue(False)
			self.queryType = 'queryServer'


	def selectXeno(self, event):
		'''Checkbox method for Xenopus Laevis algorithm vs. Top three hits selection'''
		if self.xenoCheck.IsChecked():
			self.xenopusAlgo = True
		else:
			self.xenopusAlgo = False


	def openREADME(self, event):
		webbrowser.open('README.txt')


	def closewindow(self, event):
		try:
			self.newCatalog.exit()
			self.Destroy()
		except:
			self.Destroy()
			self.newCatalog.exit()


if __name__ == '__main__':
	app = wx.App(False) #filename = 'GENe Error Log.txt')
	frame = GENeGUI(parent=None, id=-1)
	frame.Show()
	app.MainLoop()
