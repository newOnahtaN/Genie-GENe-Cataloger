import wx
from GENe import GENe

class GENeGUI(wx.Frame):

	def __init__(self,parent,id):
		wx.Frame.__init__(self,parent,id,'GENe RNA Cataloger', size=(485,300), style=wx.MINIMIZE_BOX|wx.SYSTEM_MENU|wx.CAPTION|wx.CLOSE_BOX|wx.CLIP_CHILDREN)
		self.Center()
		self.panel = wx.Panel(self)
		self.openFileName = ''
		self.saveFileName = ''

		catalogButton = wx.Button(self.panel, label = "Run Cataloger", pos = (320,210), size = (150, 50))
		openFileButton = wx.Button(self.panel, label = "Open Excel File", pos = (10,10), size = (100, 30))
		saveFileButton = wx.Button(self.panel, label = "Choose Save Destination", pos = (10,50), size = (140,30))
		self.Bind(wx.EVT_BUTTON, self.saveExcelFile, saveFileButton)
		self.Bind(wx.EVT_BUTTON, self.openExcelFile, openFileButton)
		self.Bind(wx.EVT_BUTTON, self.runCataloger, catalogButton)
		self.Bind(wx.EVT_CLOSE, self.closewindow)


		self.searchType = None

		self.newCatalog = GENe()


	def runCataloger(self,event):


		newCatalog = GENe()
		sequences = newCatalog.readBook()
		newCatalog.queryLocal(sequences)


	def openExcelFile(self, event):
		dialog = wx.FileDialog(self, message= "Open an Excel File", style = wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			self.openFileName = dialog.GetPath()

		openFileText = wx.StaticText(self, label = self.openFileName, pos = (120, 16))
		self.newCatalog.fileToOpen = self.openFileName

		dialog.Destroy()


	def saveExcelFile(self,event):
		dialog = wx.FileDialog(self, message= "Choose save Destination", style = wx.SAVE)
		dialog.SetFilename('GENe Results.xls')
		if dialog.ShowModal() == wx.ID_OK:
			self.saveFileName = dialog.GetPath()

		saveFileText = wx.StaticText(self, label = self.saveFileName, pos = (155, 57))
		self.newCatalog.saveAs = self.saveFileName

		dialog.Destroy()


	def closewindow(self, event):
		self.Destroy()

if __name__ == '__main__':
	app = wx.App(True)
	frame = GENeGUI(parent=None, id=-1)
	frame.Show()
	app.MainLoop()
