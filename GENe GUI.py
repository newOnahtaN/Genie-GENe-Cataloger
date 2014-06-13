import wx
from GENeObject import GENe

class GENeGUI(wx.Frame):

	def __init__(self,parent,id):
		wx.Frame.__init__(self,parent,id,'GENe RNA Cataloger', size=(485,300), style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER)
		panel = wx.Panel(self)
		catalogButton = wx.Button(panel, label = "Run Cataloger", pos = (10,10), size = (150, 50))
		self.Bind(wx.EVT_BUTTON, self.runCataloger, catalogButton)
		self.Bind(wx.EVT_CLOSE, self.closewindow)
		newCatalog = GENe()

	def runCataloger(self,event):
		newCatalog = GENe()
		sequences = newCatalog.readBook()
		newCatalog.queryServer(sequences)

	def closewindow(self, event):
		self.Destroy()

if __name__ == '__main__':
	app = wx.App(True)
	frame = GENeGUI(parent=None, id=-1)
	frame.Show()
	app.MainLoop()
