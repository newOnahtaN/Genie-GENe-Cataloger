from distutils.core import setup
import py2app

setup(
    app = ['GENe GUI.py'],
    data_files = ['README.txt'],
    options = {'py2app': {"iconfile" : "Lamp.icns"}  }
      )
