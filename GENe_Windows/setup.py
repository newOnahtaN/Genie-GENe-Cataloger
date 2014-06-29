from distutils.core import setup
import py2exe

setup(console=['GENe GUI.py'], data_files = ['README.txt'], 
	windows = [{"script" : "GENe GUI.py" , 
	            "icon_resources" : [(1, "GENe Icon.ico")]}])