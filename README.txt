GENe RNA Catalogger
===================

The Biology Department of the College of William and Mary is conducting research using many samples of RNA from the species Xenopus Laevus, otherwise known as the African clawed frog. The research group is headed by Dr. Saha, and previously they were using a software tool known as BLAST provided by the National Center for Biotechnolgoy Information in order to catalog each of their unidentified RNA sequences. Each sequence of RNA needed to be submitted to this tool individually, and one of the researchers needed to identify the best candidate of the many that BLAST provides. After initially attempting to do this by hand with their 32,000 samples, Dr. Saha decided she would like this process automated. This will be my first extracurricular software development project as a rising sophomore at the College of William and Mary, and I hope to keep frequent updates.


Outline
=======
The program will need to perform three tasks;

1) Implement the decision making process that the biologists themselves would use for candidate selection. This will be achieved using a decision tree. All relevant information to the candidates chosen will be recored in a csv/excel file. This README should eventually contain a visual diagram in the form of a flow chart that represents the decison tree for the user's reference.

2) Remotely send requests to the BLAST servers for query. This will be acheived using the python library Biopython. This libary will also aid in the parsing of the returned candidates for data. 

3) Represent itself to the user as a comprehensive GUI, complete with loading bar for inevitable long run times, file selection gui, and column selection options - the program only needs a list of sequences to run, and if the user submits a large excel file with many columns where only one is a list of sequences, this program should be able to handle that-. The implementation of this portion of the project should be relatively simple, but it is listed here as a major feature because it will be the first time that I have created a GUI. I will use wxPython for this.

*4) Potential features depending on time: Toggle switches for certain features like setting statiscal significance e-values, number of candidates recorded, an entry bar to replace the default search words 'xenopus laevus' and 'xenopus tropicalus' with other search words, and options for the types of data to be recorded.  An ideal but difficult feature to implement would be to have some way for the users to define and change the existing data structure, in the case that their research interests change or another research team would like to use this software program

Version One
-----------
Version one will only have tasks one and two above implemented.



