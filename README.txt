GENe BLAST Automation and Gene Cataloger - June 2014
Nathan Owen, ncowen@email.wm.edu
========================================

This program was originally created for use by the Biology Department of the College of William and Mary. The GENe program was developed knowing the difficulty and tediousness of having to Blast many thousands of sequences and then having to sort through the hits of these Blasts by hand. The program confronts that challenge and is fully capable of handling as many sequences as necessary. GENe serves as a middleman between the user and NCBI's Blast servers or NCBI's local Blast tool, BLAST+, and accomplishes both of these tasks with the help of the opensource Biopython module. This program works first and foremost as a recipient of excel files that contain a list of sequences that need to be BLASTed. Users may choose to Blast these sequences either locally or on NCBI's servers one by one, but both options have drawbacks that will be discussed in more detail in this README.  


Version Compatability
=====================
Software for use:
-----------------
NCBI-BLAST+-2.2.29 

Software for development:
-------------------------
NCBI-BLAST+-2.2.29
Biopython, numPy
wxPython
Python Excel : wlwt , wlrd



Table of Contents
=================

1. Common Use and Basic Instructions
2. NCBI Server Blast
3. BLAST+ Local Blast
4. Xenopus Laevis Algorithm / Top three Blast hits
5. Backup Data and Error Checking
6. Closing the program (computer locking up)



1) Common Use and Basic Instructions
====================================

The first operations that are necessary for use of this program are found in the top left corner of the window. The 'Open Excel File' button has the user choose an excel workbook to read sequences out of.  It is possible that the user may choose a CSV file or variant that can be read by Excel and therefore has an Excel icon, but in order for GENe to read the file, it has to be saved as an Excel workbook. After the user selects an excel file to read from, the scroller to the right labeled 'Column Containing Sequences' should be adjusted so that the number showing is representative of the column in the excel spreadsheet that contains all of the sequences that are to be BLASTed.  It should be kept in mind that the leftmost column is considered 'Column 0' in this arrangement.

The user must then choose a directory that they wish to save the results to.  It is imperative that this file end with .xls - if the user does not type it at the end of the filename they choose, it will be added for them.  They must not, however, choose a filename that ends in some other extension than .xls or the program will error. If the user has mistakenly misnamed their save directory and did not realize until the end of a long calculation, a backup of their data can be retrieved. More information about this is in section 5 of this readme.

At this point, the user must decide whether to run a local or NCBI server Blast, whether to Use the Xenopus Laevis algorithm or not, a maximum e-value, the type of Blast they would like to perform, and over which database they would like to do it. The defaults for these settings are server Blast because it requires no installation other than this program, and the rest are set to accomplish the task of querying the specific set of sequences that this program was originally created for. 

The details of each of these settings besides the e-value maximum will be discussed in sections 2,3, and 4. 

When all the settings are set to the user's satisfaction, the Run GENe button should be pushed.  In most cases, the program should run appropriately until finished, but if the program errors or terminates early, the reasons why should be displayed on the python console. 

2) NCBI Server Blast
====================
This option is the defualt ticked because once this program is installed on a computer, this is the only option that is ready to run. Using this option is only advised for lists of sequences that are short or that did not return any results when Blasted over a local database.  This is because using the server blast is incredibly slow and processes each sequence one by one.  The process time of a sequence is entirely dependent on the status of NCBI's servers and whether or not other users are querying it at the same time.  During the evening and on weekends, an individual sequence query time can be as short as 5-15 seconds, but during peak hours it can reach as high as 1000 seconds or higher.  That being said, this option is definitely more reliable than running blast locally, and likely will crash less often. I highly advise that if a list of sequences is any longer than 500 entries, a local Blast be used.

The Blast search types for both local and server queries are the same - 'blastn, blastp, blastx, tblastn, tblastx'.  Information about these types of searches can be found here: http://blast.ncbi.nlm.nih.gov/Blast.cgi.  The database that these searches run over, however, differ greatly between server and local queries.  For this reason, only the database that is specified for the type of search that GENe is performing when you use it will be considered when it runs. The server databases are described in detail here: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#db . While many of the databases that appear on that page are represented in the GENe program, I cannot guarantee that all will work.  The only server database that I can be sure works is the 'nr' database. The ability to search over the rest of the databases relies upon whether or not BioPython supports them or not. Users can try their luck.

