GENe BLAST Automation and Gene Cataloger - June 2014
Nathan Owen, ncowen@email.wm.edu
Github and Source Files - https://github.com/newOnahtaN/Genie-GENe-Cataloger
========================================

This program was originally created for use by the Biology Department of the 
College of William and Mary. The GENe program was developed knowing the 
difficulty and tediousness of having to Blast many thousands of sequences and 
then having to sort through the hits of these Blasts by hand. The program 
confronts that challenge and is fully capable of handling as many sequences as 
necessary. GENe serves as a middleman between the user and NCBI's Blast servers 
or NCBI's local Blast tool, BLAST+, and accomplishes both of these tasks with 
the help of the opensource Biopython module. This program works first and 
foremost as a recipient of excel files that contain a list of sequences that 
need to be BLASTed. Users may choose to Blast these sequences either locally 
or on NCBI's servers one by one, but both options have drawbacks that will be 
discussed in more detail in this README.  


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



1) Common Use and Basic Instructions
====================================

Installation: Installation is as simple as downloading the appropriate installer 
for your operating system and running it. GENe will install itself according your 
operating system and you can then run the application directly.  If 
you are viewing this file from GitHub, you can either download this entire 
repository or click on either the Mac or Windows installer and click on ‘view raw’,
then you should run the GENe installer. If you have trouble opening GENe 
after installing it on Windows, attempt to run it as an administrator. Once you 
have done so, you can delete the rest of the files present if you are not 
planning on developing further on the source code. 

The first operations that are necessary for use of this program are found in the 
top left corner of the window. The 'Open Excel File' button has the user choose 
an excel workbook to read sequences out of.  It is possible that the user may 
choose a CSV file or variant that can be read by Excel and therefore has an 
Excel icon, but in order for GENe to read the file, it has to be saved as an 
Excel workbook. After the user selects an excel file to read from, the scroller 
to the right labeled 'Column Containing Sequences' should be adjusted so that 
the number showing is representative of the column in the excel spreadsheet that 
contains all of the sequences that are to be BLASTed.  It should be kept in mind 
that the leftmost column is considered 'Column 0' in this arrangement.

The user must then choose a directory that they wish to save the results to.  It 
is imperative that this file end with .xls - if the user does not type it at the 
end of the filename they choose, it will be added for them.  They must not, 
however, choose a filename that ends in some other extension than .xls or the 
program will error. If the user has mistakenly misnamed their save directory and 
did not realize until the end of a long calculation, a backup of their data can 
be retrieved. More information about this is in section 5 of this readme.

At this point, the user must decide whether to run a local or NCBI server Blast, 
whether to Use the Xenopus Laevis algorithm or not, a maximum e-value, the type 
of Blast they would like to perform, and over which database they would like to 
do it. The defaults for these settings are server Blast because it requires no 
installation other than this program, and the rest are set to accomplish the 
task of querying the specific set of sequences that this program was originally 
created for. 

The details of each of these settings besides the e-value maximum will be 
discussed in sections 2,3, and 4. 

When all the settings are set to the user's satisfaction, the Run GENe button 
should be pushed.  In most cases, the program should run appropriately until 
finished, but if the program errors or terminates early, the reasons why will 
be written into an error log that exists in a folder called 'GENe' that should 
be in your 'User' directory, also known as your Home Directory. 


The Excel file that is the result of this program is organized to show to the 
user three 'hits' for each sequence that were the results of the blast. How 
hese hits are chosen is discussed in section 4. The leftmost column is the 
sequence that was blasted, and each hit has it's name (or title) listed, the 
-value it scored, the acession number it has been assigned by NCBI, and, 
hopefully, its shortened gene name. Collecting the shortened gene name as 
metadata proved unsuccesful so I created a simple heuristic to scrape it from 
the full name of the gene. It is almost certain that only about 3/4 of hits will 
have a short gene name recorded, and among these it is likely that some are not 
entirely correct. The two rightmost columns are updated for each gene if a 
server Blast is being performed to give the user an idea of how long each query 
is taking, and if a local blast is being performed, then only the very last row 
will have this information and it will describe the lenth of the entire 
operation. 

If the user would like to check on their results before GENe has completed, they 
should copy the file that they choose to write the results to an open up the 
copy.  They can also view the file directly, but if they do so, GENe will pause 
until the file is closed again. 



2) NCBI Server Blast
====================

This option is the defualt ticked because once this program is installed on a 
computer, this is the only option that is ready to run. Using this option is 
only advised for lists of sequences that are short or that did not return any 
results when Blasted over a local database.  This is because using the server 
blast is incredibly slow and processes each sequence one by one.  The process 
time of a sequence is entirely dependent on the status of NCBI's servers and 
whether or not other users are querying it at the same time.  During the evening 
and on weekends, an individual sequence query time can be as short as 5-15 
seconds, but during peak hours it can reach as high as 1000 seconds or higher.  
That being said, this option is definitely more reliable than running blast 
locally, and likely will crash less often. On top of that, the results from this 
option are almost always more complete than the local blast option becauase the 
queries are run over a much greater portion of NCBI's databases by default, 
whereas the local blast's capability to yield a result for each sequence is 
dependent on how much of NCBI's database the user downloads to their machine. 
Even so, I highly advise that if a list of sequences is any longer than 500 
entries, a local Blast be used.

The Blast search types for both local and server queries are the same 
- 'blastn, blastp, blastx, tblastn, tblastx'.  Information about these types of 
searches can be found here: http://blast.ncbi.nlm.nih.gov/Blast.cgi.  The 
databases that these searches run over, however, differ greatly between server 
and local queries.  For this reason, only the database that is specified for the 
type of search that GENe is performing (server or local) will be considered when 
it runs. The server databases are described in detail here: 
http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#db . 
While many of the databases that appear on that page are represented in the GENe 
program, I cannot guarantee that all will work.  The only server database that I 
can be sure works is the 'nr' database. The ability to search over the rest of 
the databases relies upon whether or not BioPython supports them or not. Most 
should, but there are likely a few that do not work.



3) BLAST+ Local Blast
=====================

Blasting sequences locally, while requiring some set up, pays off exponentially 
(quite literally) as the list of sequences that need to be Blasted grows. Local 
blast, unlike the server blast, process all sequences in a batch request, and is 
much faster than the server option for this reason. The main drawback is that 
while the server query option can quietly run in the background of a system 
sending of requests slowly but surely, BLAST+ is very much a brute force method 
and will consume all of the system's memory and processing power until finished. 
This means that while this program much faster, it leaves the system unusable 
for other purposes until it has finished. That cost is well worth the payoff 
however as a sample of 32,000 sequences can be BLASTed locally in just over 24 
hours whereas 32,000 sequences with the server option may take weeks, even months.
Depending on the databases selected for local blast, it is likely that percentage
of sequences will return no results.  It is in this case that I advise that those
sequences that did not return anything from a local blast then be blasted on the 
server, in order to obtain complete data. 

VERY IMPORTATNT: If the user wishes to close out of GENe while a local Blast is 
executing in order to use their system, closing GENe WILL NOT close out the 
local blast process. If GENe terminates before the local blast does, the user needs
to open their task manager and local the blast process manually in order to end it. 

Much of the information about the BLAST+ Local Command Line Applications can be 
found at this web address: http://www.ncbi.nlm.nih.gov/books/NBK1763/ , but I wi
ll describe here only what needs to be set up in order to use GENe. The BLAST+ 
suite is a set of tools that allow users to perform regular blast queries using 
their operating system's command line.  For users that are uncomfortable using 
the command line, GENe is a great GUI alternative.  

In order to use the BLAST+ option, the user must install it from the webpage in
the paragraph above. There are instructions that are specific to each operating 
system, and they must be followed carefully. When BLAST+ is installed, settings 
known as environment variables must be appropriately adjusted so that the 
perating system's command line can access the BLAST+ executables and also the 
databases that it is supposed to run over. Specific instructions on how to do 
this are again found in the link above in the installation section.

In addition to having BLAST+ installed, a set of databases must be downloaded an 
kept up to date, or custom databases may be created. This webpage 
- ftp://ftp.ncbi.nlm.nih.gov/blast/db/ - contains all of NCBI's most current 
databases. Descriptions of each database on that page are found here: 
ftp://ftp.ncbi.nlm.nih.gov/blast/db/README. Each databases is split up into many 
.##.tar.gz files that are each about a gigabyte large. The prefix before the 
.##.tar.gz extenstion is the name of the database, and in order for that database 
to be used, all of the files from 0 to the topmost ## must be downloaded and 
extracted using an extracting tool like WinZip. The folder that contains the 
extracted portions of the database(s) you have chosen needs to have a path 
ariable set to it so that the command line can directly access it. On windows, 
a new system variable needs to be created named BLASTDB that directs to this 
folder's directory. In MacOS and Unix systems, instructions on how to do this 
can be found in the installation instructions referenced earlier.

When a database has either been downloaded or a custom database has been set up 
and BLAST+ is fully configured, GENe is ready to run.  Before running, the user 
must type their database(s) into the 'Local Database:' box. If the user is only 
using one database, then it shoud be the name of that database as it appears
before the .##.tar.gz extension. For example, if you want to use the nt database,
just write nt in the box.  Same thing for any other database. If the user would 
like to use multiple databases, then it is as simple as writing each one out with 
a single space between each, no commas necessary. For example, in order 
to blast the nr, sprot, and trembl databases at the same time, one would 
write - nr sprot trembl - in the box (no hyphens). Keep in mind that when blasting
multiple databases, all databases must be of the same molecule type or BLAST+ will error. 

If a user creates a custom database, instructions on how to blast it can be 
found in the manual above. It should be as simple as just writing out the entire path
directory to wherever you have saved the database, ie. C:\Local Disk\.......\Database

4) Xenopus Laevis Algorithm / Top Three Blast Hits
==================================================

As stated before, this program was originally written for a team of researchers 
at the College of William and Mary's Biology Department whose interest was in 
cataloging several thousands sequences from the organism Xenopus Laevis, the 
African Clawed Frog. The original intention for this program that did not get 
implemented to due to time constraints was to allow the user to customize their 
own algortihm for sorting through the hits of a BLAST.  A specialized algorithm 
was instead developed for this team of researchers, and it works as follows:




Always record three hits:

If there exists a hit that is from the organism Xenopus Laevis, let it be 
recorded first. Do not record any more of this type.

If there exists a hit that is from the organism Xenopus Tropicalis, let it be 
recorded second. Do not record any more of this type.

After searching to find hits that are either from Xenopus Laevis or Tropicalis, 
fill any blank spaces with the rest of the hits ordered by lowest e-value first.




This algorithm will always record at least three hits unless there are less than
three hits in the first place.  It gives priority to hits from the expected 
organisms, but also accounts for off-the-wall hits that have low e-values. 


A user or developer has two alternatives to using this algorithm: they may 
either simply uncheck the box that is labeled "Use Xenopus Laevis algorithm" to 
just record the top three hits, or they may access this programs source code and 
construct a new algorithm.  I have left the GENe code wide open for this specific 
development - all that one would need to do is to local the .filterNames() 
method in the GENe class located in GENe.py, and using the documentation from 
Biopython that describes their BlastRecord class (found here in chapter 7.3 
http://biopython.org/DIST/docs/tutorial/Tutorial.html) create a new sorting 
algorithm. If you would like me to personally create a new sorting algorithm for 
you, please feel free to contact me.  I can't guarantee that I'll have the time 
to help, but I will be eager to help if I do.



5) Backup Data and Error Checking
=================================

If an excel file is lost, there will always be a backup of the most recent GENe 
operation performed saved as an excel file that will be present in the GENe folder 
which is located in your 'User' directory, also known as your Home Directory. 
This is only the case if you are on a Windows machine unfortunately because the
Mac will not allow programs like this one to create directories. Instead, the Mac
version saves a copy to this program’s working directory, which is very difficult 
to access.

If the program terminates early, or errors for any reason, the error message is 
written into a file called 'GENe Error Log.txt' that can also be found in the 
same directory as the backup file if you are on Windows, and if you are on MacOS 
then you can open it directly from the GUI.  If GENe errors, the error log will 
not update until you close GENe, so when it does error on Mac, close it down, 
open it up again, and then check the error log button If an error persists and 
cannot be resolved, please feel free to contact me at ncowen@emailwm.edu
