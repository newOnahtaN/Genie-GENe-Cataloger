from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

sequenceSample = "CCTCCACAGAAGCTCATCACCACAATAGACTTCCTGGGTCCCCACATGCTGAATCCTCCTCCCTGAGAGTTCCTCCTAATGGGAATATTGCATCAGTGCTTCCAGTGGTTGCCTCTTCAAAACTATCTCCTCCGCTGCTCTCCTCCATGGCTTCTCTCTCCGCATTCCCCTTCTCTTTTGGATCTTTCCATCTATTGTCACCCAACTCTCTCAGCCCTACGACACCCACTCCATCAGGCAAGCCCTACAGGCCTTGGGGCACAGAAATCGGAGCCTTCTAAGTGAGAACTGATAAGCTTTTTTTTCGTGCAAGTAAAAAGCAGGATTGGGTGAGGGGGGTTGTAACGTCCAGCTGTGCTGGCCTTTTATTAATGACTTTTACATTGTATTTGCCAACCAGTGTGATCCAATGGTACCATTGCAAATTATTATTNGTCTTTNGTTTTTTTTAAACATAACATTTTTTTTCCAAAAAAAAAAAAATCTGGTACATTAGCAGCGTAACACAGAACTGTTCC"

print "Working...."

resultHandle = NCBIWWW.qblast("blastn","nr", sequenceSample)

blastRecord = NCBIXML.read(resultHandle)

print sequenceSample

print "\n\n\n"

maxE = .0

first = True
second = True
count = 0

def findAcession(string):
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

def cleanTitle(string):
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


for description in blastRecord.descriptions:
	if "laevis" in description.title.lower() and first:
		print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
		first = False
		count += 1
	if "tropicalis" in description.title.lower() and second:
		print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
		second = False
		count += 1
	if "laevis" not in description.title.lower() and "tropicalis" not in description.title.lower():
		print "Title: ", cleanTitle(description.title), "\ne-value: ", description.e , "Acession: ", findAcession(description.title), "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
		count += 1
	if count == 3:
		break






