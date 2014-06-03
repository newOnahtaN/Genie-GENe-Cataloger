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

for description in blastRecord.descriptions:
	if "laevis" in description.title.lower() and first:
		print "Title: ", description.title, "\n   e-value: ", description.e , "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
		first = False
		count += 1
	if "tropicalis" in description.title.lower() and second:
		print "Title: ", description.title, "\n   e-value: ", description.e , "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
		second = False
		count += 1
	if "laevis" not in description.title.lower() and "tropicalis" not in description.title.lower():
		print "Title: ", description.title, "\n   e-value: ", description.e , "  Score: ", description.score, " Number of Alignments: ", description.num_alignments, "\n"
		count += 1
	if count == 3:
		break






