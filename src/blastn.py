from Bio.Blast import NCBIWWW

result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
print(result_handle)