from Bio import SeqIO

fastafile = "/group/bioi1/harrietd/git/STR-pipelines/simulate_reads/reference-data/hg19.STRdecoys.fasta"
outfasta = "hg19.STRdecoys.sorted.fasta"
allfasta = SeqIO.parse(fastafile, "fasta")
ids = sorted([rec.id for rec in allfasta])

record_index = SeqIO.index(fastafile, "fasta")
handle = open(outfasta, "wb")
for i in ids:
    handle.write(record_index.get_raw(i))
handle.close()
