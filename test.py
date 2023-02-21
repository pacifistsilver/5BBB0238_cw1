from Bio import SeqIO
from orffinder import orffinder

sequence = SeqIO.read("gene.fasta", "fasta")
orffinder.getORFs(sequence, minimum_length=75, remove_nested=True)