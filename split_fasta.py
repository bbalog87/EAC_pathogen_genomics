#!/usr/bin/env python3

from Bio import SeqIO

input_file = "output.fasta"

for record in SeqIO.parse(input_file, "fasta"):
    header = record.id
    filename = f"{header}"
    
    with open(filename, "w") as output_file:
        SeqIO.write(record, output_file, "fasta")

