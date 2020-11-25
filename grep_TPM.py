"""
Created on Nov 25, 2020
Author: Jiali
Usage: grep_TPM.py <input count gtf> <output tpm>
"""

import sys
file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as f, open(file_out,"w") as out:
    for line in f:
        if "TPM" in line:
            content = line.split("\t")
            gene_info = content[8].split(";")
            transcript_id = gene_info[1].replace('"','').replace(" transcript_id ","") # remove quotes and only keep transcript name
            tpm = gene_info[5].replace('"','').replace(" TPM ","")
            out.write(transcript_id + "\t" + tpm+"\n")