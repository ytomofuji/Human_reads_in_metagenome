#!/usr/bin/env python

import sys
import re
args = sys.argv
INPUT=args[1]
OUTPUT=args[2]

result=['SNP\tCHROM\tPOS\tREF\tREF_CT\tALT\tALT_CT\tDP']

f=open(INPUT)
lines=f.readlines()

for i in range(len(lines)):
    METS=lines[i].rstrip('\n').split("\t")
    CTS=METS[7].split(",")
    ALLELES=[METS[5]]+METS[6].split(",")
    if(len(CTS)!=len(ALLELES)):
        print("ERROR")
        sys.exit()
    REF=METS[1]
    REF_CT=CTS[0]
    ALT=METS[2]
    if(ALT in ALLELES):
        ALT_CT=CTS[ALLELES.index(ALT)]
    else:
        ALT_CT=0
    result.append('\t'.join([METS[0],METS[3],METS[4],REF,str(REF_CT),ALT,str(ALT_CT),METS[8]]))

with open(OUTPUT,mode='w') as f:
    f.writelines('\n'.join(result)+'\n')