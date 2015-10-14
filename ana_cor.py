__author__ = 'sunxin'

#!/usr/bin/python

'''
This script analyse the result of correlation output.
And calculate additional correlation position.
Output is the modified snp information.
'''

from take_snp import COR
from corr import *
import sys

def read_cor (COR_FILE) :
    COR_FH = open(COR_FILE, 'r')

    COR_STORE = {}
    COR_ORDER = []

    while 1:
        COR_L = COR_FH.readline()

        if len(COR_L) == 0 :
            break

        cor_l = COR(COR_L)
        COR_STORE[cor_l.pos] = cor_l

        COR_ORDER.append(cor_l.pos)

    if COR_FH :
        COR_FH.close()

    return (COR_STORE, COR_ORDER)

def main() :

    print "take in "

    if len(sys.argv) == 1 :
        exit(0)

    COR_FILE = sys.argv[1]
    BAM_FILE = sys.argv[2]
    OUTPUT_FILE = sys.argv[3]

    global COR_STORE, COR_ORDER
    [COR_STORE, COR_ORDER] = read_cor(COR_FILE)


    MP_FH = read_mp(BAM_FILE)
    OUTPUT_FH = open(OUTPUT_FILE,'w')

    SNP_OUT = {}

    base_b1 = ""
    base_b2 = ""
    pos_b = 0
    for i in range(0,len(COR_ORDER)) :
        cor_i = COR_STORE[i]

        if cor_i.n_cor == 2:
            base_11 = cor_i.cor.keys()[0].split("/")[0]
            base_12 = cor_i.cor.keys()[0].split("/")[1]
            base_21 = cor_i.cor.keys()[1].split("/")[0]
            base_22 = cor_i.cor.keys()[1].split("/")[1]

            if base_11 != base_21 && base_12 != base_22 :
                if base_b1 == "" :
                    SNP_OUT[cor_i.pos.split(":")[0]] = [base_11, base_21]
                    base_b1 = base_12
                    base_b2 = base_22
                else :
                    if base_b1 == base_11 && base_b2 == base_21 :
                        SNP_OUT[cor_i.pos.split(":")[0]] = [base_b1, base_b2]
                        base_b1 = base_12
                        base_b2 = base_22

                    elif base_b1 == base_21 && base_b2 == base_11 :
                        SNP_OUT[cor_i.pos.split(":")[0]] = [base_11, base_21]
                        base_b1 = base_22
                        base_b2 = base_12
                    else :
                        ## mark for break in sequence
                        pos_b = cor_i.pos.split(":")[0]
                        SNP_OUT[cor_i.pos.split(":")[0]] = ["?"]
            elif 










