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

    print "take in COR_FILE, BAM_FILE, SNP_FILE, OUTPUT_FILE"

    if len(sys.argv) == 1 :
        exit(0)

    COR_FILE = sys.argv[1]
    BAM_FILE = sys.argv[2]
    SNP_FILE = sys.argv[3]
    OUTPUT_FILE = sys.argv[4]

    global COR_STORE, COR_ORDER, MP_FH
    [COR_STORE, COR_ORDER] = read_cor(COR_FILE)

    global SNP_STORE, SNP_ORDER

    [SNP_STORE,SNP_ORDER] = read_snp(SNP_FILE)


    MP_FH = read_mp(BAM_FILE)
    OUTPUT_FH = open(OUTPUT_FILE, 'w')
    MCOR_FH = open("mod_" + COR_FILE, 'w')

    SNP_OUT = {}

    base_b1 = ""
    base_b2 = ""
    pos_b = -10
    for i in range(0,len(COR_ORDER)) :
        cor_i = COR_STORE[COR_ORDER[i]]

        if cor_i.n_cor == 2:
            base_11 = cor_i.cor.keys()[0].split("/")[0]
            base_12 = cor_i.cor.keys()[0].split("/")[1]
            base_21 = cor_i.cor.keys()[1].split("/")[0]
            base_22 = cor_i.cor.keys()[1].split("/")[1]

            if base_b1 == "" :
                    SNP_OUT[cor_i.pos.split(":")[0]] = [base_11, base_21]
                    base_b1 = base_12
                    base_b2 = base_22

                    if pos_b != -10 :
                        ## output cor relation between two position
                        [corr_o_hash,corr_o_sum] = cor_snp(str(pos_b), str(int(float(cor_i.pos.split(":")[0]))))
                        corr_o_line = str(i+1) + ":" + str(i+2) + '\t'

                        for i in corr_o_hash.keys() :
                            if corr_o_hash[i] >= 0.1 :
                                corr_o_line = corr_o_line + i + '\t' + str(corr_o_hash[i]) + '\t'

                        corr_o_line += str(corr_o_sum)

                        print >> MCOR_FH, corr_o_line

                        pos_b = -10


            elif base_11 != base_21 and base_12 != base_22 :
                if base_b1 == base_11 and base_b2 == base_21 :
                    SNP_OUT[cor_i.pos.split(":")[0]] = [base_b1, base_b2]
                    base_b1 = base_12
                    base_b2 = base_22

                    if pos_b != -10 :
                        ## output cor relation between two position
                        [corr_o_hash,corr_o_sum] = cor_snp(str(pos_b), str(int(float(cor_i.pos.split(":")[0]))))
                        corr_o_line = str(i+1) + ":" + str(i+2) + '\t'

                        for i in corr_o_hash.keys() :
                            if corr_o_hash[i] >= 0.1 :
                                corr_o_line = corr_o_line + i + '\t' + str(corr_o_hash[i]) + '\t'

                        corr_o_line += str(corr_o_sum)

                        print >> MCOR_FH, corr_o_line

                        pos_b = -10


                elif base_b1 == base_21 and base_b2 == base_11 :
                    SNP_OUT[cor_i.pos.split(":")[0]] = [base_b1, base_b2]
                    base_b1 = base_22
                    base_b2 = base_12

                    if pos_b != -10 :
                        ## output cor relation between two position
                        [corr_o_hash,corr_o_sum] = cor_snp(str(pos_b), str(int(float(cor_i.pos.split(":")[0]))))
                        corr_o_line = str(i+1) + ":" + str(i+2) + '\t'

                        for i in corr_o_hash.keys() :
                            if corr_o_hash[i] >= 0.1 :
                                corr_o_line = corr_o_line + i + '\t' + str(corr_o_hash[i]) + '\t'

                        corr_o_line += str(corr_o_sum)

                        print >> MCOR_FH, corr_o_line

                        pos_b = -10

                else :
                    SNP_OUT[cor_i.pos.split(":")[0]] = ["?"]
                    base_b1 = base_12
                    base_b2 = base_22
                    ## mark for break in sequence
                    if pos_b == -10 :
                        pos_b = int(float(COR_STORE[COR_ORDER[i - 1]].pos.split(":")[0]))



            elif base_11 == base_21 :

                SNP_OUT[cor_i.pos.split(":")[0]] = [base_b1,base_b2]
                base_b1 = base_12
                base_b2 = base_22

                if pos_b == -10 :
                    pos_b = int(float(COR_STORE[COR_ORDER[i - 1]].pos.split(":")[0]))


            else :
                if base_b1 == base_11 and base_b2 == base_21 :
                        SNP_OUT[cor_i.pos.split(":")[0]] = [base_b1, base_b2]
                        base_b1 = base_12
                        base_b2 = base_22

                        if pos_b != -10 :
                            ## output cor relation between two position
                            [corr_o_hash,corr_o_sum] = cor_snp(str(pos_b), str(int(float(cor_i.pos.split(":")[0]))))
                            corr_o_line = str(i+1) + ":" + str(i+2) + '\t'

                            for i in corr_o_hash.keys() :
                                if corr_o_hash[i] >= 0.1 :
                                    corr_o_line = corr_o_line + i + '\t' + str(corr_o_hash[i]) + '\t'

                            corr_o_line += str(corr_o_sum)

                            print >> MCOR_FH, corr_o_line

                            pos_b = -10


                elif base_b1 == base_21 and base_b2 == base_11 :
                        SNP_OUT[cor_i.pos.split(":")[0]] = [base_b1, base_b2]
                        base_b1 = base_22
                        base_b2 = base_12

                        if pos_b != -10 :
                            ## output cor relation between two position
                            [corr_o_hash,corr_o_sum] = cor_snp(str(pos_b), str(int(float(cor_i.pos.split(":")[0]))))
                            corr_o_line = str(i+1) + ":" + str(i+2) + '\t'

                            for i in corr_o_hash.keys() :
                                if corr_o_hash[i] >= 0.1 :
                                    corr_o_line = corr_o_line + i + '\t' + str(corr_o_hash[i]) + '\t'

                            corr_o_line += str(corr_o_sum)

                            print >> MCOR_FH, corr_o_line

                            pos_b = -10



                else :
                    SNP_OUT[cor_i.pos.split(":")[0]] = [base_b1,base_b2]
                    base_b1 = base_12
                    base_b2 = base_22

                    if pos_b == -10 :
                        pos_b = int(float(COR_STORE[COR_ORDER[i - 1]].pos.split(":")[0]))

        elif cor_i.n_cor == 1:
            base_11 = cor_i.cor.keys()[0].split("/")[0]
            base_12 = cor_i.cor.keys()[0].split("/")[1]
            SNP_OUT[cor_i.pos.split(":")[0]] = [base_b1,base_b2]
            base_b1 = base_12
            base_b2 = base_12

            if pos_b == -10 :
                pos_b = int(float(COR_STORE[COR_ORDER[i - 1]].pos.split(":")[0]))

        else :
            ## cor_i.n_cor == 3 or else
            SNP_OUT[cor_i.pos.split(":")[0]] = ["?"]
            base_b1 = ""
            base_b2 = ""
            if pos_b == -10 :
                pos_b = int(float(COR_STORE[COR_ORDER[i - 1]].pos.split(":")[0]))


    for i in range(0,len(COR_ORDER)) :
        SNP_key = COR_ORDER[i].split(":")[0]

        SNP_list = SNP_OUT[SNP_key]
        SNP_o_line = SNP_key
        for j in range(0,len(SNP_list)) :
            SNP_o_line += ('\t' + SNP_list[j])

        print >> OUTPUT_FH, SNP_o_line


if __name__ == '__main__' :
    main()













