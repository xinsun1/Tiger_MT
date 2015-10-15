__author__ = 'sunxin'

#!/usr/bin/python

'''
This script take in the output from ana_cor.py.
Sort the SNP sites for two haplotype.
'''
import sys
from take_snp import COR

def read_ana (ANA_FILE) :
    ANA_FH = open(ANA_FILE,'r')

    ANA_STORE = {}
    ANA_ORDER = []
    ANA_F = {}
    while 1:
        ANA_L = ANA_FH.readline().strip()
        if len(ANA_L) == 0 :
            break

        ANA_l = ANA_L.split('\t')

        ANA_ORDER.append(ANA_l[0])
        ANA_F[ANA_l[0]] = 1 ## 1 : no flip; -1 : flip
        if ANA_l[1] == "?" :
            ANA_STORE[ANA_l[0]] = ANA_l[1]
        else :
            ANA_STORE[ANA_l[0]] = [ANA_l[1],ANA_l[2]]

    if ANA_FH :
        ANA_FH.close()

    return [ANA_STORE,ANA_ORDER,ANA_F]



def main () :
    print "take in MOD_FILE, ANA_FILE"

    if len(sys.argv) == 1 :
        exit(0)

    MOD_FILE = sys.argv[1]
    ANA_FILE = sys.argv[2]

    global ANA_STORE, ANA_F, ANA_ORDER

    [ANA_STORE,ANA_ORDER,ANA_F] = read_ana(ANA_FILE)

    MOD_FH = open(MOD_FILE, 'r')

    MOD_NEW_FH = open("mod2_" + MOD_FILE, 'w')
    ANA_NEW_FH = open("mod2_" + ANA_FILE, 'w')

    while 1 :
        MOD_L = MOD_FH.readline()

        if len(MOD_L) == 0 :
            break

        mod_cor = COR(MOD_L)
        if mod_cor.n_cor == 2 :
            base_11 = mod_cor.cor.keys()[0].split("/")[0]
            base_12 = mod_cor.cor.keys()[0].split("/")[1]
            base_21 = mod_cor.cor.keys()[1].split("/")[0]
            base_22 = mod_cor.cor.keys()[1].split("/")[1]

            if base_11 != base_21 and base_12 != base_22 :
                pos_1 = mod_cor.pos.split(":")[0]
                pos_2 = mod_cor.pos.split(":")[1]

                bp11 = ANA_STORE[pos_1][0]
                bp12 = ANA_STORE[pos_1][1]
                bp21 = ANA_STORE[pos_2][0]
                bp22 = ANA_STORE[pos_2][1]

                if bp11 == base_11 :
                    if bp21 == base_12 :
                        continue
                    else :
                        for i in ANA_ORDER[ANA_ORDER.index(pos_2) : ] :
                            ANA_F[i] *= -1

                else : ## bp11 = base_21
                    if bp21 == base_22 :
                        continue
                    else :
                        for i in ANA_ORDER[ANA_ORDER.index(pos_2) : ] :
                            ANA_F[i] *= -1

            else :
                print >> MOD_NEW_FH, MOD_L

        else :
            print >> MOD_NEW_FH, MOD_L

    #==== output new ANA ====
    for i in range(0, len(ANA_ORDER)) :
        if ANA_F[i] == 1 :
            if ANA_STORE[i][0] == "?" :
                print >> ANA_NEW_FH, i + '\t' + ANA_STORE[i][0]
            else :
                print >> ANA_NEW_FH, i + '\t' + ANA_STORE[i][0] + '\t' + ANA_STORE[i][1]
        else :
            if ANA_STORE[i][0] == "?" :
                print >> ANA_NEW_FH, i + '\t' + ANA_STORE[i][0]
            else :
                print >> ANA_NEW_FH, i + '\t' + ANA_STORE[i][1] + '\t' + ANA_STORE[i][0]


    if MOD_FH :
        MOD_FH.close()
    if ANA_NEW_FH :
        ANA_NEW_FH.close()
    if MOD_NEW_FH :
        MOD_NEW_FH.close()



if __name__ == '__main__' :
    main()
