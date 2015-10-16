__author__ = 'sunxin'


import pysam as ps
import sys

from corr import read_snp

'''
This script generate the haplotype MT sequence using
 the analysed correlation for each SNP.
'''


def main () :

    print "take in Ref_FILE, MOD2_FILE, SNP_FILE, OUT_FILE"

    if len(sys.argv) == 1 :
        exit(0)

    REF_FILE = sys.argv[1]
    MOD2_FILE = sys.argv[2]
    SNP_FILE = sys.argv[3]
    OUT_FILE = sys.argv[4]

    REF_FH = ps.FastaFile(REF_FILE)
    ref = REF_FH.fetch(REF_FH.references[0],0,REF_FH.get_reference_length(REF_FH.references[0]))
    con_1 = ref
    con_2 = ref

    global SNP_STORE, SNP_ORDER

    [SNP_STORE, SNP_ORDER] = read_snp(SNP_FILE)

    MOD2_FH = open(MOD2_FILE, 'r')

    while 1 :
        MOD_L = MOD2_FH.readline()

        if len(MOD_L) == 0 :
            break

        mod_l = MOD_L.strip().split('\t')

        if len(mod_l[1]) != 1 or len(mod_l[2]) != 1 :
            len_gt = len(SNP_STORE[mod_l[0]].gt)
            con_1 = con_1[ : (mod_l[0] - 1)] + mod_l[1] + con_1[(mod_l[0] + len_gt - 1) : ]
            con_2 = con_2[ : (mod_l[0] - 1)] + mod_l[2] + con_2[(mod_l[0] + len_gt - 1) : ]
        else :
            con_1 = con_1[ : (mod_l[0] - 1)] + mod_l[1] + con_1[mod_l[0] : ]
            con_2 = con_2[ : (mod_l[0] - 1)] + mod_l[2] + con_2[mod_l[0] : ]

    OUT_FH = open(OUT_FILE, 'w')

    print >> OUT_FILE, ">consensus_1"
    print >> OUT_FILE, con_1
    print >> OUT_FILE, ">consensus_2"
    print >> OUT_FILE, con_2

    if MOD2_FH :
        MOD2_FH.close()
    if OUT_FH :
        OUT_FH.close()



if __name__ == '__main__' :
    main()