__author__ = 'sunxin'

#!/usr/bin/python
from corr import *
import sys

#==== main ====

def main() :
    print "take in by order : SNP_file(.txt), BAM_file(indexed), output name"

    if len(sys.argv) == 1 :
        exit(0)

    SNP_FILE = sys.argv[1]
    BAM_FILE = sys.argv[2]
    OUTPUT_FILE = sys.argv[3]

    read_snp(SNP_FILE)

    read_mp(BAM_FILE)

    global SNP_ORDER

    OUTPUT_FH = open(OUTPUT_FILE,'w')
    for i in range(0,len(SNP_ORDER) - 1 ) :
        corr_o_hash = cor_snp(SNP_ORDER[i],SNP_ORDER[i + 1])[0]

        corr_o_sum = cor_snp(SNP_ORDER[i],SNP_ORDER[i + 1])[1]

        corr_o_line = ""

        for i in corr_o_hash.keys() :
            if corr_o_hash[i] != 0 :
                corr_o_line = corr_o_line + i + '\t' + str(corr_o_hash[i]) + '\t'

        corr_o_line += corr_o_sum

        print >> OUTPUT_FH, corr_o_line

    if OUTPUT_FH :
        OUTPUT_FH.close()




if __name__ == '__main__':
    main()



