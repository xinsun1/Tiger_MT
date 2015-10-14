__author__ = 'sunxin'


'''
This script calculate correlation between each SNP site.
'''

from take_snp import SNP

# SNP_FILE = ""


#==== read in SNP file ====

def read_snp(SNP_FILE) :

    global SNP_STORE
    SNP_STORE = {}
    SNP_FH = open(SNP_FILE)

    global SNP_ORDER
    SNP_ORDER = []

    while 1:
        SNP_L = SNP_FH.readline()

        if len(SNP_L) == 1 :  ##need to be debug
            break

        SNP_STORE[SNP(SNP_L).s] = SNP(SNP_L)    ## the pos of SNP should be
                                                  ## unique

        SNP_ORDER.append(SNP(SNP_L).s)

    if SNP_FH :
        SNP_FH.close()




#==== read in mapping result ====

import pysam as ps

# MAP_FILE = ""  ## take bam file and should be in
# dexed

def read_mp (MAP_FILE) :
    global MP_FH
    MP_FH = ps.AlignmentFile(MAP_FILE,"rb")




#==== calculate correlation between two SNPs or Indel ====

def cor_snp(pos_A,pos_B) :
    re_snp = {}
    ## initialize for two SNP
    for i in (SNP_STORE[pos_A].gt + SNP_STORE[pos_A].pt) :
        for j in (SNP_STORE[pos_B].gt + SNP_STORE[pos_B].pt) :
            re_snp[i + "/" + j] = 0


    pu = MP_FH.pileup(MP_FH.references, pos_A - 1, pos_B) ## the Ref is 0 based

    ## get reads mapped to posA,B
    for i in pu :
        if i.reference_pos == pos_A - 1 :
            L_A = i.pileups
        if i.reference_pos == pos_B - 1 :
            L_B = i.pileups


    ## take the shared reads
    L_Share = []

    for i in L_A :
        for j in L_B :
            if i.alignment == j.alignment :
                ## record base information
                # if SNP_STORE[pos_A].type == "snp" && SNP_STORE[pos_B].type == "snp" :
                #     base_A = i.alignment.query_sequence[i.query_position]
                #     base_B = j.alignment.query_sequence[j.query_position]
                #     re_snp[base_A + "/" + base_B] += 1
                #     break
                #
                if SNP_STORE[pos_A].type == "indel" :
                    if i.indel != 0 :
                        base_A = SNP_STORE[pos_A].pt[0]
                    else :
                        base_A = SNP_STORE[pos_A].gt
                else :
                    base_A = i.alignment.query_sequence[i.query_position]

                if SNP_STORE[pos_B].type == "indel" :
                    if i.indel != 0 :
                        base_B = SNP_STORE[pos_B].pt[0]
                    else :
                        base_B = SNP_STORE[pos_B].gt
                else :
                    base_B = j.alignment.query_sequence[j.query_position]

                re_snp[base_A + "/" + base_B] += 1
                break


    if len(L_Share) <= 10 :
        ## consider pair end information
        for i in L_A :
            i_mate = MP_FH.mate(i.alignment)
            for j in L_B :
                if i_mate == j.alignment :
                    ## record base information
                    # base_A = i.alignment.query_sequence[i.query_position]
                    # base_B = j.alignment.query_sequence[j.query_position]
                    # re_snp[base_A + "/" + base_B] += 1
                    # break

                    if SNP_STORE[pos_A].type == "indel" :
                        if i.indel != 0 :
                            base_A = SNP_STORE[pos_A].pt[0]
                        else :
                            base_A = SNP_STORE[pos_A].gt
                    else :
                        base_A = i.alignment.query_sequence[i.query_position]

                    if SNP_STORE[pos_B].type == "indel" :
                        if i.indel != 0 :
                            base_B = SNP_STORE[pos_B].pt[0]
                        else :
                            base_B = SNP_STORE[pos_B].gt
                    else :
                        base_B = j.alignment.query_sequence[j.query_position]

                    re_snp[base_A + "/" + base_B] += 1
                    break

    if len(L_Share) < 1 :
        print "Warning : No read across :", pos_A, "and ",pos_B
        return [re_snp,0]
    else :
        ## calculate correlation between posA & B
        re_sum = 0
        for i in re_snp.keys() :
            re_sum += re_snp[i]

        for i in re_snp.keys() :
            re_snp[i] = re_snp[i] / re_sum

        return [re_snp,re_sum]

# ######
# TACACGTACACACGTACACACGTACACACGTACACACGTA
# TACACGTACACACGTACACACGTACACACGTACACACGTACACACGTA
#
# f_pu = f.pileup('gi|313768475|ref|NC_014770.1|',200,400)
# for i in f_pu :
#     if i.reference_pos == 267 :
#         a = i.pileups
#
# for i in a :
#     if i.indel != 0 :
#         print i.indel
# ######


