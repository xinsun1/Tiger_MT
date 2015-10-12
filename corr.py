__author__ = 'sunxin'


'''
This script calculate correlation between each SNP site.
'''

from take_snp import SNP

SNP_FILE = ""


#====read in SNP file====


SNP_STORE = {}
def read_snp(SNP_FILE) :
    SNP_FH = open(SNP_FILE)

    while 1:
        SNP_L = SNP_FH.readline()

        if len(SNP_L) == 0:  ##need to be debug
            break

        SNP_STORE[SNP(SNP_FH).s] = SNP(SNP_FH)    ## the pos of SNP should be
                                                  ## unique
    if SNP_FH :
        SNP_FH.close()

#====read in mapping result====






