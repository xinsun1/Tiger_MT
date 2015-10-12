__author__ = 'sunxin'


#!/usr/bin/python


'''
This script take in the reference file, sam file (already mapped to the reference),
 the SNP index file (from samtools).
Read the index and output each SNP site information in each read.
'''


__import__(#package used to read sam file#)
__import__(#package used to read fastq file#)

#====read in files====
def readinfile(Ref,Sam,IndexFile) :
    Ref_FH = open(Ref)
    Index_FH = open(IndexFile)
    Sam_FH = open(Sam)




#====Class for SNP info====

class SNP:
    def __init__(self, s=0, type="snp", gt='', pt='', n_pt=2) :
        self.s = self.strip().split('\t')[0]
        self.gt = self.strip().split('\t')[1]
        if len(self.gt) > 1: self.type = "indel"
        self.pt = self.strip().split('\t')[2]
        self.n_pt = len(self.pt.split(","))




#====Class for PE read info====
class PE:
    def __init__(self) :
        self.readline






if __name__ == '__main__':
    ###




