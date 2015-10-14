__author__ = 'sunxin'


#!/usr/bin/python


'''
This script take in the reference file, sam file (already mapped to the reference),
 the SNP index file (from samtools).
Read the index and output each SNP site information in each read.
'''

def main () :
    exit(0)



#==== Class for SNP info ====

class SNP (str):
    def __init__(self, s=0, type="snp", gt=None, pt=None, n_pt=2) :
        super(SNP,self).__init__()
        self.s = self.strip().split('\t')[0]
        self.gt = self.strip().split('\t')[1]
        if len(self.gt) > 1: self.type = "indel"
        self.pt = self.strip().split('\t')[2].split(",")
        self.n_pt = len(self.pt)





if __name__ == '__main__':
    main()




