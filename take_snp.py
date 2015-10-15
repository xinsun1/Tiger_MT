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
    def __init__(self, s=0, type=None, gt=None, pt=None, n_pt=2) :
        super(SNP, self).__init__()
        self.s = self.strip().split('\t')[0]
        self.gt = self.strip().split('\t')[1]
        if len(self.gt) > 1:
            self.type = "indel"
        else :
            self.type = "snp"
        self.pt = self.strip().split('\t')[2].split(",")
        self.n_pt = len(self.pt)



#==== Class for cor info ====
class COR (str) :
    def __init__(self, pos=None, n_cor = None, n_read = None, cor = None) :
        super(COR, self).__init__()
        a = self.strip().split('\t')
        self.pos = a[0]
        self.n_read = int(float(a[-1]))
        self.n_cor = len(a[1:-1]) / 2
        self.cor = {}
        for i in range(0, self.n_cor -1) :
            self.cor[a[2*i + 1]] = a[2*i + 2]














if __name__ == '__main__':
    main()




