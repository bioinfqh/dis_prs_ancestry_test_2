from snps import SNPs
import sys

s = SNPs(sys.argv[1])
saved_snps = s.save(sys.argv[2], vcf=True)
#saved_snps = s.save("/media/quirin/INTENSO/new_ancestry/tmp/converted.vcf", vcf=True)
