# -*- coding: UTF-8 -*-
__author__ = 'konradjk'

'''将vcf文件中的vep注释提取并整理出来。没有的注释全部认为是NA'''

import argparse
import gzip
import re
import sys

print >> sys.stderr, "mixed indentation must use python 2.2 below"

def main(args):
    f = gzip.open(args.vcf) if args.vcf.endswith('.gz') else open(args.vcf)
    vep_field_names = None
    header = None

    # xin print the header I need
    # print >> sys.stdout, "CHROM","POS","REF","ALT","Allele","Consequence","Gene","Feature"
    vep_f = args.vep_field

    for line in f:
        line = line.strip()

        # Reading header lines to get VEP and individual arrays
        if line.startswith('#'):
            line = line.lstrip('#')
            if line.find('ID=CSQ') > -1:  #find如果找到则返回找到的位置，如果没找到返回-1
                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
		if vep_f is None:
			vep_f = vep_field_names
		else:
			vep_f = [vep_f]
		print >> sys.stdout, "CHROM\t","POS\t","REF\t","ALT\t","\t".join(vep_f)
            if line.startswith('CHROM'):
                header = line.split()
                header = dict(zip(header, range(len(header))))
            continue

        if vep_field_names is None:
            print >> sys.stderr, "VCF file does not have a VEP header line. Exiting."
            sys.exit(1)
        if header is None:
            print >> sys.stderr, "VCF file does not have a header line (CHROM POS etc.). Exiting."
            sys.exit(1)

        # Pull out annotation info from INFO and ALT fields
        fields = line.split('\t')
        info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[header['INFO']])]) #re.split是对有多个分割符的字符串进行分割,支持正则表达式,是re包内的

        # Only reading lines with an annotation after this point
        if 'CSQ' not in info_field: continue
        annotations = [dict(zip(vep_field_names, x.split('|'))) for x in info_field['CSQ'].split(',') if len(vep_field_names) == len(x.split('|'))]
        # lof_annotations = [x for x in annotations if x['LoF'] == 'HC']

        # Code to process annotations and VCF line goes here...

        # annotations = list of dicts, each corresponding to a transcript-allele pair
        # (each dict in annotations contains keys from vep_field_names)
        # lof_annotations = list of dicts, only high-confidence LoF
        # if len(lof_annotations) > 0 can determine if at least one allele is LoF for at least one transcript
        # fields = vcf line, can be accessed using:
        # fields[header['CHROM']] for chromosome,
        # fields[header['ALT']] for alt allele,
        # or samples using sample names, as fields[header['sample1_name']]
	
	for x in annotations:
		#consequence_field = x['Consequence'].split('&')
		#for y in consequence_field:
		#no need to split, will process later
			# print >> sys.stdout, fields[header['CHROM']],fields[header['POS']],fields[header['REF']],fields[header['ALT']],x['Allele'],x['Consequence'],x['Gene'],x['Feature']
			for y in vep_f:
				if x[y] is None or x[y] is "":
					x[y] = 'NA'
			print >> sys.stdout, fields[header['CHROM']],"\t",fields[header['POS']],"\t",fields[header['REF']],"\t",fields[header['ALT']],"\t","\t".join([x[y] for y in vep_f])

    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file (from VEP+LoF); may be gzipped', required=True)
    parser.add_argument('--vep_field', '-v', help='Field of vep to output', required=False)
    args = parser.parse_args()
    main(args)
