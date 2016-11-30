#!/usr/bin/python3
#import sys
#infile = sys.stdin
#infile = file('/mnt/userdata/kinzie/RP_SSC_Parents_Mosaic/ssc_sample_5k.txt', 'r')

#outfile = file('/mnt/userdata/kinzie/RP_SSC_Parents_Mosaic/ssc_sample_5k_mvparentsclean.txt', 'w')
import os
os.chdir('C:\\Users\\plab\\Documents\\Kyra')
infile = file('ssc_sample_5k.txt')
from scipy import stats
#for the 24tr pc: module add anaconda


#Dictionaries
col_map_dict = {}
sample_map_dict = {}
parent_dict = {}
child_dict = {}
family_dict = {}

all_unique_dict = {}
unique_inherited_dict = {}
unique_not_inherited_dict = {}
#unique_chr_dict = {}

bi_test_dict = {}
binom_result = []

unique_id = ''
chr_position = ''

for line in infile:
    line = line.rstrip()
    if line.startswith('##'):
        #print(line) #Header
        continue
    else:
        line = line.split('\t')
        if line[0][0] == '#':
            for i in range(9, len(line)):
                sample = line[i]                                        #ie - 11000.fa
                col_map_dict[i] = sample                                #ie - 9, 11000.fa (columns:sample index)
                sample_map_dict[sample] = i                             #ie - 11000.fa, 9 (sample index: columns)
                if sample.endswith('fa') or sample.endswith('mo'):
                    parent_dict[i] = sample                             #ie - (column: parent index)
                else:
                    child_dict[i] = sample                              #ie - (column: child index)
                for parent in parent_dict.values():
                    for keys, child in child_dict.items():
                        if parent[:-3] == child[:-3]:
                            if parent not in family_dict:
                                family_dict[parent] = set()
                            family_dict[parent].add(keys)
                            #print family_dict             #ie - (parent index: children's columns in set)
            print('\t'.join(line))
            continue
        else:
            parent_with_unique_allele = False
            #unique_id = ''                                      #Hashed out purposefully.see above under dictionaries
            alt_alleles = line[4].split(',')
            chr_position = line[1]

            for num_alt_alleles in range(len(alt_alleles)):
                unique_allele = False
                alt_allele_count = 0
                for i in range(9, len(line)):
                    if i not in parent_dict:
                        continue
                    genotype = line[i]
                    alt_alleles_ids = num_alt_alleles + 1
                    genotype_info = line[i].split(':')
                    genotype_number = genotype_info[0]
                    if genotype_number != '0/0/0/0/0' and genotype_number != '././././.':
                        present = str(alt_alleles_ids) in genotype_number
                        if present:
                            alt_allele_count += 1
                            unique_id = col_map_dict[i]
                            unique_allele = alt_alleles_ids
                if alt_allele_count == 1:
                    parent_with_unique_allele = True
                if parent_with_unique_allele:
                    for i in range(9, len(line)):
                        if i in parent_dict:
                            continue
                        genotype = line[i]
                        alt_alleles_ids = num_alt_alleles + 1
                        genotype_info = line[i].split(':')
                        genotype_number = genotype_info[0]
                        present = str(alt_alleles_ids) in genotype_number
                        if present:
                            try:
                                if i not in family_dict[unique_id]:
                                    continue
                            except KeyError:
                                pass
                            try:
                                if i in family_dict[unique_id]:
                                    #if i not in unique_inherited_dict:
                                    if chr_position not in unique_inherited_dict:
                                        unique_inherited_dict[sample_map_dict[unique_id]] = set()
                                        unique_inherited_dict[sample_map_dict[unique_id]].add(chr_position)
                                        #unique_inherited_dict[sample_map_dict[unique_id]] = chr_position
                                        #print unique_inherited_dict.keys()
                                        #print 'Look Here'
                                    #if i not in all_unique_dict:
                                    if chr_position not in all_unique_dict:
                                        all_unique_dict[sample_map_dict[unique_id]] = set()
                                        all_unique_dict[sample_map_dict[unique_id]].add(chr_position)
                                        #print all_unique_dict
                            except KeyError:
                                pass

                        if present == False:
                            if chr_position not in unique_not_inherited_dict:
                                unique_not_inherited_dict[sample_map_dict[unique_id]] = set()
                                unique_not_inherited_dict[sample_map_dict[unique_id]].add(chr_position)
                                print unique_not_inherited_dict
                                print 'Hello2'
                            if chr_position in all_unique_dict:
                                all_unique_dict[sample_map_dict[unique_id]].add(chr_position)
                            if chr_position not in all_unique_dict:
                                all_unique_dict[sample_map_dict[unique_id]] = set()
                                all_unique_dict[sample_map_dict[unique_id]].add(chr_position)
                                print all_unique_dict
                                #print 'Look HERE'
                                #print unique_inherited_dict
                            #if i not in all_unique_dict:
                                #all_unique_dict[i] = set()
                                #all_unique_dict[i].add(chr_position)-
                    #print all_unique_dict


           if parent_with_unique_allele:
                for i in unique_not_inherited_dict and for i in unique_inherited_dict:
                    chr_position = line[1]
                    genotype = line [i]
                    genotype_info = line[i].split(':')
                    reads = genotype_info[1].split(',')
                    ref_reads = int(reads[0])
                    all_alt_reads = map(int, reads[1:])
                    alt_read = int(reads[unique_allele])
                    total_reads = ref_reads + all_alt_reads
                    binom_result = stats.binom_test(alt_read, total_reads) #should be alt_reads versus total_reads. before it was ref_reads versus total_reads.
                    #print(unique_id)
                    #print(binom_result)
                    if binom_result >= 0.05:
                        continue
                    else:
                        bi_test_dict[unique_id] = binom_result
                        print('\t'.join([chr_position, unique_id, str(binom_result)]))
                        #print('\t'.join([chr_position, str(bi_test_dict)]))
                        continue
#For line 108: you want to apply the binomial test only to unique alleles within the parents, so you want to store the
#unique alleles in another variable and then apply them in the binomial test so you're only looking at the integer of the
# reads of the unique alleles for the parents.      '''
