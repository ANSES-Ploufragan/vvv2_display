#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# This file is part of the vvv2_display distribution (https://github.com/ANSES-Ploufragan/vvv2_display).
# Copyright (c) 2025 Fabrice Touzain.
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# FT - add vcfout output, a vcf file conatining only significant variant when requested by user: September 5th 2025
# FT - correct regexp part for description parsing: simplified, now do not miss variant: March 27th 2024
# FT - correct gene identif/name (position provided instead sometimes): June 15th 2023
# FT - last modification: September 19th 2022 to replace pyvcf by pysam
# AF - last modification: August 21st 2018
# Script name : convert_vcffile_to_readablefile2.py
#
# Description:  using a json file and a vcf file, return a tab delimited file
#               in order to feed a R script for SNP visualization
#               and another tab-delimited file which resume the reference and
#               alternative nucleotides at one position.
#
################################################################################

## ~ start script ~ ##

#############
## Modules ##
#############

# import vcf ## replaced by
from pysam import VariantFile
import os
from os import path
import numpy
import json
import argparse
import sys
import shutil
import re

# to be able to report line number in error messages
import inspect
frame = inspect.currentframe()

###############
## Debug var ##
###############
b_verbose = False
prog_tag = '[' + os.path.basename(__file__) + ']'


###############
## Functions ##
###############


def find_key_genes(dico_json, genomepos):
    """ This function retrieves the gene name and the start and end position of
this gene.
    """
    gene = "intergene"
    base_inf_allgenes = ''
    base_sup_allgenes = ''

    # print("find_key_genes call pos ("+str(genomepos)+")")
    for typ in dico_json:
        # print(f"dico_json key treated:{key}")

        if typ == "genes":
            for key in dico_json[typ]:
                base_inf , base_sup = dico_json[typ][key][0] , dico_json[typ][key][1]

                # if gene in reverse orientation, change coordinates to make base_inf < base_sup
                if base_sup < base_inf:
                    base_tmp = base_sup
                    base_sup = base_inf
                    base_inf = base_tmp

                if (base_inf <= genomepos)and(genomepos <= base_sup):
                    # print("GENE base_inf:"+str(base_inf)+" <= genomepos:"+str(genomepos)+" <= base_sup:"+str(base_sup))
                    if gene == 'intergene':
                        gene = key
                        base_inf_allgenes = str(base_inf)
                        base_sup_allgenes = str(base_sup)
                    else:
                        gene = gene+','+key
                        base_inf_allgenes = base_inf_allgenes+','+str(base_inf)
                        base_sup_allgenes = base_sup_allgenes+','+str(base_sup)
                #     break
                # else:
                #     print(f"NOT( base_inf:{base_inf} <= genomepos:{genomepos} <= base_sup:{base_sup} )")
                #     continue
    # print(f"return gene:{gene}\tbase_inf:{base_inf_allgenes}\tbase_sup:{base_sup_allgenes}")

    return gene, base_inf_allgenes, base_sup_allgenes

def find_key_proteins(dico_json, genomepos, gene_res):
    """ This function retrieves the protein name and the start and end position of
this protein.
    """
    protein = 'untranslated RNA'
    base_inf_allproteins = ''
    base_sup_allproteins = ''
    # print("find_key_proteins call pos ("+str(genomepos)+", "+gene_res+")")    
    for typ in dico_json:
        # print(f"dico_json key treated:{key}")

        if typ == "proteins":
            for key in dico_json[typ]:
                base_inf , base_sup = dico_json[typ][key][0] , dico_json[typ][key][1]

                 # if gene in reverse orientation, change coordinates to make base_inf < base_sup
                if base_sup < base_inf:
                    base_tmp = base_sup
                    base_sup = base_inf
                    base_inf = base_tmp

                if (base_inf <= genomepos)and(genomepos <= base_sup):
                    # print("PROT base_inf:"+str(base_inf)+" <= genomepos:"+str(genomepos)+" <= base_sup:"+str(base_sup))
                    if gene_res == 'intergene':
                        protein = 'intergene'
                        base_inf_allproteins = str(base_inf)
                        base_sup_allproteins = str(base_sup)
                    elif protein == 'untranslated RNA':
                        protein = key
                        base_inf_allproteins = str(base_inf)
                        base_sup_allproteins = str(base_sup)
                    else:
                        protein = protein+','+key
                        base_inf_allproteins = base_inf_allproteins+','+str(base_inf)
                        base_sup_allproteins = base_sup_allproteins+','+str(base_sup)
                #     break
                # else:
                #     print(f"NOT( base_inf:{base_inf} <= genomepos:{genomepos} <= base_sup:{base_sup} )")
                #     continue
    # print(f"return protein:{protein}\tbase_inf:{base_inf_allproteins}\tbase_sup:{base_sup_allproteins}")
    return protein, base_inf_allproteins, base_sup_allproteins


def write_line(temp_count, pos, gene, protein, dico):
    """ This function write the line that will be added in the final output file.
    'position\tSNP\tref\talt\tvariant_percent\tadd_ref\tadd_alt\tgene_id\tprotein_id\tsize_point\tlseq\trseq\thomopolymer\n'
    """
    if 'INTERGEN' in gene: # to replace the names INTERGEN1 -2 -3 etc...
        gene = 'intergen'
    global threshold
    global A

    if temp_count == 0:
        line = str(pos) + "\tNA\tNA\tNA\tNA\tNA\tNA\t"+gene+"\t"+protein+"\t1\tNA\tNA\tNA\tNA\n"
    elif temp_count == 1:
        ref, alt = dico[str(pos)][0], dico[str(pos)][1]
        freq, snp = dico[str(pos)][2], dico[str(pos)][3]
        lseq, rseq = dico[str(pos)][4], dico[str(pos)][5]

        if float(freq) >= threshold:
            A+=1
            is_homo = is_homopolymer(lseq, rseq, ref, alt)
            if is_homo is True:
                line = "\t".join( (str(pos), str(snp), ref, alt, freq, ref, alt, gene, protein, "1", str(A), lseq, rseq, "yes") ) + "\n"
            else:
                line = "\t".join( (str(pos), str(snp), ref, alt, freq, ref, alt, gene, protein, "1", str(A), lseq, rseq, "no") ) + "\n"
        else:
            line = "\t".join( (str(pos), str(snp), ref, alt, freq, "NA\tNA", gene, protein, "1\tNA\tNA\tNA\tNA\n") )
    return line

def is_homopolymer(lseq, rseq, ref, alt):
    """This function test if the region is homopolymeric.
    Return True if the region is considerd to be homopolymeric,
    otherwise, it will return false.
    The region is considered to homopolymeric if there is at least 3 identical
    nucleotides.
    """
    if ">" not in alt: # check if alt allele is not <DEL>
        if len(ref) < len(alt):# e.g. A vs AT
            if alt[-1] == rseq[0] and alt[-1] == rseq[1]: return True
        elif len(ref) > len(alt):# e.g. AT vs A
            if ref[-1] == rseq[0] and ref[-1] == rseq[1]: return True
        else:
            if ref[0] == rseq[0] and ref[0] == lseq[-1]: return True
            else: return False

    elif len(ref) == 1:
        if lseq[-1] == ref and rseq[0] == ref:
            return True
        elif lseq[-1] == ref and lseq[-2] == ref:
            return True
        elif rseq[0] == ref and rseq[1] == ref:
            return True
        else:
            return False

##########
## Main ##
##########

## call for the command line arguments ##
parser = argparse.ArgumentParser()
parser.add_argument('--vcfs', type = str, help='the vcf file to count')
parser.add_argument('--json', type = str, help = "json gene position file")
parser.add_argument('--out', type = str, help = "the output file")
parser.add_argument('--outs', type = str, help = "the output summary file")
parser.add_argument('--threshold', type = str, help='the threshold you use to cut the data', default = 0.1)
parser.add_argument('--vcfo', type = str, help='the vcf output file summarizing signfificant variants only')

args = parser.parse_args()

## check for the presence of all command line arguments ##
if not args.vcfs or not args.json or not args.out or not args.outs or not args.threshold:           
    parser.print_help()
    sys.exit(prog_tag + "\nAn error occured while entering the arguments.\nPlease, read the help section above.\n")
else:
    # recover genomesize data
    threshold = float(args.threshold)
    filin = open(args.json,"r")
    dico_json = json.load(filin)
    if b_verbose:
        print("dico_json:"+str(dico_json))
    
    genomesize = dico_json["genomesize"]
    filin.close()

    # lists initializing
    genomeposition = []
    temp_counts = []

    # dataframe
    snp_positions = numpy.zeros(genomesize+1, dtype=bool)
    
    # recover the snp position inside the vcf file
#    file_vcf = open(args.vcfs, "r")         # for pyvcf deprecated

    try:
        file_vcf = VariantFile(args.vcfs)
    except ValueError:
        print(prog_tag + " No variant found")
        shutil.copy(args.vcfs, args.out)
        print(prog_tag + ' ' + args.out +" file created")
        sys.exit()
        
#    for record in vcf.Reader(file_vcf):     # use pyvcf deprecated
#        snp_positions[record.POS-1] = True  # use pyvcf deprecated
    for record in file_vcf.fetch():
        snp_positions[record.pos-1] = True
        # for each position written inside the vcf file,
        # the value in the dataframe is changed from False to True
    file_vcf.close()

    # turn the False or True values from the dataframe into 0 or 1
    for n in range(genomesize):
        temp_counts.append(numpy.count_nonzero(snp_positions[n]))
        genomeposition.append(n+1)

    # print("genomeposition:"+str(genomeposition)+", line "+str(frame.f_lineno))

    # look for 3 groups: SNP_position , REF , ALT
    regex1 = r'[a-zA-Z0-9\._]+[\t]([0-9]+)[\t][a-zA-Z0-9\._-]+[\t]([ATCGKMSWRYBDHVN\.]+)\t(([ATCGKMSWRYBDHVN\,\.]+)||<DUP>||<DEL>||<INV>||<INS>||<FUS>)\t[0-9]+\t[A-Za-z0-9\.\;\,]+\tSAMPLE='
    # look for 1 group: Variant frequency
    regex2 = r'[0-9\/\,\.]+:[0-9\/\,\.]+:[0-9\/\,\.]+:[0-9\/\,\.]+:([0-9\/\,\.]+):[0-9\/\,\.]+:[0-9\.\/\,\.%]+'

    ## regex compilation
    reg1 , reg2 = re.compile(regex1) , re.compile(regex2)
    # get the left sequence and the right sequence from the reference base
    #           part added to left regexp:---------                       and right:--- and not finished by ;
    lseq, rseq = re.compile(r";LSEQ=([A-Z0\.\:0-9\-]+);"), re.compile(r";RSEQ=([A-Z0SNV]+)")
    # lseq, rseq = re.compile(r";LSEQ=([A-Z0]+);"), re.compile(r";RSEQ=([A-Z0]+);")

    # dictionnary initialization to store data
    dico = {} # forward in the script key = snp_location

    # parse vcf file and harvest snp data
    with open(args.vcfs, "r") as file_vcf:
        for line in file_vcf:
            if "#" in line: # skip the header section
                continue
            else: # data section
                new_list = [] # add in this order = ref, alt, freq, snp_number, lseq, rseq

                # print(line)

                # search for the motif and capture the wanted group
                SNP_position = (reg1.search(line)).group(1)
                ref = (reg1.search(line)).group(2) ; new_list.append(ref)
                alt = (reg1.search(line)).group(3) ; new_list.append(alt)
                freq = (reg2.search(line)).group(1) ; new_list.append(freq)
                snp = len(alt.split(",")) ; new_list.append(snp)
                try:
                    rightseq = (rseq.search(line)).group(1); 
                except AttributeError:
                    sys.exit("")
                # added 2024 12 05
                # special case of vcf vardict result usually occuring at le last positions of the virus, no seq are provided
                # for both leftseq (ACCNR:coordinate instead, regexp changed for this) and rightseq (SNV)
                if rightseq == 'SNV':
                    rightseq = ''
                else:
                    leftseq = (lseq.search(line)).group(1); 
                new_list.append(rightseq)
                new_list.append(leftseq)
                dico[SNP_position] = new_list

    # Once harvested, data need to be recorded in a file

    # TODO: get the vcf header to write it in the output file
    if b_create_summary_vcf:
        filoutvcf = open(vcfo, "w")

    A = 0 # this flag is used in the write_line function in order to add indices to line where variants are upper than threshold
    with open(args.out, "w") as filout:
        summary_list = []
        regex = r'([0-9]+)\t1\t([A-Z\>\<]+)\t([A-Z\>\<]+)\t([0-9\.]+)\t[A-Z\>\<]+\t[A-Z\>\<]+\t([A-Z0-9a-z\-_\, /]+)\t([^\t]*)\t1\t([0-9]+)\t([A-Z]+\t[A-Z]+\t(no|yes))\n'
        REGEX = re.compile(regex)    
        filout.write("position\tSNP\tref\talt\tvariant_percent\tadd_ref\tadd_alt\tgene_id\tprotein_id\tsize_point\tindice\tlseq\trseq\tisHomo\n")
        a = 0 # flag to find the first genomic region after initialization of i

        # create a vcf file of only significant variant for SnpEff and NextClade downlstream analyses
        if b_create_summary_vcf:
            filoutvcf.write(vcf_header)

        # this part is required to write the correct name of the genomic region
        for i in range(len(genomeposition)):
            pos = genomeposition[i]
            gene_id, base_inf, base_sup = find_key_genes(dico_json, pos)
            protein_id, pbase_inf, pbase_sup = find_key_proteins(dico_json, pos, gene_id)
            # print("pos gene_id protein_id:"+str(pos)+"\t"+gene_id+"\t"+protein_id)
            line = write_line(temp_counts[i], pos, gene_id, protein_id, dico)
            filout.write(line)

            # generate another file, which resume the variation
            #print(line)

            if REGEX.match(line):
                indice = (REGEX.search(line)).group(7)
                gene = (REGEX.search(line)).group(5)
                protein = (REGEX.search(line)).group(6)
                ref = (REGEX.search(line)).group(2)
                alt = (REGEX.search(line)).group(3)
                position = (REGEX.search(line)).group(1)
                freq = (REGEX.search(line)).group(4)
                homo = (REGEX.search(line)).group(8)
                new_line = str("\t".join( (indice, position, ref, alt, freq, gene, protein, homo) ))
                summary_list.append(new_line + "\n")

                # write in the vcf file of only significant variant for SnpEff and NextClade downlstream analyses
                filoutvcf.write(line)

    with open(args.outs, "w") as filout:
        filout.write("indice\tposition\tref\talt\tfreq\tgene\tprot\tlseq\trseq\tisHomo*\n")
        for line in summary_list:
            filout.write(line)
        filout.write("""
*NB: an homopolymer region is set to 'yes' if there is a succession of at least 3 identical nucleotides.
     it looks like a restrictive measure, but Ion Torrent and Nanopore sequencing are very bad on such region, so make sure you verify these variants.""")
    print(prog_tag + ' '+ args.out +" file created")
    print(prog_tag + ' '+ args.outs +" file created")
    if b_create_summary_vcf:
        print(prog_tag + ' '+ args.vcfo +" file created")
## ~ end of script ~ ##
