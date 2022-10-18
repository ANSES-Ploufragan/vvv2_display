#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###
# USE PYTHON3
# from tbl file of vadr annotator (2 inputs), creates a json output file of
# annotations (output)
###
import argparse, os, sys, csv, re, warnings
from os import path
import subprocess

# to be able to report line number in error messages
import inspect
frame = inspect.currentframe()

# debug
b_test_convert_tbl2json = False # ok 2022 04 14
b_test = False

prog_tag = '[' + os.path.basename(__file__) + ']'


pass_annot_f = '' # in, pass tbl file from vadr
fail_annot_f = '' # in, fail tbl file from vadr
seq_stat_f   = '' # in, seq stat file from vadr
json_annot_f = '' # out
# bed_annot_f  = '' # out
bed_vardict_annot_f  = '' # out same as previous but positions on contigs not shifted according to previous contigs

material_type_list = ['gene',
                      'misc_feature',
                      'mat_peptide']
material_type = '|'.join(material_type_list)

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pass_annot_f", dest='pass_annot_f',
                    help="tabular file of vadr annotations, with pass status",
                    metavar="FILE")
parser.add_argument("-f", "--fail_annot_f", dest='fail_annot_f',
                    help="tabular file of vadr annotations, with fail status",
                    metavar="FILE")
parser.add_argument("-s", "--seq_stat_f", dest='seq_stat_f',
                    help="seq stat file of vadr annotator",
                    metavar="FILE")
parser.add_argument("-j", "--json_out_f", dest='json_annot_f',
                    help="json output file of annotations",
                    metavar="FILE")
# parser.add_argument("-b", "--bed_out_f", dest='bed_annot_f',
#                     help="bed output file of annotations",
#                     metavar="FILE")
parser.add_argument("-c", "--bed_vardict_out_f", dest='bed_vardict_annot_f',
                    help="bed output file of annotations for vardict",
                    metavar="FILE")
parser.add_argument("-z", "--test_all", dest='b_test_all',
                    help="[Optional] run all tests. Additionally, with --load_ncbi_tax_f, allow to download ncbi ete3 tax db the first time you use the script",
                    action='store_true')
parser.add_argument("-v", "--verbose", dest='b_verbose',
                    help="[Optional] To have details on records when running",
                    action='store_true')
parser.set_defaults(b_test_all=False)
parser.set_defaults(b_verbose=False)

# get absolute path in case of files
args = parser.parse_args()

# -------------------------------------------
# check arguments
b_test_all = args.b_test_all

if b_test_all:
    b_test_convert_tbl2json = True
    b_test = True
else:
    b_test = b_test_convert_tbl2json

if ((not b_test)and
    ((len(sys.argv) < 9) or (len(sys.argv) > 11))):
    parser.print_help()
    print(prog_tag + "[Error] we found "+str(len(sys.argv)) +
          " arguments, exit line "+str(frame.f_lineno))
    sys.exit(0)

# print('args:', args)
if args.pass_annot_f is not None:
    pass_annot_f = os.path.abspath(args.pass_annot_f)
elif(not b_test):
    sys.exit("[Error] You must provide pass_annot_f")
if args.fail_annot_f is not None:
    fail_annot_f = os.path.abspath(args.fail_annot_f)
elif(not b_test):
    sys.exit("[Error] You must provide fail_annot_f")
if args.seq_stat_f is not None:
    seq_stat_f = os.path.abspath(args.seq_stat_f)
elif(not b_test):
    sys.exit("[Error] You must provide seq_stat_f")
if args.json_annot_f is not None:
    json_annot_f = os.path.abspath(args.json_annot_f)
elif(not b_test):
    sys.exit("[Error] You must provide json_annot_f")
# if args.bed_annot_f is not None:
#     bed_annot_f = os.path.abspath(args.bed_annot_f)
# elif(not b_test):
#     sys.exit("[Error] You must provide bed_annot_f")
if args.bed_vardict_annot_f is not None:
    bed_vardict_annot_f = os.path.abspath(args.bed_vardict_annot_f)
elif(not b_test):
    sys.exit("[Error] You must provide bed_vardict_annot_f")

if args.b_verbose is not None:
    b_verbose = args.b_verbose


if b_test_convert_tbl2json:
    test_dir = 'test_convert_tbl2json/'
    test_name = [
        'res'
        # 'res2'
    ]
    for resn in test_name:
        pass_annot_f = f"{test_dir}{resn}.vadr.pass.tbl"
        fail_annot_f = f"{test_dir}{resn}.vadr.fail.tbl"
        json_annot_f = f"{test_dir}{resn}.vadr.json"
#        bed_annot_f  = f"{test_dir}{resn}.vadr.bed"
        bed_vardict_annot_f  = f"{test_dir}{resn}.vadr.4vardict.bed"                
        seq_stat_f   = f"{test_dir}{resn}.vadr.seqstat"
        cmd = ' '.join(['./convert_tbl2json.py',
                        f"-p {pass_annot_f}",
                        f"-f {fail_annot_f}",
                        f"-s {seq_stat_f}",                    
                        f"-j {json_annot_f}",
#                       f"-b {bed_annot_f}",
                        f"-c {bed_vardict_annot_f}"                        
                        ])
        if b_verbose:
            cmd = f"{cmd} -v"
        print(f"cmd:{cmd}")
        print("START")
        os.system(cmd)
        print("END")
    sys.exit()

b_next_is_gene = False
b_next_is_product = False
b_next_is_protein_id = False
b_next_is_note = False
b_additional = False

# final var
# closest_ref = ''

# Returns all indexes of an item in a list or a string"
def indexlist(item2find, list_or_string):
  return [n for n,item in enumerate(list_or_string) if item==item2find]

# ----------------------------------------------------------------------
# reads seqstat file to get info on genome_length and on contigs lengths
# ----------------------------------------------------------------------
contig_names   = []
contig_lengths = []

with open(seq_stat_f) as ssf:
    for line in ssf:
        if re.match(r"Total # residues:    ", line):
            line_fields =  line.split(' ')
            # print("line_fields:"+','.join(line_fields))
            genome_length = line.split(' ')[6]
            genome_length = genome_length.rstrip()
            if b_verbose:
                print(f"genome_length: {genome_length}")
        elif re.match(r"= ", line):
            # get contigs size to deduce value of contig_pos_shift for each: example of file:
            # = 1                            13356 
            # = 4                              541 
            # = 5                            13761
            line_fields =  line.split()
            contig_name = line_fields[1]
            contig_length = line_fields[2]
            contig_length.rstrip()
            # print("line_fields:"+','.join(line_fields))
            contig_names.append(contig_name)
            contig_lengths.append(contig_length)
            if b_verbose:
                print(f"contig {contig_name} length: {contig_length}")

ssf.close()
# ----------------------------------------------------------------------

    
# tmp var
contig = ''
gene_name = ''
gene_start = ''
gene_end = ''
cds_start = ''
cds_end = ''
product = ''
protein_id = ''
misc_feature_start = ''
misc_feature_end = ''
note = ''
curr_type = ''  # 'CDS' or 'mat_peptide'
ori_end = -1    # store end of CDS before correction due to detected STOP codon

# var to count of nt to shift position for contigs those are not first contigs in the assembly
contig_pos_shift = 0 

# stored info
names = []
types = []
starts= []
ends  = []
starts_vardict = []
ends_vardict   = []
chrs  = []

# pattern to remove non alphanumeric character (ie: <> character befire start end positions of genes, cds, etc...)
non_alphanum = re.compile("\W")

# for each annotation file created by vadr
for annot_f in [pass_annot_f, fail_annot_f]:

    # read annotation file
    print(f"read {annot_f} file")
    with open(annot_f) as paf:
        for line in paf:
            try:
                line_fields = line.split()

                if b_verbose:
                    print(' '.join(["line_fields:",
                                    ','.join(line_fields)
                                    ]))
                    print(f"line_fields 2: {line_fields[2]}")

                # -------------------------------------------------------------
                # treat 'next' lines when something is anounced the line before
                if b_next_is_gene:
                    gene_name = line_fields[1]
                    print(' '.join(['gene',
                                   gene_start,
                                   gene_end,
                                   gene_name]))
                    b_next_is_gene = False
                    # store info
                    chrs.append(contig)
                    names.append(gene_name)
                    types.append('gene')
                    starts.append(int(gene_start) + contig_pos_shift)
                    ends.append(int(gene_end) + contig_pos_shift)
                    starts_vardict.append(int(gene_start))
                    ends_vardict.append(int(gene_end))

                elif b_next_is_product:
                    if line_fields[0] == "product":
                        product = line_fields[1]
                        b_next_is_product = False
                        b_next_is_protein_id = True
                    elif len(line_fields) == 2:
                        # means only coordinates because gene longer than CDS
                        print("CDS shorter than gene, skip useless coordinates")
                    # -----------------------------------------------------------------------------------
                    # this part until else MUTS be useless, but bug found only when b_verbose is True...
                    elif (len(line_fields) == 3) and line_fields[2] == 'misc_feature':
                        misc_feature_start = line_fields[0]
                        misc_feature_end = line_fields[1]
                        b_next_is_note = True
                        b_next_is_product = False
                    elif (len(line_fields) == 3) and line_fields[2] == 'gene':
                        gene_start = line_fields[0]
                        gene_end = line_fields[1]
                        b_next_is_gene = True
                        b_next_is_product = False
                    # -----------------------------------------------------------------------------------
                    else:
                        print("line_fields >= 3 ("+str(len(line_fields))+"):"+','.join(line_fields))
                        sys.exit("Case not encountered line "+ str(sys._getframe().f_lineno) )
                        
                elif b_next_is_protein_id:
                    if line_fields[0] == 'protein_id':
                        protein_id = line_fields[1]
                        b_next_is_protein_id = False
                        print(' '.join([curr_type,
                                        cds_start,
                                        cds_end,
                                        f"{product} ({protein_id})" ]))
                        # store info
                        chrs.append(contig)
                        names.append(f"{product} ({protein_id})")
                        types.append('cds')
                        starts.append(int(cds_start) + contig_pos_shift)
                        ends.append(int(cds_end) + contig_pos_shift)
                        starts_vardict.append(int(cds_start))
                        ends_vardict.append(int(cds_end))
                    elif line_fields[0] == 'product':
                        product = product + ' ' + ' '.join(line_fields)
                    elif line_fields[0] == 'exception':
                        # add info between brackets
                        product = product + ' ['+' '.join(line_fields)+']'
                    else:
                        print("expected protein_id")
                        print("line_fields:"+','.join(line_fields))                        
                        sys.exit("Case not encountered line "+ str(sys._getframe().f_lineno) )

                elif b_next_is_note:
                    if line_fields[0] == 'note':
                        note = ' '.join(line_fields[1:])
                        b_next_is_note = False
                        print(' '.join(['misc_feature',
                                        misc_feature_start,
                                        misc_feature_end,
                                        note]))
                        # store info
                        chrs.append(contig)
                        names.append(note)
                        types.append('misc_feature')
                        starts.append(int(misc_feature_start)  + contig_pos_shift)
                        ends.append(int(misc_feature_end)  + contig_pos_shift)
                        starts_vardict.append(int(misc_feature_start))
                        ends_vardict.append(int(misc_feature_end))
                        
                    else: # means new start stop for misc_features
                        misc_feature_start = line_fields[0]
                        misc_feature_end = line_fields[1]
                        # print(f"next still note for line '{line}'")
                        # print(f"lien_fields 0:{line_fields[0]}")
                        # print(f"lien_fields 1:{line_fields[1]}")
                        # print(f"lien_fields 2:{line_fields[2]}")
                        # sys.exit("NOTE STOP")

                elif b_additional:
                    # treatment of exception (stop codon in frame)
                    if 'CDS_HAS_STOP_CODON' in line:
                        if b_verbose:
                            print(f"STOP in CDS: {line}")
                            
                        # get corrected end coordonates
                        m = re.search(r'\(CDS:([^\)]+)\) .*? seq-coords:(\d+)', line)
                        # first match: CDS name
                        cds_name_2correct = m.group(1)
                        # second match: position just before ORF end
                        cds_stop_corrected = int(m.group(2)) - 1  + contig_pos_shift 
                        cds_stop_corrected_str = str( cds_stop_corrected )
                        if b_verbose:
                            print(f"we will correct '{cds_name_2correct}' with new end {cds_stop_corrected_str}")
                            
                        # get name index of the gene to modify
                        try:
                           # if found, it is a gene
                           index2correct = names.index(cds_name_2correct)
                           print(f"'{cds_name_2correct}' with new end {cds_stop_corrected_str} (replace {ends[index2correct]})")
                           ori_end = ends[index2correct]
                           ends[index2correct] = cds_stop_corrected
                        except ValueError:
                            # get name index of the misc_feature to modify
                            try:
                                # if found, it is misc_feature
                                index2correct = names.index(f"similar to {cds_name_2correct}")
                                print(f"'similar to {cds_name_2correct}' with new end {cds_stop_corrected_str} (replace {ends[index2correct]})")
                                ori_end = ends[index2correct]
                                ends[index2correct] = cds_stop_corrected
                            except ValueError:
                                sys.exit(f"'{cds_name_2correct}' or 'similar to {cds_name_2correct}' not found in names to correct end position of CDS with STOP codon into, line "+str(sys._getframe().f_lineno))
                                
                        try:
                            print(f"search for start:{starts[index2correct]}, line "+str(sys._getframe().f_lineno) )
                            # if found, it is a misc_feature, we search first a gene with the same limits to modify too
                            # index_gene_start2correct = starts.index(starts[index2correct])
                            indexes_of_similar_starts = indexlist(starts[index2correct], # item2find
                                                                  starts)                # list_or_string
                            print("indexes of similar starts:"+str(indexes_of_similar_starts)+", line "+str(sys._getframe().f_lineno) )
                            # index_gene_end2correct = ends.index(ends[index2correct])
                            
                            for index2check in indexes_of_similar_starts:
                                # check if end is also the same for the given index (compared to value of index to correct)
                                if( (ori_end <= ends[index2check])and
                                    (types[index2check] in material_type_list) ):
                                    print(f"'{names[index2check]}' with new end {cds_stop_corrected} (replace {ends[index2check]})")
                                    ends[index2check] = cds_stop_corrected
                                    corrected_name = names[index2check]
                                else:
                                    print(f"{prog_tag} [Warn] material_type ({types[index2check]}) is not in {material_type} or ori_end ({ori_end} > found end ({ends[index2check]})) line "+str(sys._getframe().f_lineno) )
                                       
                        except ValueError:
                            warnings.warn(f"Warn '{cds_name_2correct}' has no {material_type} with the same start ({starts[index2correct]}) and end ({starts[index2correct]}), line "+ str(sys._getframe().f_lineno) )

                            # try:
                            #    # then we modify the original misc_feature too
                            #    index2correct = names.index(f"similar to {cds_name_2correct}")
                            #    print(f"'similar to {cds_name_2correct}' with new end {cds_stop_corrected} (replace {ends[index2correct]})")                               
                            #    ends[index2correct] = cds_stop_corrected
                            # except ValueError:
                            #    sys.exit(f"Error 'similar to {cds_name_2correct}' not found in names")
                # -------------------------------------------------------------

                # elif re.search(r"^>[A-Za-z]", line_fields[0]):
                if re.search(r">Feature", line_fields[0]):                                        
                    b_additional = False
                    print(f"Treating {line_fields[0]} contig")
                    m = re.search(r">Feature ([A-Za-z0-9\.]+)", line)
                    contig = m.group(1)
                    # print(f"contig found:{contig} for line_fields 0:{line_fields[0]} in line {line}")
                    # get the index of current contig
                    contig_index = contig_names.index(contig)
                    # if not the first contig, we must shift all start end position recorded according
                    # to the sum of all previous found contigs
                    if contig_index != 0:
                        contig_pos_shift += int(contig_lengths[contig_index - 1])
                    
                elif line_fields[2] == 'gene':
                    gene_start = line_fields[0]
                    gene_end   = line_fields[1]
                    gene_start = re.sub(non_alphanum, '', gene_start)
                    gene_end   = re.sub(non_alphanum, '', gene_end)
                    b_next_is_gene = True

                elif line_fields[2] == 'CDS':
                    cds_start = line_fields[0]
                    cds_end = line_fields[1]
                    cds_start = re.sub(non_alphanum, '', cds_start)
                    cds_end   = re.sub(non_alphanum, '', cds_end)
                    b_next_is_product = True
                    curr_type = 'CDS'

                elif line_fields[2] == 'mat_peptide':
                    cds_start = line_fields[0]
                    cds_end = line_fields[1]
                    cds_start = re.sub(non_alphanum, '', cds_start)
                    cds_end   = re.sub(non_alphanum, '', cds_end)
                    b_next_is_product = True    
                    curr_type = 'mat_peptide'
                    
                elif line_fields[2] == 'misc_feature':
                    misc_feature_start = line_fields[0]
                    misc_feature_end = line_fields[1]
                    misc_feature_start = re.sub(non_alphanum, '', misc_feature_start)
                    misc_feature_end   = re.sub(non_alphanum, '', misc_feature_end)
                    b_next_is_note = True

                elif "Additional" in line_fields[0]:
                    b_additional = True

            except IndexError:
                pass # means no other line

    paf.close()

# create json output from recorded info

# example of expected json format:
# {"genomesize": 169296, "virus_id": "SvRSV-RNA2-I2_MH883760.1", "BNRF1 tegument protein": [1736, 5692], "BCRF0": [5851, 6369], "BCRF1": [9680, 10192], "BHRF1": [41939, 42514], "BFRF1A": [46092, 46499], "BFRF1": [46459, 47469], "BFRF2": [47376, 49151], "BFRF3": [49075, 49605], "BORF1": [62792, 63886], "BORF2": [63961, 66441], "BaRF1": [66454, 67362], "BMRF1": [67453, 68667], "BMRF2": [68672, 69745], "BLRF1": [76099, 76407], "BRRF1": [92636, 93568], "BKRF2": [97412, 97825], "BKRF3": [97807, 98574], "BKRF4": [98588, 99241], "BBRF1": [101658, 103499], "BBRF2": [103402, 104238], "BBRF3": [106591, 107808], "BcRF1": [124915, 127167], "BTRF1": [127154, 128368], "BXRF1": [132309, 133055], "BVRF2": [135376, 137193], "BARF1": [164720, 165385], "UTR5": [1, 1735], "UTR3": [165386, 169296]}

## PUT LATER after gene selection in lists (en AFTER sorting: verify no CREATED PB)
# print(f"creates {json_annot_f} file")
# with open(json_annot_f, 'w+') as f:
#     f.write("{\"genomesize\": "+ genome_length +", ")
#     if len(chrs) == 1:
#         f.write(f"\"virus_id\": \"{contig}\", ")
#     else:
#         f.write(f"\"virus_id\": \"unkown\", ")
#     b_first_cds_found = False
#     for i in range(len(names)):
#         if (b_first_cds_found) and (types[i] == 'gene'):
#             f.write(", ")
#         if types[i] == 'gene':
#             f.write(f"\"{names[i]}\": [{starts[i]}, {ends[i]}]")
#             b_first_cds_found = True
#     f.write("}")
# f.close()

# print(f"len starts:{len(starts)}")
# print(f"len ends  :{  len(ends)}")
# print(f"len chrs  :{  len(chrs)}")
# print(f"len names :{ len(names)}")

# -------------------------------------------------------------
# keep only gene to avoid redundant start and stop (and missorting)
genes_indices =  [index for (index, item) in enumerate(types) if item == "gene"]
starts =         [ starts[i]         for i in genes_indices ]
ends   =         [ ends[i]           for i in genes_indices ]
starts_vardict = [ starts_vardict[i] for i in genes_indices ]
ends_vardict   = [ ends_vardict[i]   for i in genes_indices ]
names  =         [ names[i]          for i in genes_indices ]
chrs   =         [ chrs[i]           for i in genes_indices ]
# -------------------------------------------------------------

# sort the list of gene according to their start positions
starts, ends, starts_vardict, ends_vardict, names, chrs = (list(t) for t in zip(*sorted(zip(starts, ends, starts_vardict, ends_vardict, names, chrs))))


print(f"creates {json_annot_f} file")
with open(json_annot_f, 'w+') as f:
    f.write("{\"genomesize\": "+ genome_length +", ")
    if len(chrs) == 1:
        f.write(f"\"virus_id\": \"{contig}\", ")
    else:
        f.write(f"\"virus_id\": \"unkown\", ")
    b_first_cds_found = False
    for i in range(len(names)):
        if b_first_cds_found:
            f.write(", ")
 
        f.write(f"\"{names[i]}\": [{starts[i]}, {ends[i]}]")
        b_first_cds_found = True
    f.write("}")
f.close()



# # no used because keep original contig start position for vardict and real position used
# # are obtained from json
# print(f"creates {bed_annot_f} file")
# with open(bed_annot_f, 'w+') as f:
#     for i in range(len(names)):
#             f.write(f"{chrs[i]}\t{starts[i]}\t{ends[i]}\t{names[i]}\n")
# f.close()

print(f"creates {bed_vardict_annot_f} file")
with open(bed_vardict_annot_f, 'w+') as f:
    # no header otherwise bug in vardict
    for i in range(len(names)):
        f.write(f"{chrs[i]}\t{starts_vardict[i]}\t{ends_vardict[i]}\t{names[i]}\n")
f.close()
