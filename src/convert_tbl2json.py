#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of the vvv2_display distribution (https://github.com/ANSES-Ploufragan/vvv2_display).
# Copyright (c) 2023 Fabrice Touzain.
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
b_check_gene_prot_rec = False
# to have protein id found for match (default: deactivated because takes a lot of place)
b_include_protein_id = False
b_display_not_treated_cases = False

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
        pass_annot_f = test_dir + resn +".vadr.pass.tbl"
        fail_annot_f = test_dir + resn +".vadr.fail.tbl"
        json_annot_f = test_dir + resn +".vadr.json"
#        bed_annot_f  = test_dir + resn +".vadr.bed"
        bed_vardict_annot_f  = test_dir + resn +".vadr.4vardict.bed"                
        seq_stat_f   = test_dir + resn +".vadr.seqstat"
        cmd = ' '.join(['./convert_tbl2json.py',
                        "-p "+pass_annot_f,
                        "-f "+fail_annot_f,
                        "-s "+seq_stat_f,                    
                        "-j "+json_annot_f,
#                       "-b "+bed_annot_f,
                        "-c "+bed_vardict_annot_f                        
                        ])
        if b_verbose:
            cmd = cmd+" -v"
        print("cmd:"+cmd)
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
                print("genome_length: "+genome_length)
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
                print("contig "+contig_name+" length: "+str(contig_length))

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
cpt_gene = 0

# pattern to remove non alphanumeric character (ie: <> character befire start end positions of genes, cds, etc...)
non_alphanum = re.compile(r"\W") # [\'\(\)\W]")

# for each annotation file created by vadr
for annot_f in [pass_annot_f, fail_annot_f]:

    # initialize at the beginning of the file to avoid to keep variables of previous file (vadr pass when reading fail)
    b_next_is_gene = False
    b_next_is_product = False
    b_next_is_protein_id = False
    b_next_is_note =  False

    # read annotation file
    print(prog_tag + " read "+annot_f+" file")
    with open(annot_f) as paf:
        for line in paf:
            try:
                line = line.strip()
                # print("strip passed")
                line_fields = re.split(r"[\t\s]+", line) # PREVIOUS
                # line_fields = re.split(r"[\t]+", line)                
                # print("split passed")

                if line == '':
                    continue
                
                if b_verbose:
                    print("line_fields:"+ str(line_fields)+", line "+str(frame.f_lineno) + "\n\n")                

                # if re.search(r"gene", line):
                #     print("\t".join([
                #         "GENE ON LINE",
                #         "b_next_is_gene      :"+str(b_next_is_gene),
                #         "b_next_is_product   :"+str(b_next_is_product),
                #         "b_next_is_protein_id:"+str(b_next_is_protein_id),
                #         "b_next_is_note      :"+str(b_next_is_note)
                #     ]))

                # -------------------------------------------------------------
                # treat 'next' lines when something is anounced the line before
                if b_next_is_gene:
                    if b_verbose:
                        print(prog_tag + " NXT_is_GENE "+str(b_next_is_gene)+", line "+str(frame.f_lineno))

                    if line_fields[0] == 'gene':

                        # record gene info: correct previous recorded gene name (gene_n)
                        if re.search(r'^gene_', gene_name) or gene_name == '': 
                            gene_name = ' '.join(line_fields[1:])
                        elif gene_name == '':
                            gene_name = ' '.join(line_fields[1:])
                        else:
                            gene_name = gene_name + ' ' + ' '.join(line_fields[1:])

                        # print(' '.join(['gene',
                        #                gene_start,
                        #                gene_end,
                        #                gene_name]))

                       # store info
                        chrs.append(contig)
                        names.append(gene_name)
                        types.append('gene')
                        starts.append(int(gene_start) + contig_pos_shift)
                        ends.append(int(gene_end) + contig_pos_shift)
                        starts_vardict.append(int(gene_start))
                        ends_vardict.append(int(gene_end))
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag,
                                "RECORD: name:" + gene_name,
                                "type:gene",
                                "starts:"+gene_start,
                                "end:"+gene_end,
                                "contig:"+contig+", line "+str(frame.f_lineno)
                            ]))

                        b_next_is_gene = False                            
                        b_next_is_product = False
                        curr_type = 'gene'
                        gene_name = ''
                        gene_start = ''
                        gene_end = ''

                        if len(names) != len(types):
                            sys.exit(prog_tag+"[Error] names len:"+str(len(names))+" != types len:"+str(len(types)))

                        continue
                    
                    elif re.search(r'^>Feature', line_fields[0]):
                        b_additional = False
                        m = re.search(r">Feature (\S+)", line)
                        contig = m.group(1)
                        if b_verbose:
                            print(prog_tag + " Treating contig "+contig)

                        # get the index of current contig
                        contig_index = contig_names.index(contig)
                        # if not the first contig, we must shift all start end position recorded according
                        # to the sum of all previous found contigs
                        if contig_index != 0:
                            contig_pos_shift += int(contig_lengths[contig_index - 1])
                        continue

                    elif line_fields[2] == 'CDS':

                        # # added 2024 10 01, getting info on the line is more sure, previous gene info not always available
                        # if gene_start == '': # get gene pos only if not already buffered
                        #     gene_start = re.sub(non_alphanum, '', line_fields[0])
                        #     gene_end   = re.sub(non_alphanum, '', line_fields[1])
                        cds_start = re.sub(non_alphanum, '', line_fields[0])
                        cds_end   = re.sub(non_alphanum, '', line_fields[1])

                        # ----------------------------------------------------------
                        # check if previous gene is the same with a better name
                        last_gene_index = None
                        if(len(starts) > 0):
                            for i in reversed( range(len(starts)) ):
                                if types[i] == 'gene':
                                    last_gene_index = i
                                    break

                        if( (last_gene_index is not None) and 
                            (gene_start != "") and
                            (starts[last_gene_index] == int(gene_start) + contig_pos_shift ) and # same start pos
                            (ends[last_gene_index]   == int(gene_end  ) + contig_pos_shift ) # same start pos
                        ):
                            if not re.match('gene_', gene_name):
                                if re.match('gene_', names[last_gene_index]):
                                    # if so, we replace previously recorded name by current one
                                    names[i] = gene_name
                                elif re.match('similar to ', names[last_gene_index]):
                                    # if so, we replace previously recorded name by current one
                                    names[i] = gene_name + ' ' + names[last_gene_index]

                            if b_verbose or b_check_gene_prot_rec:
                                # if so, no need to record again
                                print("\t".join([
                                    prog_tag,
                                    "PASS (already RECORDED with "+names[last_gene_index]+")",
                                    "type:gene",
                                    "starts:"+gene_start,
                                    "end:"+gene_end,
                                    "contig:"+contig+", line "+str(frame.f_lineno)
                                ]))                       
                                
                        else:
                            # ----------------------------------------------------------
                            # record gene info only if found previously
                            if gene_start != '':
                                cpt_gene += 1
                                gene_name = 'gene_'+str(cpt_gene)
                                #    types.append('gene')

                                # print("CDS for line "+line)

                                if b_verbose:
                                    print(' '.join(['gene',
                                                    gene_start,
                                                    gene_end,
                                                    gene_name])
                                        )
                                
                                # store info of previous line \d+ \d+, gene
                                chrs.append(contig)
                                names.append(gene_name)
                                types.append('gene')
                                # types.append('cds') # changed 2024 10 08
                                starts.append(int(gene_start) + contig_pos_shift)
                                ends.append(int(gene_end) + contig_pos_shift)
                                starts_vardict.append(int(gene_start))
                                ends_vardict.append(int(gene_end))
                                if b_verbose or b_check_gene_prot_rec:
                                    print("\t".join([prog_tag,
                                                    "RECORD: name:"+gene_name,
                                                    "type:"+types[len(types)-1],
                                                    "starts:"+gene_start,
                                                    "end:"+gene_end,
                                                    "contig:"+contig+ ", line "+str(frame.f_lineno)]))
                            
                        if len(names) != len(types):
                            sys.exit(prog_tag + "[Error] names len:"+str(len(names))+" != types len:"+str(len(types)))
                        
                        b_next_is_gene = False
                        b_next_is_product = True
                        curr_type = 'CDS'
                        gene_name = ''
                        gene_start = ''
                        gene_end = ''
                        continue

                    # NEW
                    elif line_fields[2] == 'mat_peptide':
                        if b_verbose:
                            print(' '.join(['mat_peptide',
                                        gene_start,
                                        gene_end,
                                        gene_name])
                                )
                        b_next_is_gene = False
                        b_next_is_product = True                        
                        
                        # ask to treat current mat_peptide info
                        cds_start = line_fields[0]
                        cds_end = line_fields[1]
                        cds_start = re.sub(non_alphanum, '', cds_start)
                        cds_end   = re.sub(non_alphanum, '', cds_end)
                        curr_type = 'mat_peptide'
                        gene_name = ''
                        continue
                    
                    # section added 2024 09 24, when 2 gene section are following each other
                    # without any other description
                    elif line_fields[2] == 'gene':

                        gene_start = re.sub(non_alphanum, '', line_fields[0])
                        gene_end   = re.sub(non_alphanum, '', line_fields[1])

                        # ----------------------------------------------------------
                        # check if previous gene is the same with a better name
                        if(len(starts) > 0):
                            last_gene_index = None
                            for i in reversed( range(len(starts)) ):
                                if types[i] == 'gene':
                                    last_gene_index = i
                                    break

                            if( (last_gene_index is not None) and 
                                (starts[last_gene_index] == int(gene_start) + contig_pos_shift ) and # same start pos
                                (ends[last_gene_index]   == int(gene_end  ) + contig_pos_shift ) # same start pos
                            ):
                                if not re.match('gene_', gene_name):
                                    if re.match('gene_', names[last_gene_index]):
                                        # if so, we replace previously recorded name by current one
                                        names[i] = gene_name
                                    elif re.match('similar to ', names[last_gene_index]):
                                        # if so, we replace previously recorded name by current one
                                        names[i] = gene_name + ' ' + names[last_gene_index]

                                if b_verbose or b_check_gene_prot_rec:
                                    # if so, no need to record again
                                    print("\t".join([
                                        prog_tag,
                                        "PASS (already RECORDED with "+names[last_gene_index]+"): name:" + gene_name,
                                        "type:gene",
                                        "starts:"+gene_start,
                                        "end:"+gene_end,
                                        "contig:"+contig+", line "+str(frame.f_lineno)
                                    ]))
                                b_next_is_gene = False                            
                                b_next_is_product = False
                                curr_type = 'gene'
                                gene_name = ''
                                gene_start = ''
                                gene_end = ''
                                
                                if len(names) != len(types):
                                    sys.exit(prog_tag+"[Error] names len:"+str(len(names))+" != types len:"+str(len(types)))
                                continue
                        # ----------------------------------------------------------

                        # record gene info
                        cpt_gene += 1
                        gene_name = 'gene_'+str(cpt_gene)

                        # # record gene info: correct previous recorded gene name (gene_n)
                        # if re.search(r'^gene_', gene_name) or gene_name == '': 
                        #     gene_name = ' '.join(line_fields[1:])
                        # elif gene_name == '':
                        #     gene_name = ' '.join(line_fields[1:])
                        # else:
                        #     gene_name = gene_name + ' ' + ' '.join(line_fields[1:])

                        # print(' '.join(['gene',
                        #                gene_start,
                        #                gene_end,
                        #                gene_name]))

                       # store info
                        chrs.append(contig)
                        names.append(gene_name)
                        types.append('gene')
                        starts.append(int(gene_start) + contig_pos_shift)
                        ends.append(int(gene_end) + contig_pos_shift)
                        starts_vardict.append(int(gene_start))
                        ends_vardict.append(int(gene_end))
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag,
                                "RECORD: name:" + gene_name,
                                "type:gene",
                                "starts:"+gene_start,
                                "end:"+gene_end,
                                "contig:"+contig+", line "+str(frame.f_lineno)
                            ]))

                        gene_start = line_fields[0]
                        gene_end   = line_fields[1]
                        gene_start = re.sub(non_alphanum, '', gene_start)
                        gene_end   = re.sub(non_alphanum, '', gene_end)
                        
                        if b_verbose:
                            print("\t".join([
                                "TMP_GENE:"+line_fields[2],
                                gene_start,
                                gene_end,
                                "for line "+str(line_fields)
                            ]))

                        b_next_is_gene = False                                                        
                        b_next_is_product = False
                        
                        curr_type = 'gene'
                        gene_name = ''

                        if len(names) != len(types):
                            sys.exit(prog_tag+"[Error] names len:"+str(len(names))+" != types len:"+str(len(types)))

                        continue
                    # -----------------------------------------------------------------------------------
                    # added 2024 09 26 because misc_feature found after note (therefore b_next_gene is True)
                    elif (line_fields[2] == 'misc_feature'):
                        misc_feature_start = re.sub(non_alphanum, '', line_fields[0])
                        misc_feature_end   = re.sub(non_alphanum, '', line_fields[1])

                        if b_next_is_gene: # means gene recorded just above, if same coordinate
                                           # we need to record gene
                            if( (gene_start == misc_feature_start)and
                                (gene_end   == misc_feature_end) ):
                            
                                # record gene info
                                cpt_gene += 1
                                gene_name = 'gene_'+str(cpt_gene)

                                # # record gene info: correct previous recorded gene name (gene_n)
                                # if re.search(r'^gene_', gene_name) or gene_name == '': 
                                #     gene_name = ' '.join(line_fields[1:])
                                # elif gene_name == '':
                                #     gene_name = ' '.join(line_fields[1:])
                                # else:
                                #     gene_name = gene_name + ' ' + ' '.join(line_fields[1:])

                                # print(' '.join(['gene',
                                #                gene_start,
                                #                gene_end,
                                #                gene_name]))

                                # store info
                                chrs.append(contig)
                                names.append(gene_name)
                                types.append('gene')
                                starts.append(int(gene_start) + contig_pos_shift)
                                ends.append(int(gene_end) + contig_pos_shift)
                                starts_vardict.append(int(gene_start))
                                ends_vardict.append(int(gene_end))
                                if b_verbose or b_check_gene_prot_rec:
                                    print("\t".join([
                                        prog_tag,
                                        "RECORD: name:" + gene_name,
                                        "type:gene",
                                        "starts:"+gene_start,
                                        "end:"+gene_end,
                                        "contig:"+contig+", line "+str(frame.f_lineno)
                                    ]))

                        else:
                            if b_verbose or b_check_gene_prot_rec:
                                print("\t".join([
                                    prog_tag,
                                    "MISC_FEATURE:",
                                    "starts:"+misc_feature_start,
                                    "end:"+misc_feature_end,
                                    "contig:"+contig + ", line "+str(frame.f_lineno)
                                ]))
                        
                        b_next_is_gene     = False # added 2024 10 01
                        b_next_is_note     = True
                        b_next_is_product  = False
                        continue
                    # added 2024 11 12 handle stem_loop after protein, expecting gene
                    elif line_fields[2] == 'stem_loop':
                        misc_feature_start = line_fields[0]
                        misc_feature_end   = line_fields[1]
                        misc_feature_start = re.sub(non_alphanum, '', misc_feature_start)
                        misc_feature_end   = re.sub(non_alphanum, '', misc_feature_end)

                        b_next_is_product = False
                        b_next_is_gene = False
                        b_next_is_note = True                    
                        curr_type = 'stem_loop'
                        if b_verbose:
                            print("\t".join([
                                "TMP_STEM_LOOP:"+line_fields[2],
                                misc_feature_start,
                                misc_feature_end
                            ]))
                    # added 2024 09 26 handle note after misc_feature
                    elif line_fields[0] == 'note':
                        note = ' '.join(line_fields[1:])
                        if b_verbose:
                            print("initial note: "+note+", line "+ str(sys._getframe().f_lineno) )    

                        # print("line_fields 1..3:"+str(line_fields[1:3]))
                        
                        # -----------------------------                        
                        # means that we face a protein similarity we must display as a protein / cds
                        if((' '.join(line_fields[1:3])) == 'similar to'):
                        # if(re.search(r'^similar to',line_fields[1] ):                            

                            m = re.search(r'(similar to [^\t]+)', line)
                            note = re.sub('\'','', m.group(1))

                            if b_verbose:
                                print("note set to "+note+", line "+ str(sys._getframe().f_lineno) )
                                print(' '.join(['cds',
                                                misc_feature_start,
                                                misc_feature_end,
                                                note]))
                                print("found "+note+", line "+ str(sys._getframe().f_lineno) )    
                            types.append('cds')
                            # types.append(curr_type)                            
                        # -----------------------------
                        else: # otherwise, records classical misc_feature

                            if b_verbose:
                                print(' '.join(['misc_feature',
                                                misc_feature_start,
                                                misc_feature_end,
                                                note]))
                            # types.append('misc_feature')
                            types.append(curr_type)                   
                            if b_verbose:         
                                print("found misc_feature "+note+", line "+ str(sys._getframe().f_lineno) )
                                    
                        # store info
                        chrs.append(contig)
                        names.append(note)
                        # type recorded above
                        starts.append(int(misc_feature_start)  + contig_pos_shift)
                        ends.append(int(misc_feature_end)  + contig_pos_shift)
                        starts_vardict.append(int(misc_feature_start))
                        ends_vardict.append(int(misc_feature_end))
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag, 
                                "RECORD: name:"+note,
                                "type:"+types[-1],
                                "starts:"+misc_feature_start,
                                "end:"+misc_feature_end,
                                "contig:"+contig+", line "+str(frame.f_lineno)
                            ]))
                        note = ''
                        b_next_is_note = False
                        b_next_is_gene = False
                        
                        if len(names) != len(types):
                            sys.exit("names and types not with the same number of elements, line "+str(len(types)))
                        continue          

                    # added 2024 10 01 because in some cases, we can find protein_id when gene is expected
                    elif line_fields[0] == 'protein_id':
                        protein_id = line_fields[1]
                        tmp_name = protein_id
                        b_next_is_protein_id = False

                        # escape cotes, because misinterprated by R
                        product    = re.sub(r"([\'])", '', product)
                        protein_id = re.sub(r"([\'])", '', protein_id)
                        # remove spaces at beginning and end
                        product    = product.strip()
                        protein_id = protein_id.strip()                        
                        
                        # store info
                        chrs.append(contig)
                        if b_include_protein_id and (product != ''):
                            tmp_name = product + " " +protein_id
                        elif product == '': # in some case, only protein_id is provided, we must keep in this case
                            # even if b_include_protein_id is False
                            tmp_name = protein_id
                        else:
                            tmp_name = product
                        names.append(tmp_name)

                        # print(' '.join([curr_type,
                        #                 cds_start,
                        #                 cds_end,
                        #                 tmp_name ]))
                            
                        types.append('cds')
                        starts.append(int(cds_start) + contig_pos_shift)
                        ends.append(int(cds_end) + contig_pos_shift)
                        starts_vardict.append(int(cds_start))
                        ends_vardict.append(int(cds_end))
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag,
                                "RECORD: name:"+tmp_name,
                                "type:cds",
                                "starts:"+gene_start,
                                "end:"+gene_end+", line "+str(frame.f_lineno)
                            ]))
                        product = ''
                        protein_id = ''
                        b_next_is_protein_id = False
                        continue

                    # added 2024 09 26 handle "Additional information", no gene record afterthat, only possible corrections
                    # on start/stop info
                    elif re.match(r"^Additional", line_fields[0]) :
                        b_next_gene          = False
                        b_next_is_note       = False
                        b_next_is_product    = False
                        b_next_is_protein_id = False 
                        b_additional         = True
                        if b_verbose:
                            print("Additional section detected")
                        continue
                    elif not b_additional:

                        # print("No CDS field2 or gene field0, line "+str(frame.f_lineno))
                        gene_name = line_fields[1]

                        # print(' '.join(['gene',
                        #                gene_start,
                        #                gene_end,
                        #                gene_name]))
                        b_next_is_gene = False
                        # store info
                        chrs.append(contig)
                        names.append(gene_name)
                        types.append('gene')
                        starts.append(int(gene_start) + contig_pos_shift)
                        ends.append(int(gene_end) + contig_pos_shift)
                        starts_vardict.append(int(gene_start))
                        ends_vardict.append(int(gene_end))
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag,
                                "RECORD: name:{gene_name}",
                                "type:gene",
                                "starts:"+gene_start,
                                "end:"+gene_end,
                                "contig:"+contig + ", line "+str(frame.f_lineno)
                            ]))
                        if len(names) != len(types):
                            sys.exit(prog_tag +"[Error] names len:"+str(len(names))+" != types len:"+str(len(types)))
                        gene_name  = ''
                        gene_start = ''
                        gene_end   = ''
                        continue
                    
                elif b_next_is_product:
                    if b_verbose:
                        print(prog_tag + " NXT is PRODUCT "+str(b_next_is_product)+", line "+str(frame.f_lineno))
                    
                    if re.search(r'^product', line_fields[0]):
                        product = ' '.join(line_fields[0].split(' ')[1:]) + ' ' + ' '.join(line_fields[1:])
                        b_next_is_product = False
                        b_next_is_protein_id = True
                        if b_verbose:
                            print("TMP_PRODUCT:"+product+", line "+str(frame.f_lineno))
                        continue
                    elif len(line_fields) == 2:
                        # means only coordinates because gene longer than CDS
                        if b_verbose:
                            print(prog_tag + " CDS shorter than gene, skip useless coordinates")
                        continue
                    # -----------------------------------------------------------------------------------
                    # this part until else MUST be useless, but bug found only when b_verbose is True...
                    elif (len(line_fields) == 3) and (line_fields[2] == 'misc_feature'):
                        misc_feature_start = re.sub(non_alphanum, '', line_fields[0])
                        misc_feature_end   = re.sub(non_alphanum, '',bline_fields[1])
                        b_next_is_note     = True
                        b_next_is_product  = False
                        continue
                    elif (len(line_fields) == 3) and (line_fields[2] == 'gene'):
                        gene_start = re.sub(non_alphanum, '', line_fields[0])
                        gene_end   = re.sub(non_aplhanum, '', line_fields[1])
                        b_next_is_gene    = True
                        b_next_is_product = False
                        continue
                    elif (len(line_fields) == 2) and (line_fields[0] == 'gene'):

                        gene_name =  ' '.join( line_fields[0].split(' ')[1:] )+ ' '+ ' '.join(line_fields[1:])

                        if types[-1] == 'gene':
                            names[-1] = gene_name
                        else:
                            sys.exit(prog_tag + " [Error] found 'gene genename' supposed to correct recorded name at previous line, but previous type is "+types[-1]+", not normal, line "+str(frame.f_lineno))
                                                
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag, 
                                "CORRECT: name:" + gene_name,
                                "type:gene",
                                "starts:"+gene_start,
                                "end:"+gene_end,
                                "contig:"+contig+", line "+str(frame.f_lineno)
                            ]))
                        if len(names) != len(types):
                            sys.exit(prog_tag+"[Error] names len:"+str(len(names))+" != types len:"+str(len(types)))
                        
                        curr_type = 'gene'
                        b_next_is_product = False
                        b_next_is_gene = False
                        continue
                    
                    # -----------------------------------------------------------------------------------
                    elif (len(line_fields) == 3) and (line_fields[2] == 'CDS'):
                        cds_start = re.sub(non_alphanum, '', line_fields[0])
                        cds_end   = re.sub(non_alphanum, '', line_fields[1])
                        curr_type = 'CDS'
                        b_next_is_gene = False
                        b_next_is_product = True
                        continue
                    else:
                        print("line_fields:"+str(line_fields)+", line "+ str(sys._getframe().f_lineno) )
                        sys.exit(prog_tag + "[Error] Case not encountered line "+ str(sys._getframe().f_lineno) )
                        
                elif b_next_is_protein_id:
                    if b_verbose:
                        print(prog_tag + " NXT is PROTEIN "+str(b_next_is_protein_id)+", line "+str(frame.f_lineno))                    
                    if line_fields[0] == 'protein_id':
                        protein_id = line_fields[1]
                        tmp_name = protein_id
                        
                        # escape cotes, because misinterprated by R
                        product    = re.sub(r"([\'])", '', product)
                        protein_id = re.sub(r"([\'])", '', protein_id)
                        # remove spaces at beginning and end
                        product    = product.strip()
                        protein_id = protein_id.strip()                        
                        
                        # check if cds already recorded
                        # if yes, keep more relevant name, do not add new cds
                                               # ----------------------------------------------------------
                        # check if previous gene is the same with a better name
                        last_cds_index = None
                        if(len(starts) > 0):
                            for i in reversed( range(len(starts)) ):
                                if types[i] == 'cds':
                                    last_cds_index = i
                                    break

                        if( (last_cds_index is not None) and 
                            (starts[last_cds_index] == int(cds_start) + contig_pos_shift ) and # same start pos
                            (ends[last_cds_index]   == int(cds_end  ) + contig_pos_shift ) # same start pos
                        ):
                            if not re.match('gene_', product):
                                if re.match('gene_', names[last_cds_index]):
                                    # if so, we replace previously recorded name by current one
                                    names[i] = product
                                elif re.match('similar to ', names[last_cds_index]):
                                    # if so, we replace previously recorded name by current one
                                    names[i] = product + ' ' + names[last_cds_index]

                            if b_verbose or b_check_gene_prot_rec:
                                # if so, no need to record again
                                print("\t".join([
                                    prog_tag,
                                    "PASS (already RECORDED with "+names[last_cds_index]+")",
                                    "type:cds",
                                    "starts:"+cds_start,
                                    "end:"+cds_end,
                                    "contig:"+contig+", line "+str(frame.f_lineno)
                                ]))                       
                            product = ''
                            protein_id = ''
                            b_next_is_protein_id = False    
                            b_next_is_gene = True
                            continue
                        else:
                            # otherwise
                            # store info
                            chrs.append(contig)
                            if b_include_protein_id and (product != ''):
                                tmp_name = product + " " +protein_id
                            elif product == '': # in some case, only protein_id is provided, we must keep in this case
                                # even if b_include_protein_id is False
                                tmp_name = protein_id
                            else:
                                tmp_name = product
                            names.append(tmp_name)

                            # print(' '.join([curr_type,
                            #                 cds_start,
                            #                 cds_end,
                            #                 tmp_name ]))
                                
                            types.append('cds')
                            starts.append(int(cds_start) + contig_pos_shift)
                            ends.append(int(cds_end) + contig_pos_shift)
                            starts_vardict.append(int(cds_start))
                            ends_vardict.append(int(cds_end))
                            if b_verbose or b_check_gene_prot_rec:
                                print("\t".join([
                                    prog_tag,
                                    "RECORD: name:"+tmp_name,
                                    "type:cds",
                                    "starts:"+cds_start,
                                    "end:"+cds_end+", line "+str(frame.f_lineno)
                                ]))
                            product = ''
                            protein_id = ''
                            b_next_is_protein_id = False    
                            b_next_is_gene = True # added 2024 10 01                    
                            continue

                    if line_fields[0] == 'ncRNA_class':
                        protein_id = line_fields[1]
                        tmp_name = protein_id
                        b_next_is_protein_id = False

                        # escape cotes, because misinterprated by R
                        product    = re.sub(r"([\'])", '', product)
                        protein_id = re.sub(r"([\'])", '', protein_id)
                        # remove spaces at beginning and end
                        product    = product.strip()
                        protein_id = protein_id.strip()                        
                        
                        # store info
                        chrs.append(contig)
                        if b_include_protein_id and (product != ''):
                            tmp_name = product + " " +protein_id
                        elif product == '': # in some case, only protein_id is provided, we must keep in this case
                            # even if b_include_protein_id is False
                            tmp_name = protein_id
                        else:
                            tmp_name = product
                        names.append(tmp_name)

                        # print(' '.join([curr_type,
                        #                 cds_start,
                        #                 cds_end,
                        #                 tmp_name ]))
                            
                        # types.append('cds')
                        types.append('rna')                        
                        starts.append(int(cds_start) + contig_pos_shift)
                        ends.append(int(cds_end) + contig_pos_shift)
                        starts_vardict.append(int(cds_start))
                        ends_vardict.append(int(cds_end))
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag,
                                "RECORD: name:"+tmp_name,
                                "type:rna",
                                "starts:"+gene_start,
                                "end:"+gene_end+", line "+str(frame.f_lineno)
                            ]))
                        product = ''
                        protein_id = ''
                        continue
                    
                    elif re.match('product', line_fields[0]):
                        product_str = ' '.join(line_fields[0].split(' ')[1:]) + ' ' + ' '.join(line_fields[1:])
                        product = product + ' ' + product_str
                        b_next_is_protein_id = True
                        b_next_is_product    = False
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag,
                                "TMP_PRODUCT: product:"+product,
                                ", line "+str(frame.f_lineno)
                            ]))
                        continue
                    elif line_fields[0] == 'exception':
                        # add info between brackets
                        product = product + ' ['+' '.join(line_fields)+']'
                        continue
                    else:
                        print("expected protein_id")
                        print("line_fields:"+','.join(line_fields))                        
                        sys.exit("Case not encountered line "+ str(sys._getframe().f_lineno) )

                elif b_next_is_note:
                    if line_fields[0] == 'note':
                        note = ' '.join(line_fields[1:])
                        # print("initial note: "+note+", line "+ str(sys._getframe().f_lineno) )    
             
                        # -----------------------------                        
                        # means that we face a protein similarity we must display as a protein / cds
                        if((' '.join(line_fields[1:3])) == 'similar to'):
                        # if(re.search(r'^similar to',line_fields[1] ):                            

                            m = re.search(r'(similar to [^\t]+)', line)
                            note = re.sub('\'','', m.group(1))

                            if b_verbose:
                                print("note set to "+note+", line "+ str(sys._getframe().f_lineno) )
                                print(' '.join(['cds',
                                                misc_feature_start,
                                                misc_feature_end,
                                                note]))
                            types.append('cds')
                            # types.append(curr_type)                            
                        # -----------------------------
                        else: # otherwise, records classical misc_feature

                            if b_verbose:
                                print(' '.join(['misc_feature',
                                                misc_feature_start,
                                                misc_feature_end,
                                                note]))
                            # types.append('misc_feature')
                            types.append(curr_type)       
                            if b_verbose:                     
                                print("found misc_feature "+note+", line "+ str(sys._getframe().f_lineno) )
                                    
                        # store info
                        chrs.append(contig)
                        names.append(note)
                        # type recorded above
                        starts.append(int(misc_feature_start)  + contig_pos_shift)
                        ends.append(int(misc_feature_end)  + contig_pos_shift)
                        starts_vardict.append(int(misc_feature_start))
                        ends_vardict.append(int(misc_feature_end))
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag, 
                                "RECORD: name:"+note,
                                "type:"+types[-1],
                                "starts:"+misc_feature_start,
                                "end:"+misc_feature_end,
                                "contig:"+contig+", line "+str(frame.f_lineno)
                            ]))
                        note = ''
                        b_next_is_note = False
                        b_next_is_gene = True # added 2024 09 19
                        
                        if len(names) != len(types):
                            sys.exit("names and types not with the same number of elements, line "+str(len(types)))
                        continue                        
                        
                    else: # means new start stop for misc_features
                        misc_feature_start = line_fields[0]
                        misc_feature_end = line_fields[1]
                        # b_next_is_note = False                        
                        b_next_is_note = True # changed 2024 10 01                        
                        # sys.exit("NOTE STOP")
                        continue

                elif b_additional:
                    # treatment of exception (stop codon in frame)
                    if 'CDS_HAS_STOP_CODON' in line:
                        if b_verbose:
                            print("STOP in CDS: "+line)
                            
                        # get corrected end coordonates
                        m = re.search(r'\(CDS:([^\)]+)\) .*? seq-coords:(\d+)', line)
                        # first match: CDS name
                        cds_name_2correct = m.group(1)
                        # second match: position just before ORF end
                        cds_stop_corrected = int(m.group(2)) - 1  + contig_pos_shift 
                        cds_stop_corrected_str = str( cds_stop_corrected )
                        if b_verbose:
                            print(prog_tag + " We will correct "+cds_name_2correct+" with new end "+cds_stop_corrected_str)
                            
                        # get name index of the gene to modify
                        try:
                           # if found, it is a gene
                           index2correct = names.index(cds_name_2correct)

                           if b_verbose:
                               print(prog_tag + ' '+ cds_name_2correct+" with new end "+cds_stop_corrected_str+" (replace "+ends[index2correct])
                           ori_end = ends[index2correct]
                           ends[index2correct] = cds_stop_corrected
                        except ValueError:
                            # get name index of the misc_feature to modify
                            try:
                                # if found, it is misc_feature
                                index2correct = names.index("similar to "+cds_name_2correct)
                                if b_verbose:
                                    print("'similar to "+cds_name_2correct+"' with new end "+cds_stop_corrected_str+" (replace "+str(ends[index2correct])+")")
                                ori_end = ends[index2correct]
                                ends[index2correct] = cds_stop_corrected
                            except ValueError:
                                sys.exit("'"+cds_name_2correct+"' or 'similar to "+cds_name_2correct+"' not found in names to correct end position of CDS with STOP codon into, line "+str(sys._getframe().f_lineno))
                                
                        try:
                            if b_verbose:
                                print(prog_tag + " Search for start:"+str(starts[index2correct])+", line "+str(sys._getframe().f_lineno) )
                            # if found, it is a misc_feature, we search first a gene with the same limits to modify too
                            # index_gene_start2correct = starts.index(starts[index2correct])
                            indexes_of_similar_starts = indexlist(starts[index2correct], # item2find
                                                                  starts)                # list_or_string
                            if b_verbose:
                                print("indexes of similar starts:"+str(indexes_of_similar_starts)+", line "+str(sys._getframe().f_lineno) )
                            # index_gene_end2correct = ends.index(ends[index2correct])
                            
                            for index2check in indexes_of_similar_starts:
                                # check if end is also the same for the given index (compared to value of index to correct)
                                if( (ori_end <= ends[index2check])and
                                    (types[index2check] in material_type_list) ):
                                    if b_verbose:
                                        print("'"+names[index2check]+"' with new end "+str(cds_stop_corrected)+" (replace "+str(ends[index2check])+")")
                                    ends[index2check] = cds_stop_corrected
                                    corrected_name = names[index2check]
                                else:
                                    print(prog_tag+"[Warn] material_type ("+types[index2check]+") is not in "+material_type+" or ori_end ("+str(ori_end)+" > found end ("+str(ends[index2check])+"), line "+str(sys._getframe().f_lineno) )
                                    continue # can be found in others
                                       
                        except ValueError:
                            warnings.warn(prog_tag + "[Warn] '"+cds_name_2correct+"' has no "+material_type+" with the same start ("+starts[index2correct]+") and end ("+starts[index2correct]+"), line "+ str(sys._getframe().f_lineno) )

                            # try:
                            #    # then we modify the original misc_feature too
                            #    index2correct = names.index("similar to "+cds_name_2correct)
                            #    print(f"'similar to {cds_name_2correct}' with new end {cds_stop_corrected} (replace {ends[index2correct]})")                               
                            #    ends[index2correct] = cds_stop_corrected
                            # except ValueError:
                            #    sys.exit(f"Error 'similar to {cds_name_2correct}' not found in names")
                    

                    # treatment of exception (frameshift)
                    if 'FRAMESHIFT_HIGH_CONF' in line:
                        if b_verbose:
                            print("FRAMESHIFT_HIGH_CONF detected: "+line)
                            
                        # get corrected end coordonates
                        m = re.search(r'\(CDS:([^\)]+)\) .*? seq-coords:(\d+)\.\.(\d+)\:([\+\-]).*? mdl-coords:(\d+)\.\.(\d+)\:([\+\-])', line)
                        # first match: CDS name
                        cds_name_with_frameshift = m.group(1)
                        # second match: position of frameshift start
                        # third  match: position of frameshift end
                        fs_start = int(m.group(2)) - 1  + contig_pos_shift 
                        fs_end   = int(m.group(3)) - 1  + contig_pos_shift 
                        fs_strand=     m.group(4)
                        mdl_start = int(m.group(5)) - 1  + contig_pos_shift 
                        mdl_end   = int(m.group(6)) - 1  + contig_pos_shift 
                        mdl_strand=     m.group(7)
                        if b_verbose:
                            print(prog_tag + " Frameshift in "+cds_name_with_frameshift+" at pos "+fs_start+':'+fs_end)

                if b_verbose:
                    print(prog_tag + " line_fields:"+str(line_fields)+", line "+str(sys._getframe().f_lineno))
                if re.search(r'^>Feature', line_fields[0]):
                    b_additional = False
                    m = re.search(r">Feature (\S+)", line)
                    contig = m.group(1)
                    if b_verbose:
                        print(prog_tag + " Treating contig "+contig)

                    # get the index of current contig
                    contig_index = contig_names.index(contig)
                    # if not the first contig, we must shift all start end position recorded according
                    # to the sum of all previous found contigs
                    if contig_index != 0:
                        contig_pos_shift += int(contig_lengths[contig_index - 1])
                    continue

                # added 2024 09 24: now need to replace name of previous record, testing if start end are the same
                # otherwise record as set here
                elif line_fields[0] == 'gene':

                        
                        # record gene info: correct previous recorded gene name (gene_n)
                        if re.match(r'^gene_', gene_name) or gene_name == '': 
                            gene_name = ' '.join(line_fields[1:])
                            if b_verbose:
                                print("gene_name set to "+gene_name+", line "+str(frame.f_lineno))
                            
                        # ----------------------------------------------------------
                        # check if previous gene is the same with a better name
                        if(len(starts) > 0):
                            last_gene_index = None
                            for i in reversed( range(len(starts)) ):
                                if types[i] == 'gene':
                                    last_gene_index = i
                                    break

                            if( (last_gene_index is not None) and 
                                (starts[last_gene_index] == int(gene_start) + contig_pos_shift ) and # same start pos
                                (ends[last_gene_index]   == int(gene_end  ) + contig_pos_shift ) # same start pos
                            ):
                                if not re.match('gene_', gene_name):
                                    if re.match('gene_', names[last_gene_index]):
                                        # if so, we replace previously recorded name by current one
                                        names[last_gene_index] = gene_name
                                        if b_verbose:
                                            print("previous gene_name renamed to "+gene_name+", line "+str(frame.f_lineno))
                                    # elif names[last_gene_index] == '':
                                    #     # if so, we replace previously missing name by current one
                                    #     names[last_gene_index] = gene_name
                                    elif re.match('similar to ', names[last_gene_index]):
                                        # if so, we replace previously recorded name by current one
                                        names[last_gene_index] = gene_name + ' ' + names[last_gene_index]
                                        if b_verbose:
                                            print("previous gene_name renamed to "+gene_name+", line "+str(frame.f_lineno))

                                if b_verbose or b_check_gene_prot_rec:
                                    # if so, no need to record again
                                    print("\t".join([
                                        prog_tag,
                                        "PASS (already RECORDED new_name '"+names[last_gene_index]+"')",
                                        "type:gene",
                                        "starts:"+gene_start,
                                        "end:"+gene_end,
                                        "contig:"+contig+", line "+str(frame.f_lineno)
                                    ]))
                                b_next_is_gene = False                            
                                b_next_is_product = False
                                curr_type = 'gene'
                                gene_name = ''
                                gene_start = ''
                                gene_end = ''
                                if len(names) != len(types):
                                    sys.exit(prog_tag+"[Error] names len:"+str(len(names))+" != types len:"+str(len(types)))
                                continue
                        # ----------------------------------------------------------

                        else:
                            gene_name = gene_name + ' ' + ' '.join(line_fields[1:])

                        # print(' '.join(['gene',
                        #                gene_start,
                        #                gene_end,
                        #                gene_name]))

                       # store info
                        chrs.append(contig)
                        names.append(gene_name)
                        types.append('gene')
                        starts.append(int(gene_start) + contig_pos_shift)
                        ends.append(int(gene_end) + contig_pos_shift)
                        starts_vardict.append(int(gene_start))
                        ends_vardict.append(int(gene_end))
                        if b_verbose or b_check_gene_prot_rec:
                            print("\t".join([
                                prog_tag,
                                "RECORD: name:" + gene_name,
                                "type:gene",
                                "starts:"+gene_start,
                                "end:"+gene_end,
                                "contig:"+contig+", line "+str(frame.f_lineno)
                            ]))

                        b_next_is_gene = False                            
                        b_next_is_product = False
                        curr_type  = 'gene'
                        gene_name  = ''
                        gene_start = ''
                        gene_end   = ''

                        if len(names) != len(types):
                            sys.exit(prog_tag+"[Error] names len:"+str(len(names))+" != types len:"+str(len(types)))

                        continue

                elif line_fields[2] == 'gene':
                    gene_start = line_fields[0]
                    gene_end   = line_fields[1]
                    gene_start = re.sub(non_alphanum, '', gene_start)
                    gene_end   = re.sub(non_alphanum, '', gene_end)
                    b_next_is_gene = True
                    b_next_is_product = False
                    if b_verbose:
                        print("\t".join([
                            "TMP_GENE:"+line_fields[2],
                            gene_start,
                            gene_end,
                            "for line "+str(line_fields)
                        ]))

                elif line_fields[2] == 'CDS':
                    cds_start = line_fields[0]
                    cds_end   = line_fields[1]
                    cds_start = re.sub(non_alphanum, '', cds_start)
                    cds_end   = re.sub(non_alphanum, '', cds_end)
                    b_next_is_product = True
                    b_next_is_gene = False                    
                    curr_type = 'CDS'
                    if b_verbose:
                        print("\t".join([
                            "TMP_CDS:"+line_fields[2],
                            cds_start,
                            cds_end
                        ]))
                        
                elif line_fields[2] == 'ncRNA':
                    cds_start = line_fields[0]
                    cds_end   = line_fields[1]
                    cds_start = re.sub(non_alphanum, '', cds_start)
                    cds_end   = re.sub(non_alphanum, '', cds_end)
                    b_next_is_product = True
                    b_next_is_gene = False
                    curr_type = 'rna' # 'CDS'
                    if b_verbose:
                        print("\t".join([
                            "TMP_NC_RNA:"+line_fields[2],
                            cds_start,
                            cds_end
                        ]))

                elif line_fields[2] == 'misc_feature':
                    misc_feature_start = line_fields[0]
                    misc_feature_end   = line_fields[1]
                    misc_feature_start = re.sub(non_alphanum, '', misc_feature_start)
                    misc_feature_end   = re.sub(non_alphanum, '', misc_feature_end)
                    b_next_is_note = True
                    b_next_is_gene = False
                    if b_verbose:
                        print("\t".join([
                            "TMP_MISC_FEATURE:"+line_fields[2],
                            misc_feature_start,
                            misc_feature_end
                        ]))

                elif line_fields[2] == 'mat_peptide':
                    cds_start = line_fields[0]
                    cds_end   = line_fields[1]
                    cds_start = re.sub(non_alphanum, '', cds_start)
                    cds_end   = re.sub(non_alphanum, '', cds_end)
                    b_next_is_product = True
                    b_next_is_gene = False
                    curr_type = 'CDS'
                    if b_verbose:
                        print("\t".join([
                            "TMP_MAT_PEPTIDE:"+line_fields[2],
                            cds_start,
                            cds_end
                        ]))

                elif line_fields[2] == 'stem_loop':
                    misc_feature_start = line_fields[0]
                    misc_feature_end   = line_fields[1]
                    misc_feature_start = re.sub(non_alphanum, '', misc_feature_start)
                    misc_feature_end   = re.sub(non_alphanum, '', misc_feature_end)

                    b_next_is_product = False
                    b_next_is_gene = False
                    b_next_is_note = True                    
                    curr_type = 'stem_loop'
                    if b_verbose:
                        print("\t".join([
                            "TMP_STEM_LOOP:"+line_fields[2],
                            misc_feature_start,
                            misc_feature_end
                        ]))
                        
                elif "Additional" in line_fields[0]:
                    b_additional         = True
                    b_next_is_gene       = False
                    b_next_is_note       = False
                    b_next_is_product    = False
                    b_next_is_protein_id = False 

                else:
                    if b_display_not_treated_cases:
                        print(prog_tag + " Case not treated for line "+line+", line "+str(sys._getframe().f_lineno))

            except IndexError:
                print("Exception IndexError for line '"+line+"', line "+str(sys._getframe().f_lineno))
                raise # means no other line
            except:
                print("Exception for line "+line+", line "+str(sys._getframe().f_lineno))
                print("Unexpected error:", sys.exc_info()[0])
                raise
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

if b_verbose and b_check_gene_prot_rec:
    print(prog_tag+ " check genes/prot records BEFORE cds/gene selection--------------")
    for i in range(len(names)):
        print("name:"+names[i]+"\ttype:"+types[i]+" ["+str(starts[i])+", "+str(ends[i])+"]")
    print(prog_tag+ "---------------------------------------")
    
# -------------------------------------------------------------
# keep only gene / prot to avoid redundant start and stop (and missorting)

genes_indices  = []

starts_ori         = starts
ends_ori           = ends
starts_vardict_ori = starts_vardict
ends_vardict_ori   = ends_vardict
names_ori          = names
chrs_ori           = chrs
types_ori          = types

starts         = []
ends           = [] 
starts_vardict = [] 
ends_vardict   = [] 
names          = [] 
chrs           = []
types          = []

for i in range(len(types_ori)):
    if((types_ori[i] == "cds")or
       (types_ori[i] == "gene")or
       (types_ori[i] == "rna")or
       (types_ori[i] == "stem_loop")       
       ):
        genes_indices.append(i)
        starts.append(         starts_ori[ i ])
        ends.append(           ends_ori[ i ])
        starts_vardict.append( starts_vardict_ori[ i ])
        ends_vardict.append(   ends_vardict_ori[ i ])
        names.append(          names_ori[ i ])
        chrs.append(           chrs_ori[ i ])
        types.append(          types_ori[ i ])        
        

# -------------------------------------------------------------
# check gene prot record
if b_check_gene_prot_rec:
    print(prog_tag+ " check genes/prot/rna/stem_loop records --------------")
    for i in range(len(names)):
        print("name:"+names[i]+"\ttype:"+types[i]+" ["+str(starts[i])+", "+str(ends[i])+"]")
    print(prog_tag+ "---------------------------------------")
# -------------------------------------------------------------


print("creates "+json_annot_f+" file")
with open(json_annot_f, 'w+') as f:
    f.write("{\"genomesize\": "+ genome_length +", ")
    if len(chrs) == 1:
        f.write(f"\"virus_id\": \""+contig+"\", ")
    else:
        f.write("\"virus_id\": \"unkown\", ")
    b_first_cds_found = False

    # **********************************************
    # for genes
    f.write("\"genes\": {")
    b_gene_found = False
    for i in range(len(names)):
        if types[i] == 'gene':
            if b_first_cds_found:
                f.write(", ")

            f.write("\""+names[i]+"\": ["+str(starts[i])+", "+str(ends[i])+"]")
            b_first_cds_found = True
            b_gene_found = True

    # ------------------------------------------
    # special case of not recorded genes
    # means genes not recorded, like in PCV2,
    # therefore we must deduce them from proteins, calling them ORF1...ORFN
    if not b_gene_found:
        b_first_cds_found = False
        for i in range(len(names)):
            if types[i] == 'cds':
                tmp_name = 'gene_'+str(i+1)
                if b_first_cds_found:
                    f.write(", ")

                f.write("\""+tmp_name+"\": ["+str(starts[i])+", "+str(ends[i])+"]")
                b_first_cds_found = True
    # ------------------------------------------

    f.write("}")
    # **********************************************

    b_first_cds_found = False
    # for proteins
    f.write(", \"proteins\": {")
    for i in range(len(names)):
        if types[i] == 'cds':
            if b_first_cds_found:
                f.write(", ")

            f.write("\""+names[i]+"\": ["+str(starts[i])+", "+str(ends[i])+"]")
            b_first_cds_found = True
    f.write("}")

    b_first_cds_found = False
    # for RNA stem_loops
    f.write(", \"rnas\": {")
    for i in range(len(names)):
        if types[i] == 'rna':
            if b_first_cds_found:
                f.write(", ")

            f.write("\""+names[i]+"\": ["+str(starts[i])+", "+str(ends[i])+"]")
            b_first_cds_found = True
    f.write("}")

    b_first_cds_found = False
    # for RNA stem_loops
    f.write(", \"stem_loops\": {")
    for i in range(len(names)):
        if types[i] == 'stem_loop':
            if b_first_cds_found:
                f.write(", ")

            f.write("\""+names[i]+"\": ["+str(starts[i])+", "+str(ends[i])+"]")
            b_first_cds_found = True
    f.write("}")
    
    f.write("}")
f.close()



# # no used because keep original contig start position for vardict and real position used
# # are obtained from json
# print(f"creates {bed_annot_f} file")
# with open(bed_annot_f, 'w+') as f:
#     for i in range(len(names)):
#             f.write(f"{chrs[i]}\t{starts[i]}\t{ends[i]}\t{names[i]}\n")
# f.close()

print("creates "+bed_vardict_annot_f+" file")
with open(bed_vardict_annot_f, 'w+') as f:
    # no header otherwise bug in vardict
    for i in range(len(names)):
        # print("for vardict, treat BED contig "+names[i])
        f.write("\t".join([
            chrs[i],
            str(starts_vardict[i]),
            str(ends_vardict[i]),
            names[i]+"\n"
        ]))
f.close()
