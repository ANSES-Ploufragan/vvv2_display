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
# from vcf file of vardict and bed file with contig positions, creates a vcf file with
# contigs shifted position (to create a picture of variants and annotations)
###
import argparse, os, sys, csv, re, warnings
from os import path
import subprocess

# to be able to report line number in error messages
import inspect
frame = inspect.currentframe()

# debug
b_test_correct_multicontig_vardict_vcf = False # ok 2022 04 29
b_test = False

prog_tag = '[' + os.path.basename(__file__) + ']'

seq_stat_f    = '' # in, seq stat file from vadr to deduce shift foreach contig and position
vardict_vcf_f = '' # in, vcf file of vardict, but with bad position for contigs except the first
correct_vcf_f = '' # out same as previous but positions on contigs are shifted according to previous contigs
contig_limits_f = '' # store positions of limits between contigs, to be used in visualize_snp...R

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seq_stat_f", dest='seq_stat_f',
                    help="seq stat file of vadr annotator",
                    metavar="FILE")
parser.add_argument("-b", "--vardict_vcf_f", dest='vardict_vcf_f',
                    help="[input] vcf file from vardict",
                    metavar="FILE")
parser.add_argument("-c", "--correct_vcf_f", dest='correct_vcf_f',
                    help="[input] vcf file from vardict",
                    metavar="FILE")
parser.add_argument("-l", "--contig_limits_f", dest='contig_limits_f',
                    help="[input] text file with limit positions between contigs",
                    metavar="FILE")
parser.add_argument("-z", "--test_all", dest='b_test_all',
                    help="[Optional] run all tests",
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
    b_test_correct_multicontig_vardict_vcf = True
    b_test = True
else:
    b_test = b_test_correct_multicontig_vardict_vcf

if ((not b_test)and
    ((len(sys.argv) < 6) or (len(sys.argv) > 10))):
    parser.print_help()
    print(prog_tag + "[Error] we found "+str(len(sys.argv)) +
          " arguments, exit line "+str(frame.f_lineno))
    sys.exit(0)

# print('args:', args)
if args.seq_stat_f is not None:
    seq_stat_f = os.path.abspath(args.seq_stat_f)
elif(not b_test):
    sys.exit(prog_tag + "[Error] You must provide seq_stat_f")
if args.vardict_vcf_f is not None:
    vardict_vcf_f = os.path.abspath(args.vardict_vcf_f)
elif(not b_test):
    sys.exit(prog_tag + "[Error] You must provide vardict_vcf_f")
if args.correct_vcf_f is not None:
    correct_vcf_f = os.path.abspath(args.correct_vcf_f)
elif(not b_test):
    sys.exit(prog_tag + "[Error] You must provide vardict_vcf_f")
if args.contig_limits_f is not None:
    contig_limits_f = os.path.abspath(args.contig_limits_f)
elif(not b_test):
    sys.exit(prog_tag  + "[Error] You must provide contig_limits_f")

if args.b_verbose is not None:
    b_verbose = args.b_verbose


if b_test_correct_multicontig_vardict_vcf:
    test_dir = 'test_correct_multicontig_vardict_vcf/'
    test_name = [
        'res'
        # 'res2'
    ]
    for resn in test_name:
        vardict_vcf_f   = f"{test_dir}{resn}.vardict.vcf"
        correct_vcf_f   = f"{test_dir}{resn}.correct.vcf"
        contig_limits_f = f"{test_dir}{resn}.contig_limits.txt"        
        seq_stat_f      = f"{test_dir}{resn}.vadr.seqstat"
        cmd = ' '.join(['./correct_multicontig_vardict_vcf.py',
                        f"-b {vardict_vcf_f}",
                        f"-c {correct_vcf_f}",
                        f"-l {contig_limits_f}",                        
                        f"-s {seq_stat_f}"
                        ])
        if b_verbose:
            cmd = f"{cmd} -v"
        print(prog_tag + " cmd:"+cmd)
        print(prog_tag + " START")
        os.system(cmd)
        print(prog_tag + " END")
    sys.exit()

# ----------------------------------------------------------------------
# reads seqstat file to get info on genome_length and on contigs lengths
# ----------------------------------------------------------------------
contig_names   = []
contig_lengths = []
contig_shifts  = {}
with open(seq_stat_f) as ssf:
    for line in ssf:
        if re.match(r"Total # residues:    ", line):
            line_fields =  line.split(' ')
            # print("line_fields:"+','.join(line_fields))
            genome_length = line.split(' ')[6]
            genome_length = genome_length.rstrip()
            if b_verbose:
                print(prog_tag + " genome_length:"+str(genome_length))
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
                print(prog_tag + " contig:"+contig_name+" length:"+str(contig_length))

ssf.close()

cumulated_shift = 0
contig_start_pos = []
for i in range(len(contig_names)):
    if i == 0:
        contig_shifts[ contig_names[i] ] = 0
    else:
        cumulated_shift += int(contig_lengths[i-1])
        contig_start_pos.append(cumulated_shift)
        contig_shifts[ contig_names[i] ] = cumulated_shift
        if b_verbose:
            print(prog_tag + ' '+ "record shift "+str(contig_shifts[ contig_names[i] ])+" for contig "+ contig_names[i])
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# write contig_limit output file, used by visualise_snp...R to draw vertical lines
# related to contig limits on the graph
# ----------------------------------------------------------------------
print(prog_tag + " contig_start_pos:"+str(contig_start_pos))

l = open(contig_limits_f, 'w+')
for pos in contig_start_pos:
    l.write(str(pos)+"\n")
print(prog_tag + ' '+ contig_limits_f +" file created")
l.close()
# ----------------------------------------------------------------------

print(prog_tag + " Reads "+vardict_vcf_f+" file")
print(prog_tag + " Creates "+correct_vcf_f+" file")
o = open(correct_vcf_f, 'w+')

with open(vardict_vcf_f, 'r') as f:
    for line in f:
        if re.match(r'^#', line):            
            o.write(line)
            # sys.exit(f"header found:{line}")
        else:
            # sys.exit(f"fields found:{line}")
            fields = line.split()
            if b_verbose:
                print(f"fields:{fields[1]} becomes ")
            # correct POS
            fields[1] = str( contig_shifts[ fields[0] ] + int(fields[1])  )

            if b_verbose:            
                print(f"fields:{fields[1]}\n")
            o.write("\t".join(fields) + "\n")
f.close()
print(prog_tag + ' '+ correct_vcf_f + " file created")
o.close()
