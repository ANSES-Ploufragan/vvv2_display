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
b_test_cov_depth_corr_f = False #  2023 12 19
b_test = False

prog_tag = '[' + os.path.basename(__file__) + ']'

cov_depth_f    = '' # in, cov_depth file provided by samtools depth -aa file.bam
cov_depth_corr_f = '' # out same as previous but positions on contigs are shifted according to previous contigs

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--cov_depth_f", dest='cov_depth_f',
                    help="in: coverage depth file provided by samtools -aa in.bam (all contigs)",
                    metavar="FILE")
parser.add_argument("-c", "--cov_depth_corr_f", dest='cov_depth_corr_f',
                    help="out: coverage depth file corrected for position (cumulative position in case of several contigs)",
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
    b_test_cov_depth_corr_f = True
    b_test = True
else:
    b_test = b_test_cov_depth_corr_f

if ((not b_test)and
    ((len(sys.argv) < 2) or (len(sys.argv) > 10))):
    parser.print_help()
    print(prog_tag + "[Error] we found "+str(len(sys.argv)) +
          " arguments, exit line "+str(frame.f_lineno))
    sys.exit(0)

# print('args:', args)
if args.cov_depth_f is not None:
    cov_depth_f = os.path.abspath(args.cov_depth_f)
elif(not b_test):
    sys.exit(prog_tag + "[Error] You must provide cov_depth_f")
if args.cov_depth_corr_f is not None:
    cov_depth_corr_f = os.path.abspath(args.cov_depth_corr_f)
elif(not b_test):
    sys.exit(prog_tag + "[Error] You must provide cov_depth_corr_f")

if args.b_verbose is not None:
    b_verbose = args.b_verbose

if b_test_cov_depth_corr_f:
    test_dir = '../test_vvv2_display/'
    test_name = [
        'res'
        # 'res2'
    ]
    for resn in test_name:
        cov_depth_f         = test_dir + "res2_covdepth.txt"
        cov_depth_corr_f = test_dir + "res2_covdepth_corrected.txt"
        cmd = ' '.join(['./correct_covdepth_f.py',
                        "-s", cov_depth_f,
                        "-c", cov_depth_corr_f
                        ])
        if b_verbose:
            cmd = cmd + " -v"
        print(prog_tag + " cmd:"+cmd)
        print(prog_tag + " START")
        os.system(cmd)
        print(prog_tag + " END")
    sys.exit()

# ----------------------------------------------------------------------
# reads seqstat file to get info on genome_length and on contigs lengths
# ----------------------------------------------------------------------
print(prog_tag + " Reads "+cov_depth_f+" file")
print(prog_tag + " Creates "+cov_depth_corr_f+" file")
o = open(cov_depth_corr_f, 'w+')

# init var
prev_contig_nr = 1
contig_num   = 0

pos = 0
curr_shift = 0
with open(cov_depth_f) as cdf:
    for line in cdf:
        line_fields = line.split("\t")
        # print("line_fields:"+','.join(line_fields))
        contig_nr = line_fields[0]
        contig_pos = line_fields[1]
        cov_depth = line_fields[2]
        # print("line_fields:"+','.join(line_fields))
        if b_verbose:
            print(prog_tag + " contig:"+contig_name+" length:"+str(contig_length))

        # write in output file
        # if contig change, record shift to apply
        if contig_nr != prev_contig_nr:
            curr_shift = pos
        # get position WITH SHIFT
        pos = curr_shift + int(contig_pos)
        pos_str = str(pos)
        # write pos and related coverage depth
        o.write("\t".join([pos_str, cov_depth]))
        if b_verbose:            
            print("record pos:",pos_str,", cov_depth:", cov_depth)

        prev_contig_nr = contig_nr
cdf.close()
o.close()

print(prog_tag + ' '+ cov_depth_corr_f + " file created")
