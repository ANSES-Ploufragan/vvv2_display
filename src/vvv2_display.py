#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of the vvv2_display distribution (https://github.com/ANSES-Ploufragan/vvv2_display).
# Copyright (c) 2024 Fabrice Touzain.
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
# vvv2_display script: from vardict (variant calling) and vadr (annotator) results,
# creates a picture of variants alongside the detected viral genome
###
import argparse, os, sys, warnings, re, subprocess
# hashlib: no, needs python 3.11 to have sha256 on files, trigger installation problème
# we use ssytem sha256sum function instead
from os import path
import subprocess

# to be able to report line number in error messages
import inspect
frame = inspect.currentframe()

##### MAIN
def __main__():

    # debug
    b_test_vvv2_display                    = False # ok 2022 05 05 complet, partial tc
    b_test_convert_tbl2json                = False # ok 2022 04 26 complete tc,
    b_test_correct_multicontig_vardict_vcf = False # ok 2022 04 29 partial tc
    b_test_convert_vcffile_to_readable     = False # ok 2025 01 20 complete tc,
    b_test_correct_covdepth_f              = False # ok 2024 01 20 tc,
    b_test_visualize_snp_v4                = False # ok 2024 03 27 complete tc,
    b_test = False
    dir_path = os.path.dirname(os.path.abspath(__file__)) # dir of current script

    b_verbose = False
    var_significant_threshold = 7 # default value
    var_significant_threshold_str = str(var_significant_threshold)
    # allow to run tests from everywhere
    
    prog_tag = '[' + os.path.basename(__file__) + ']'

    # to record if we display cov depth in graph or not (depends on provided intputs)
    b_cov_depth_display = False

    # to get argument telling if covdepth ordinates are in log10 scale (default yes)
    b_log_scale_int   = 1
    b_log_scale_str   = str(b_log_scale_int)
    b_log_scale       = True

    # --------------------------------------
    # input files
    # --------------------------------------
    pass_annot_f    = '' # in, pass tbl file from vadr
    fail_annot_f    = '' # in, fail tbl file from vadr
    seq_stat_f      = '' # in, seq stat file from vadr
    vardict_vcf_f   = '' # in, vcf file from VarDict
    correct_vcf_f   = '' # in, vcf file corrected for pos (when multicontigs)
    cov_depth_f     = '' # in, cov depth file provided by samtools depth -aa
    cov_depth_corr_f= '' # in, cov depth corrected for pos (when pulticontigs)
    # --------------------------------------

    # --------------------------------------
    # intermediate files
    # --------------------------------------
    json_annot_f         = '' # json annotation file: out
    # bed_annot_f          = '' # bed  annotation file: out
    bed_vardict_annot_f  = '' # bed  annotation file for vardict: out
    snp_loc_f            = '' # text file for snp location
    snp_loc_summary_f    = '' # text file for snp location summary
    contig_limits_f      = '' # in, txt file with contig limits positions    
    contig_names_f       = '' # in, txt file with contig names    
    # --------------------------------------

    # --------------------------------------
    # final file(s)
    # --------------------------------------
    png_var_f    = '' # png output with figure illustrating variant
                      # proportion and annotations
    # --------------------------------------


    #########################################
    # directories fo sub programs
    #########################################
    PYTHON_SCRIPTS = dir_path + "/" # PYTHON_SCRIPTS/"
    R_SCRIPTS      = dir_path + "/" # R_SCRIPTS/"
    #########################################

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pass_tbl_f", dest='pass_annot_f',
                        help="in: tabular file of vadr annotations, with pass status",
                        metavar="FILE")
    parser.add_argument("-f", "--fail_tbl_f", dest='fail_annot_f',
                        help="in: tabular file of vadr annotations, with fail status",
                        metavar="FILE")
    parser.add_argument("-s", "--seq_stat_f", dest='seq_stat_f',
                        help="in: seq stat file of vadr annotator",
                        metavar="FILE")
    parser.add_argument("-n", "--vcf_f", dest='vardict_vcf_f',
                        help="in: vcf variant file provided by vardict",
                        metavar="FILE")
    parser.add_argument("-r", "--png_var_f", dest='png_var_f',
                        help="out: png file with variant proportions and annotations",
                        metavar="FILE")
    parser.add_argument("-w", "--var_significant_threshold", dest='var_significant_threshold',
                        help="(percentage var_significant_threshold) Define minimal proportion of a variant to be kept in significant results",
                        type=int)
    parser.add_argument("-y", "--covdepth_linear_scale", dest='b_log_scale',
                        help="[Optional] to display covepth ordinates in linear scale (default log10 scale)",
                        action='store_false')
    parser.add_argument("-o", "--cov_depth_f", dest='cov_depth_f',
                        help="[optional] in: text file of coverage depths (given by samtools depth)",
                        metavar="FILE")     
    parser.add_argument("-e", "--cov_depth_corr_f", dest='cov_depth_corr_f',
                        help="[optional] out: text file of coverage depths with cumulated position in case of several contigs, for display (tmp file, for galaxy compatibility)",
                        metavar="FILE")                       
    parser.add_argument("-t", "--snp_loc_f", dest='snp_loc_f',
                        help="[optional] out: variant description for relevant positions, txt file (if not provided, file name deduced from png name)",
                        metavar="FILE")  
    parser.add_argument("-u", "--snp_loc_summary_f", dest='snp_loc_summary_f',
                        help="[optional] out: variant description for relevant positions, txt file (if not provided, file name deduced from png name)",
                        metavar="FILE")  
    parser.add_argument("-j", "--json_f", dest='json_annot_f',
                        help="[Optional] out (tmp out file for galaxy compatibility, no need in other cases): vadr annotation converted to json",
                        metavar="FILE")
    parser.add_argument("-k", "--bed_f", dest='bed_vardict_annot_f',
                        help="[Optional] out (tmp out file for galaxy compatibility, no need in other cases): vardict variants adapted to annotation, in bed format",
                        metavar="FILE")
    parser.add_argument("-l", "--cvcf_f", dest='correct_vcf_f',
                        help="[Optional] out (tmp out file for galaxy compatibility, no need in other cases): vardict variants corrected for positions when several contigs, in vcf format",
                        metavar="FILE")
    parser.add_argument("-m", "--contig_limits_f", dest='contig_limits_f',
                        help="[Optional] out (tmp out file for galaxy compatibility, no need in other cases): txt file of contig limits",
                        metavar="FILE")
    parser.add_argument("-N", "--contig_names_f", dest='contig_names_f',
                        help="[Optional] out (tmp out file for galaxy compatibility, no need in other cases): txt file of contig names",
                        metavar="FILE")
    parser.add_argument("-z", "--test_vvv2_display", dest='b_test_vvv2_display',
                        help="[Optional] run all tests",
                        action='store_true')
    parser.add_argument("-a", "--test_convert_tbl2json", dest='b_test_convert_tbl2json',
                        help="[Optional] run test on tbl2json conversion",
                        action='store_true')
    parser.add_argument("-b", "--test_correct_multicontig_vardict_vcf", dest='b_test_correct_multicontig_vardict_vcf',
                        help="[Optional] run test on vcf correction of vardict vcf file",
                        action='store_true')
    parser.add_argument("-c", "--test_convert_vcffile_to_readable", dest='b_test_convert_vcffile_to_readable',
                        help="[Optional] run test on vcf to readable txt file conversion",
                        action='store_true')
    parser.add_argument("-d", "--test_visualize_snp_v4", dest='b_test_visualize_snp_v4',
                        help="[Optional] run test to visualize snp in a png",
                        action='store_true')
    parser.add_argument("-g", "--test_correct_covdepth_f", dest='b_test_correct_covdepth_f',
                        help="[Optional] run test to correct position in cov depth file",
                        action='store_true')
    parser.add_argument("-v", "--verbose", dest='b_verbose',
                        help="[Optional] To have details on records when running",
                        action='store_true')
    # parser.set_defaults(b_test_all=False)
    # parser.set_defaults(b_test_vvv2_display=False)
    # parser.set_defaults(b_test_convert_tbl2json=False)
    # parser.set_defaults(b_test_convert_vcffile_to_readable=False)
    # parser.set_defaults(b_test_visualize_snp_v4=False)
    parser.set_defaults(b_verbose=False)
    parser.set_defaults(var_significant_threshold=7)
    parser.set_defaults(b_log_scale=True)
    if var_significant_threshold is not None:
        var_significant_threshold_str = str(var_significant_threshold)

    # get absolute path in case of files
    args = parser.parse_args()
    print(prog_tag + " arguments obtained, checking... line "+str(frame.f_lineno))

    # -------------------------------------------
    # check arguments
    b_test_vvv2_display                    = args.b_test_vvv2_display
    b_test_convert_tbl2json                = args.b_test_convert_tbl2json
    b_test_correct_multicontig_vardict_vcf = args.b_test_correct_multicontig_vardict_vcf
    b_test_convert_vcffile_to_readable     = args.b_test_convert_vcffile_to_readable
    b_test_visualize_snp_v4                = args.b_test_visualize_snp_v4
    b_test_correct_covdepth_f              = args.b_test_correct_covdepth_f

    if b_test_vvv2_display:
        b_test_convert_tbl2json                = True
        b_test_correct_multicontig_vardict_vcf = True
        b_test_convert_vcffile_to_readable     = True
        b_test_visualize_snp_v4                = True
        b_test_correct_covdepth_f              = True
        b_test                                 = True
    else:
        b_test = (b_test_vvv2_display                    or
                  b_test_convert_tbl2json                or
                  b_test_correct_multicontig_vardict_vcf or
                  b_test_convert_vcffile_to_readable     or
                  b_test_correct_covdepth_f              or
                  b_test_visualize_snp_v4)
        # print(f"b_test:{b_test}")
        # print(f"b_test_convert_tbl2json:{b_test_convert_tbl2json}")    

    if ((not b_test)and
        ((len(sys.argv) < 9) or (len(sys.argv) > 30))):
        print("\n".join([prog_tag,
                         "Aim: Display of SNP proportions, annotations, for an assembly",
                         "in:", 
                         " - vardict variant calling output",
                         " - vadr assembly annotations",
                         "out:",
                         " - png file (image of SNP proportion alongside the assembly with CDS positions)",
                         " - txt file with variant calling summary, location in CDS and surround DNA sequence.\n"]))
        parser.print_help()
        print(prog_tag + "[Error] we found "+str(len(sys.argv)) +
              " arguments, exit line "+str(frame.f_lineno))
        sys.exit(0)

    if not b_test_vvv2_display:
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
        if args.vardict_vcf_f is not None:
            vardict_vcf_f = os.path.abspath(args.vardict_vcf_f)
        elif(not b_test):
            sys.exit("[Error] You must provide vcf_f")
        if args.png_var_f is not None:
            png_var_f = os.path.abspath(args.png_var_f)
        elif(not b_test):
            sys.exit("[Error] You must provide png_var_f name for output")
        if args.snp_loc_f is not None:
            snp_loc_f = os.path.abspath(args.snp_loc_f)
        if args.snp_loc_summary_f is not None:
            snp_loc_summary_f = os.path.abspath(args.snp_loc_summary_f)
                   
        if( (args.cov_depth_f is not None) and (os.path.isfile(args.cov_depth_f)) ):
            cov_depth_f = os.path.abspath(args.cov_depth_f)
            b_cov_depth_display = True
        else:
            b_cov_depth_display = False
        if args.b_log_scale:
            b_log_scale_int   = 1
        else:
            b_log_scale_int   = 0
        b_log_scale_str   = str(b_log_scale_int)

        # ----------------------------------------------------------------
        # optional arguments only for Galaxy compatibility
        if args.json_annot_f is not None:
            json_annot_f = os.path.abspath(args.json_annot_f)
        if args.bed_vardict_annot_f is not None:
            bed_vardict_annot_f = os.path.abspath(args.bed_vardict_annot_f)
        if args.correct_vcf_f is not None:
            correct_vcf_f = os.path.abspath(args.correct_vcf_f)


        if args.contig_limits_f is None:
            # if name of contig_limits file is not provided by user, deduce a name file with a part
            # deduced from pass_annot_f using sha256 checksum
            # to avoid bad interferences of serveral vvv2_display runs results  
            
            cmd = "/usr/bin/sha256sum "+pass_annot_f
            # print(prog_tag + " cmd:"+cmd)
            sha256_f = subprocess.getoutput(cmd).split()[0]
            contig_limits_f = str(sha256_f) + '_contig_limits.txt'
            # print(prog_tag + " contig_limits file name created:"+ contig_limits_f)
        else:
            contig_limits_f = os.path.abspath(args.contig_limits_f)

        if args.contig_names_f is None:
            # if name of contig_names file is not provided by user, deduce a name file with a part
            # deduced from pass_annot_f using sha256 checksum
            # to avoid bad interferences of serveral vvv2_display runs results  
            
            cmd = "/usr/bin/sha256sum "+pass_annot_f
            # print(prog_tag + " cmd:"+cmd)
            sha256_f = subprocess.getoutput(cmd).split()[0]
            contig_names_f = str(sha256_f) + '_contig_names.txt'
            # print(prog_tag + " contig_names file name created:"+ contig_names_f)
        else:
            contig_names_f = os.path.abspath(args.contig_names_f)
        
        if args.cov_depth_corr_f is not None:
            cov_depth_corr_f = os.path.abspath(args.cov_depth_corr_f)
        # ----------------------------------------------------------------
        
        if args.var_significant_threshold is not None:
            var_significant_threshold = args.var_significant_threshold
            var_significant_threshold_str = str(var_significant_threshold)
        if args.b_verbose is not None:
            b_verbose = args.b_verbose
        print(prog_tag + " arguments checked... line "+str(frame.f_lineno))


    # ------------------------------------------------------------------
    # TEST for vvv2_display
    # ------------------------------------------------------------------
    test_dir = dir_path + "/../test_vvv2_display"
    if b_test_vvv2_display:
        # --------------------------------------------------------------
        # COMPLETE GENOME
        # in files
        pass_annot_f  = test_dir + "/res2_vadr_pass.tbl" # from vadr results
        fail_annot_f  = test_dir + "/res2_vadr_fail.tbl" # from vadr results
        seq_stat_f    = test_dir + "/res2_vadr.seqstat"  # from vadr results
        vardict_vcf_f = test_dir + "/res2_vardict.vcf"   # from lofreq results  
        cov_depth_f   = test_dir + "/res2_covdepth.txt"  # from samtools results  
        # tmp out files
        json_annot_f  = test_dir + "/res2_vadr.json"     # from convert_tbl2json.py
        contig_limits_f= test_dir + "/contig_limits.txt"
        contig_names_f = test_dir + "/contig_names.txt"
        cov_depth_corr_f= test_dir + "/res2_covdepth_corr.txt"
        # final out file
        png_var_f     = test_dir + "/res2_vvv2.png"     # from ...
        cmd = ' '.join([ dir_path + "/vvv2_display.py",
                            "--pass_tbl_f", pass_annot_f,
                            "--fail_tbl_f", fail_annot_f,
                            "--seq_stat_f", seq_stat_f,
                            "--vcf_f", vardict_vcf_f,
                            "--contig_limits_f", contig_limits_f,   
                            "--contig_names_f", contig_names_f,   
                            "--cov_depth_f", cov_depth_f,  
                            "--cov_depth_corr_f", cov_depth_corr_f,              
                            "--png_var_f", png_var_f,
                            "--var_significant_threshold", var_significant_threshold_str
                    ])
        print(prog_tag + " START")    
        print(prog_tag + " cmd:" + cmd)
        os.system(cmd)
        print(prog_tag + " END")
        # --------------------------------------------------------------

        # --------------------------------------------------------------
        # CONTIGS
        # in files
        pass_annot_f  = test_dir + "/res_vadr_pass.tbl" # from vadr results
        fail_annot_f  = test_dir + "/res_vadr_fail.tbl" # from vadr results
        seq_stat_f    = test_dir + "/res_vadr.seqstat"  # from vadr results
        vardict_vcf_f = test_dir + "/res_vardict.vcf"  # from lofreq results    
        cov_depth_f   = test_dir + "/res_covdepth.txt"  # from samtools results
        # tmp out files
        json_annot_f    = test_dir + "/res_vadr.json"     # from convert_tbl2json.
        contig_limits_f = test_dir + "/contig_limits.txt" 
        contig_names_f = test_dir + "/contig_names.txt"  
        cov_depth_corr_f= test_dir + "/res_covdepth_corr.txt"
        # final out file
        png_var_f    = test_dir + "/res_vvv2.png"     # from ...
        cmd = ' '.join([ dir_path + "/vvv2_display.py",
                    "--pass_tbl_f", pass_annot_f,
                    "--fail_tbl_f", fail_annot_f,
                    "--seq_stat_f", seq_stat_f,
                    "--vcf_f", vardict_vcf_f,
                    "--contig_limits_f", contig_limits_f,   
                    "--contig_names_f", contig_names_f,   
                    "--cov_depth_f", cov_depth_f,                     
                    "--cov_depth_corr_f", cov_depth_corr_f,       
                    "--png_var_f", png_var_f,
                    "--var_significant_threshold", var_significant_threshold_str                    
                    ])
        print(prog_tag + " START")    
        print(prog_tag + " cmd:" + cmd)
        os.system(cmd)
        print(prog_tag + " END")
        # --------------------------------------------------------------

        # --------------------------------------------------------------
        # COMPLETE GENOME PCV2
        # in files
        pass_annot_f  = test_dir + "/res3_vadr_pass.tbl" # from vadr results
        fail_annot_f  = test_dir + "/res3_vadr_fail.tbl" # from vadr results
        seq_stat_f    = test_dir + "/res3_vadr.seqstat"  # from vadr results
        vardict_vcf_f = test_dir + "/res3_vardict.vcf"   # from lofreq results  
        cov_depth_f   = test_dir + "/res3_covdepth.txt"  # from samtools results  
        # tmp out files
        json_annot_f  = test_dir + "/res3_vadr.json"     # from convert_tbl2json.py
        contig_limits_f= test_dir + "/contig3_limits.txt"
        contig_names_f = test_dir + "/contig3_names.txt"
        cov_depth_corr_f= test_dir + "/res3_covdepth_corr.txt"
        # final out file
        png_var_f     = test_dir + "/res3_vvv2.png"     # from ...
        cmd = ' '.join([ dir_path + "/vvv2_display.py",
                            "--pass_tbl_f", pass_annot_f,
                            "--fail_tbl_f", fail_annot_f,
                            "--seq_stat_f", seq_stat_f,
                            "--vcf_f", vardict_vcf_f,
                            "--contig_limits_f", contig_limits_f,   
                            "--contig_names_f", contig_names_f,   
                            "--cov_depth_f", cov_depth_f,  
                            "--cov_depth_corr_f", cov_depth_corr_f,              
                            "--png_var_f", png_var_f,
                            "--var_significant_threshold", var_significant_threshold_str
                    ])
        print(prog_tag + " START")    
        print(prog_tag + " cmd:" + cmd)
        os.system(cmd)
        print(prog_tag + " END")
        # --------------------------------------------------------------

        sys.exit()
    # ------------------------------------------------------------------    


    # ------------------------------------------------------------------
    # convert tbl fileS of vadr annotator to json annotation file
    # ------------------------------------------------------------------
    if b_test_convert_tbl2json:
        # # COMPLETE GENOME
        # pass_annot_f = test_dir + "/res2_vadr_pass.tbl" # from vadr results
        # fail_annot_f = test_dir + "/res2_vadr_fail.tbl" # from vadr results
        # seq_stat_f   = test_dir + "/res2_vadr.seqstat"  # from vadr results
        # json_annot_f = test_dir + "/res2_vadr.json"
        ## bed_annot_f  = test_dir + "/res2_vadr.bed"
        # bed_vardict_annot_f  = test_dir + "/res2_vadr.4vardict.bed"        
        # CONTIGS
        pass_annot_f         = test_dir + "/res_vadr_pass.tbl" # from vadr results
        fail_annot_f         = test_dir + "/res_vadr_fail.tbl" # from vadr results
        seq_stat_f           = test_dir + "/res_vadr.seqstat"  # from vadr results
        json_annot_f         = test_dir + "/res_vadr.json"
        # bed_annot_f          = test_dir + "/res_vadr.bed"
        bed_vardict_annot_f  = test_dir + "/res_vadr.4vardict.bed"        

    if(json_annot_f == ''):
       json_annot_f = vardict_vcf_f
#       json_annot_f = json_annot_f.replace('.vcf', '.json')
       json_annot_f = re.sub(r'\.[^\.]+$', '.json', json_annot_f)       
    if(bed_vardict_annot_f == ''):
       bed_vardict_annot_f = vardict_vcf_f
#       bed_vardict_annot_f = bed_vardict_annot_f.replace('.vcf', '.vardict.bed')
       bed_vardict_annot_f = re.sub(r'\.[^\.]+$', '_vardict.bed', bed_vardict_annot_f)
       
    p_script = PYTHON_SCRIPTS + "convert_tbl2json.py"
    cmd = ' '.join([ p_script,
                    "--pass_annot_f", pass_annot_f,
                    "--fail_annot_f", fail_annot_f,
                    "--seq_stat_f", seq_stat_f,
                    "--json_out_f", json_annot_f,
    #               "--bed_out_f", bed_annot_f,
                    "--bed_vardict_out_f", bed_vardict_annot_f                             
                    ])
    print(prog_tag + " cmd:" + cmd)

    if b_test_convert_tbl2json:
        print(prog_tag + " [test_convert_tbl2json] START")
        print(prog_tag + " cmd:" + cmd)
        os.system(cmd)    
        print(prog_tag + " [test_convert_tbl2json] END")
        sys.exit()
    else:
        os.system(cmd)    

    # ------------------------------------------------------------------
    # correct vardict vcf file for multicontig assemblies to get the good positions (cumulative)
    # in the final created picture
    # ------------------------------------------------------------------
    if b_test_correct_multicontig_vardict_vcf:
        seq_stat_f           = test_dir + "/res_vadr.seqstat"  # from vadr results    
        vardict_vcf_f        = test_dir + "/res_vardict.vcf"  # from vardict results            
        correct_vcf_f        = test_dir + "/res_correct.vcf"  # corrected out results
        contig_limits_f      = test_dir + "/contig_limits.txt"  # contig limits
        contig_names_f       = test_dir + "/contig_names.txt"  # contig names

    if(correct_vcf_f == ''):
       correct_vcf_f = vardict_vcf_f
#       correct_vcf_f = correct_vcf_f.replace('vardict.vcf', '_correct.vcf')
       correct_vcf_f = re.sub(r'\.[^\.]+$', '_correct.vcf', correct_vcf_f)
         
    p_script = PYTHON_SCRIPTS + "correct_multicontig_vardict_vcf.py"
    print("p_script:" + p_script)
    cmd = ' '.join([p_script,
                    "--seq_stat_f", seq_stat_f,
                    "--vardict_vcf_f", vardict_vcf_f,
                    "--correct_vcf_f", correct_vcf_f,
                    "--contig_limits_f", contig_limits_f,
                    "--contig_names_f", contig_names_f
                    ])
    print(prog_tag + " cmd:" + cmd)

    if b_test_correct_multicontig_vardict_vcf:
        print(prog_tag + " [test_correct_multicontig_vardict_vcf] START")
        print(prog_tag + " cmd:" + cmd)
        os.system(cmd)    
        print(prog_tag + " [test_correct_multicontig_vardict_vcf] END")
        sys.exit()
    else:
        print(prog_tag + " cmd:" + cmd)
        os.system(cmd)

    # ------------------------------------------------------------------
    # convert variant vcf file to text file usable by R as a dataframe
    # ------------------------------------------------------------------
    if b_test_convert_vcffile_to_readable:
        # # COMPLETE GENOME
        # vardict_vcf_f = test_dir + "/res2_vardict.vcf"  # from vardict results
        # correct_vcf_f = test_dir + "/res2_correct.vcf"  # corrected vardict results                
        # json_annot_f      = test_dir + "/res2_vadr.json"
        # snp_loc_f         =  test_dir + "/res2_snp.txt"
        # snp_loc_summary_f =  test_dir + "/res2_snp_summary.txt"
        # CONTIGS
        # vardict_vcf_f = test_dir + "/res_vardict.vcf"  # from vardict results
        correct_vcf_f     = test_dir + "/res_correct.vcf"                  
        json_annot_f      = test_dir + "/res_vadr.json"
        snp_loc_f         = test_dir + "/res_snp.txt"
        snp_loc_summary_f = test_dir + "/res_snp_summary.txt"

    if(snp_loc_f == ''):
       snp_loc_f = vardict_vcf_f
       snp_loc_f = re.sub(r'\.[^\.]+$', '_snp.txt', snp_loc_f)
    if(snp_loc_summary_f == ''):   
       snp_loc_summary_f = vardict_vcf_f
       snp_loc_summary_f = re.sub(r'\.[^\.]+$', '_snp_summary.txt', snp_loc_summary_f)
      
    # vcf file from vardict
    # json annotation file deduced from vadr, later vigor4(5?)
    # json_annot_f = f"{ech}_gene_position_viral_consensus.json" 
    # snp_loc_f = f"{ech}_snp_location"
    # p_script = f"{PYTHON_SCRIPTS}convert_vcffile_to_readablefile.py" # use pyvcf
    p_script = PYTHON_SCRIPTS + "convert_vcffile_to_readablefile2.py" # use pysam
    threshold = "%.2f" % (var_significant_threshold / 100)
    cmd = ' '.join([p_script,
                    "--vcfs", correct_vcf_f,
                    "--json", json_annot_f,
                    "--out", snp_loc_f,
                    "--outs", snp_loc_summary_f,
                    "--threshold", threshold])
    print(prog_tag + " cmd:" + cmd)

    if b_test_convert_vcffile_to_readable:
        print(prog_tag + " [test_convert_vcffile_to_readable] START")
        print(prog_tag + " cmd:" + cmd)
        os.system(cmd)    
        print(prog_tag + " [test_convert_vcffile_to_readable] END")
        sys.exit()
    else:
        os.system(cmd)    
    # ------------------------------------------------------------------

    # creates corrected cov depth file
    if b_cov_depth_display:
        print("b_cov_depth_display:"+str(b_cov_depth_display)+ ", line " + str(frame.f_lineno))
        if b_test_correct_covdepth_f and (not b_test_vvv2_display):
            cov_depth_f = test_dir + "/res_vvv2_covdepth.txt"
            cov_depth_corr_f = test_dir + "/res_vvv2_covdepth_corrected.txt"
        elif cov_depth_corr_f == '':
            # create a file name if not given, deducing it from cov_depth_f
            cov_depth_corr_f = cov_depth_f
            cov_depth_corr_f = re.sub(r'\.[^\.]+$', '_corrected.txt', cov_depth_corr_f)
            # handle case with cov_depth_f without extension (file provided by user)
            if cov_depth_corr_f == cov_depth_f:
                cov_depth_corr_f = cov_depth_corr_f + '_corrected.txt'
        cmd = " ".join([
                    PYTHON_SCRIPTS + "correct_covdepth_f.py",
                    "--cov_depth_f", cov_depth_f,
                    "--cov_depth_corr_f", cov_depth_corr_f
                ]) 
        if b_test_correct_covdepth_f and (not b_test_vvv2_display):
            print(prog_tag + " [test_correct_covdepth_f] START")
            print(prog_tag + " cmd:" + cmd)
            os.system(cmd)    
            print(prog_tag + " [test_correct_covdepth_f] END")
            sys.exit()
        else:
            print(prog_tag + " cmd:" + cmd)
            os.system(cmd)

    # ------------------------------------------------------------------
    # creates png graphic of variants from snp file and threshold
    # ------------------------------------------------------------------
    if b_test_visualize_snp_v4:
        # # COMPLETE GENOME
        # snp_loc_f =  test_dir + "/res2_snp.txt"
        # snp_loc_summary_f =  test_dir + "/res2_snp_summary.txt"
        # png_var_f =  test_dir + "/res2_snp.png"
        # CONTIGS
        snp_loc_f         = test_dir + "/res_snp.txt"
        snp_loc_summary_f = test_dir + "/res_snp_summary.txt"
        json_annot_f      = test_dir + "/res_vadr.json"
        png_var_f         = test_dir + "/res_snp.png"
        contig_limits_f   = test_dir + "/contig_limits.txt"
        contig_names_f    = test_dir + "/contig_names.txt"
        b_log_scale_int   = 1
        b_log_scale       = str(b_log_scale_int)
        
    if(png_var_f == ''):
       png_var_f = snp_loc_f
       png_var_f = re.sub(r'\.[^\.]+$', '.png', png_var_f)
    
    # env r-env.yaml
    r_script = R_SCRIPTS + "visualize_snp_v4.R"
    threshold = "%.2f" % (var_significant_threshold / 100) #"0.07"
    
    # png_var_f = f"{ech}_graphic_variant.png"
    cmd = " ".join(["R --vanilla --quiet --args",
                snp_loc_f, 
                contig_limits_f,
                contig_names_f,
                threshold,
                json_annot_f,
                png_var_f,
                b_log_scale_str])

    # parameter to allow coverage depth display above variants/annotations
    if b_cov_depth_display:
        cmd = cmd + " " + cov_depth_corr_f + " "
    
    if b_verbose:
        cmd = cmd + " ".join([" < ", r_script, " > /dev/null"])
    else:
        cmd = cmd + " ".join([" < ", r_script])
    print(prog_tag + " cmd:" + cmd)

    if b_test_visualize_snp_v4:
        print(prog_tag + " [test_visualize_snp_v4] START")
        print(prog_tag + " cmd:" + cmd)
        os.system(cmd)
        print(prog_tag + " [test_visualize_snp_v4] END")
    else:
        os.system(cmd)    
    # ------------------------------------------------------------------
    print(png_var_f + " file created")
    
    # remove useless file
    if os.path.isfile(contig_limits_f):
        os.unlink(contig_limits_f)
    if os.path.isfile(contig_names_f):
        os.unlink(contig_names_f)

##### MAIN END
if __name__=="__main__":__main__()
