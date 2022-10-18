# -*- coding: utf-8 -*-
###
# USE PYTHON3
# vvv2_display script: from vardict (variant calling) and vadr (annotator) results,
# creates a picture of variants alongside the detected viral genome
###
import argparse, os, sys, warnings
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
    b_test_convert_vcffile_to_readable     = False # ok 2022 04 28 complete tc,
    b_test_visualize_snp_v4                = False # ok 2022 04 28 complete tc,
    b_test = False
    dir_path = os.path.dirname(os.path.abspath(__file__)) # dir of current script
    # allow to run tests from everywhere
    
    prog_tag = '[' + os.path.basename(__file__) + ']'

    # --------------------------------------
    # input files
    # --------------------------------------
    pass_annot_f  = '' # in, pass tbl file from vadr
    fail_annot_f  = '' # in, fail tbl file from vadr
    seq_stat_f    = '' # in, seq stat file from vadr
    vardict_vcf_f = '' # in, vcf file from VarDict
    correct_vcf_f = '' # in, vcf file corrected for pos (when multicontigs)
    # --------------------------------------

    # --------------------------------------
    # intermediate files
    # --------------------------------------
    json_annot_f         = '' # json annotation file: out
    # bed_annot_f          = '' # bed  annotation file: out
    bed_vardict_annot_f  = '' # bed  annotation file for vardict: out
    snp_loc_f            = '' # text file for snp location
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
    PYTHON_SCRIPTS = f"{dir_path}/" # PYTHON_SCRIPTS/"
    R_SCRIPTS      = f"{dir_path}/" # R_SCRIPTS/"
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
    parser.add_argument("-v", "--verbose", dest='b_verbose',
                        help="[Optional] To have details on records when running",
                        action='store_true')
    # parser.set_defaults(b_test_all=False)
    # parser.set_defaults(b_test_vvv2_display=False)
    # parser.set_defaults(b_test_convert_tbl2json=False)
    # parser.set_defaults(b_test_convert_vcffile_to_readable=False)
    # parser.set_defaults(b_test_visualize_snp_v4=False)
    parser.set_defaults(b_verbose=False)

    # get absolute path in case of files
    args = parser.parse_args()

    # -------------------------------------------
    # check arguments
    b_test_vvv2_display                    = args.b_test_vvv2_display
    b_test_convert_tbl2json                = args.b_test_convert_tbl2json
    b_test_correct_multicontig_vardict_vcf = args.b_test_correct_multicontig_vardict_vcf
    b_test_convert_vcffile_to_readable     = args.b_test_convert_vcffile_to_readable
    b_test_visualize_snp_v4                = args.b_test_visualize_snp_v4

    if b_test_vvv2_display:
        b_test_convert_tbl2json                = True
        b_test_correct_multicontig_vardict_vcf = True
        b_test_convert_vcffile_to_readable     = True
        b_test_visualize_snp_v4                = True
        b_test                                 = True
    else:
        b_test = (b_test_vvv2_display                            or
                  b_test_convert_tbl2json                or
                  b_test_correct_multicontig_vardict_vcf or
                  b_test_convert_vcffile_to_readable     or
                  b_test_visualize_snp_v4)
        # print(f"b_test:{b_test}")
        # print(f"b_test_convert_tbl2json:{b_test_convert_tbl2json}")    

    if ((not b_test)and
        ((len(sys.argv) < 9) or (len(sys.argv) > 17))):
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

    if args.b_verbose is not None:
        b_verbose = args.b_verbose


    # ------------------------------------------------------------------
    # TEST for vvv2_display
    # ------------------------------------------------------------------
    test_dir = f"{dir_path}/../test_vvv2_display"
    if b_test_vvv2_display:
        # --------------------------------------------------------------
        # COMPLETE GENOME
        # in files
        pass_annot_f  = f"{test_dir}/res2_vadr_pass.tbl" # from vadr results
        fail_annot_f  = f"{test_dir}/res2_vadr_fail.tbl" # from vadr results
        seq_stat_f    = f"{test_dir}/res2_vadr.seqstat"  # from vadr results
        vardict_vcf_f = f"{test_dir}/res2_vardict.vcf"  # from lofreq results    
        # tmp out files
        json_annot_f  = f"{test_dir}/res2_vadr.json"     # from convert_tbl2json.py
        # final out file
        png_var_f     = f"{test_dir}/res2_vvv2.png"     # from ...
        cmd = ' '.join([f"python3 {dir_path}/vvv2_display.py",
                    f"--pass_tbl_f {pass_annot_f}",
                    f"--fail_tbl_f {fail_annot_f}",
                    f"--seq_stat_f {seq_stat_f}",
                    f"--vcf_f {vardict_vcf_f}",
                    f"--png_var_f {png_var_f}"
                    ])
        print(f"{prog_tag} START")    
        print(f"cmd:{cmd}")
        os.system(cmd)
        print(f"{prog_tag} END")
        # --------------------------------------------------------------

        # --------------------------------------------------------------
        # CONTIGS
        # in files
        pass_annot_f  = f"{test_dir}/res_vadr_pass.tbl" # from vadr results
        fail_annot_f  = f"{test_dir}/res_vadr_fail.tbl" # from vadr results
        seq_stat_f    = f"{test_dir}/res_vadr.seqstat"  # from vadr results
        vardict_vcf_f = f"{test_dir}/res_vardict.vcf"  # from lofreq results    
        # tmp out files
        json_annot_f  = f"{test_dir}/res_vadr.json"     # from convert_tbl2json.py
        # final out file
        png_var_f    = f"{test_dir}/res_vvv2.png"     # from ...
        cmd = ' '.join([f"python3 {dir_path}/vvv2_display.py",
                    f"--pass_tbl_f {pass_annot_f}",
                    f"--fail_tbl_f {fail_annot_f}",
                    f"--seq_stat_f {seq_stat_f}",
                    f"--vcf_f {vardict_vcf_f}",
                    f"--png_var_f {png_var_f}"
                    ])
        print(f"{prog_tag} START")    
        print(f"cmd:{cmd}")
        os.system(cmd)
        print(f"{prog_tag} END")
        # --------------------------------------------------------------

        sys.exit()
    # ------------------------------------------------------------------    


    # ------------------------------------------------------------------
    # convert tbl fileS of vadr annotator to json annotation file
    # ------------------------------------------------------------------
    if b_test_convert_tbl2json:
        # # COMPLETE GENOME
        # pass_annot_f = f"{test_dir}/res2_vadr_pass.tbl" # from vadr results
        # fail_annot_f = f"{test_dir}/res2_vadr_fail.tbl" # from vadr results
        # seq_stat_f   = f"{test_dir}/res2_vadr.seqstat"  # from vadr results
        # json_annot_f = f"{test_dir}/res2_vadr.json"
        ## bed_annot_f  = f"{test_dir}/res2_vadr.bed"
        # bed_vardict_annot_f  = f"{test_dir}/res2_vadr.4vardict.bed"        
        # CONTIGS
        pass_annot_f         = f"{test_dir}/res_vadr_pass.tbl" # from vadr results
        fail_annot_f         = f"{test_dir}/res_vadr_fail.tbl" # from vadr results
        seq_stat_f           = f"{test_dir}/res_vadr.seqstat"  # from vadr results
        json_annot_f         = f"{test_dir}/res_vadr.json"
        # bed_annot_f          = f"{test_dir}/res_vadr.bed"
        bed_vardict_annot_f  = f"{test_dir}/res_vadr.4vardict.bed"        

    if( (json_annot_f == '') or
       (not os.path.exists(json_annot_f)) ):
       json_annot_f = vardict_vcf_f
       json_annot_f = json_annot_f.replace('.vcf', '.json')
    if( (bed_vardict_annot_f == '') or
       (not os.path.exists(bed_vardict_annot_f)) ):
       bed_vardict_annot_f = vardict_vcf_f
       bed_vardict_annot_f = bed_vardict_annot_f.replace('.vcf', '.vardict.bed')

    p_script = f"{PYTHON_SCRIPTS}convert_tbl2json.py"
    cmd = ' '.join([f"python3 {p_script}",
                    f"--pass_annot_f {pass_annot_f}",
                    f"--fail_annot_f {fail_annot_f}",
                    f"--seq_stat_f {seq_stat_f}",
                    f"--json_out_f {json_annot_f}",
    #                f"--bed_out_f {bed_annot_f}",
                    f"--bed_vardict_out_f {bed_vardict_annot_f}"                             
                    ])
    print(f"cmd:{cmd}")

    if b_test_convert_tbl2json:
        print(f"{prog_tag} [test_convert_tbl2json] START")
        print(f"cmd:{cmd}")
        os.system(cmd)    
        print(f"{prog_tag} [test_convert_tbl2json] END")
        sys.exit()
    else:
        os.system(cmd)    

    # ------------------------------------------------------------------
    # correct vardict vcf file for multicontig assemblies to get the good positions (cumulative)
    # in the final created picture
    # ------------------------------------------------------------------
    if b_test_correct_multicontig_vardict_vcf:
        seq_stat_f           = f"{test_dir}/res_vadr.seqstat"  # from vadr results    
        vardict_vcf_f        = f"{test_dir}/res_vardict.vcf"  # from vardict results            
        correct_vcf_f        = f"{test_dir}/res_correct.vcf"  # corrected out results


    if( (correct_vcf_f == '') or
       (not os.path.exists(correct_vcf_f)) ):
       correct_vcf_f = vardict_vcf_f
       correct_vcf_f = correct_vcf_f.replace('vardict.vcf', 'correct.vcf')

    p_script = f"{PYTHON_SCRIPTS}correct_multicontig_vardict_vcf.py"
    print(f"p_script:{p_script}")
    cmd = ' '.join([f"python3 {p_script}",
                    f"--seq_stat_f {seq_stat_f}",
                    f"--vardict_vcf_f {vardict_vcf_f}",
                    f"--correct_vcf_f {correct_vcf_f}"
                    ])
    print(f"cmd:{cmd}")

    if b_test_correct_multicontig_vardict_vcf:
        print(f"{prog_tag} [test_correct_multicontig_vardict_vcf] START")
        print(f"cmd:{cmd}")
        os.system(cmd)    
        print(f"{prog_tag} [test_correct_multicontig_vardict_vcf] END")
        sys.exit()
    else:
        os.system(cmd)

    # ------------------------------------------------------------------
    # convert variant vcf file to text file usable by R as a dataframe
    # ------------------------------------------------------------------
    if b_test_convert_vcffile_to_readable:
        # # COMPLETE GENOME
        # vardict_vcf_f = f"{test_dir}/res2_vardict.vcf"  # from vardict results
        # correct_vcf_f = f"{test_dir}/res2_correct.vcf"  # corrected vardict results                
        # json_annot_f  = f"{test_dir}/res2_vadr.json"
        # snp_loc_f     =  f"{test_dir}/res2_snp.txt"
        # CONTIGS
        # vardict_vcf_f = f"{test_dir}/res_vardict.vcf"  # from vardict results
        correct_vcf_f = f"{test_dir}/res_correct.vcf"                  
        json_annot_f  = f"{test_dir}/res_vadr.json"
        snp_loc_f     =  f"{test_dir}/res_snp.txt"

    if( (snp_loc_f == '') or
       (not os.path.exists(snp_loc_f)) ):
       snp_loc_f = vardict_vcf_f
       snp_loc_f = snp_loc_f.replace('vardict.vcf', 'snp.txt')

    # vcf file from vardict
    # json annotation file deduced from vadr, later vigor4(5?)
    # json_annot_f = f"{ech}_gene_position_viral_consensus.json" 
    # snp_loc_f = f"{ech}_snp_location"
    # p_script = f"{PYTHON_SCRIPTS}convert_vcffile_to_readablefile.py" # use pyvcf
    p_script = f"{PYTHON_SCRIPTS}convert_vcffile_to_readablefile2.py" # use pysam
    threshold = "0.07"
    cmd = ' '.join([f"python3 {p_script}",
                    f"--vcfs {correct_vcf_f}",
                    f"--json {json_annot_f}",
                    f"--out {snp_loc_f}",
                    f"--threshold {threshold}"])
    print(f"cmd:{cmd}")

    if b_test_convert_vcffile_to_readable:
        print(f"{prog_tag} [test_convert_vcffile_to_readable] START")
        print(f"cmd:{cmd}")
        os.system(cmd)    
        print(f"{prog_tag} [test_convert_vcffile_to_readable] END")
        sys.exit()
    else:
        os.system(cmd)    
    # ------------------------------------------------------------------


    # ------------------------------------------------------------------
    # creates png graphic of variants from snp file and threshold
    # ------------------------------------------------------------------
    if b_test_visualize_snp_v4:
        # # COMPLETE GENOME
        # snp_loc_f =  f"{test_dir}/res2_snp.txt"
        # png_var_f =  f"{test_dir}/res2_snp.png"
        # CONTIGS
        snp_loc_f =  f"{test_dir}/res_snp.txt"
        png_var_f =  f"{test_dir}/res_snp.png"

    if( (png_var_f == '') or
       (not os.path.exists(png_var_f)) ):
       png_var_f = snp_loc_f
       png_var_f = png_var_f.replace('.txt', '.png')

    # env r-env.yaml
    r_script = f"{R_SCRIPTS}visualize_snp_v4.R"
    threshold = "0.07"
    # png_var_f = f"{ech}_graphic_variant.png"
    cmd = f"R --vanilla --args {snp_loc_f} {threshold} {png_var_f} < {r_script}"
    print(f"cmd:{cmd}")

    if b_test_visualize_snp_v4:
        print(f"{prog_tag} [test_visualize_snp_v4] START")
        print(f"cmd:{cmd}")
        os.system(cmd)
        print(f"{prog_tag} [test_visualize_snp_v4] END")
    else:
        os.system(cmd)    
    # ------------------------------------------------------------------
    print(f"{png_var_f} file created")

##### MAIN END
if __name__=="__main__":__main__()
