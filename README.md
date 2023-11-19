# vvv2_display

# Description

Tools to create:
- a .png image file describing all variants (obtained from vardict-java variant caller) alongside a genome/assembly (to provide) with their proportion (ordinates), with CDS descriptions (obtained from vadr annotator).

Python/R scripts and Galaxy wrapper to use them.

It uses the results of:
- vadr >= 1.4.1 for annotation (of reference/assembly)
- vardict-java 1.8.3 for variant calling (of BAM alignement using reference/assembly and reads)

# Programs

- ```vvv2_display.py```: main script running each step of analyses
This script can be run independently, once __vvv2__ conda environment is installed and activated.
Type ```./vvv2_display.py``` then enter to get help on how to use it.

- ```PYTHON_SCRIPTS/convert_tbl2json.py```: 
Convert ```vadr``` annotation output .tbl file to json

- ```PYTHON_SCRIPTS/convert_vcffile_to_readablefile.py```: 
Convert ```vardict-java``` variant calling vcf file to human readable txt file

- ```PYTHON_SCRIPTS/correct_multicontig_vardict_vcf.py```: 
Correct ```vadr``` annotation output .tbl file for contigs positions when the assembly provided is composed of more than one contig.


<!-- - ```R_SCRIPTS/visualize_coverage_depth.R```: -->
<!-- Create a .png file showing coverage depth alongside the genome, from a bam alignment file. -->

- ```R_SCRIPTS/visualize_snp_v4.R```:
Create a .png file showing variant proportions alongside the genome/assembly and CDS positions.

# Installation

Use conda environment:
```
conda create -n vvv2_display -y
conda activate vvv2_display
mamba/conda install -c bioconda vvv2_display
```
Prefer mamba installation if completely new conda environments (faster). Do not mix mamba and conda.

Description:
```
vvv2_display.py -h
```

Typical usage:
```
vvv2_display.py -p res_vadr_pass.tsv -f res_vadr_fail.tsv -s res_vadr_seqstat.txt -n res_vardict_all.vcf -m contig_limits.txt -r res_vvv2_display.png 
```

# Galaxy wrapper

- ```vvv2_display.xml```:
Allow Galaxy integration of ```vvv2_display.py```. vvv2_display can be used in Galaxy pipelines.

# Citation

Please, if you use __vvv2_display__ and publish results, cite:
- Lai, Zhongwu, Aleksandra Markovets, Miika Ahdesmaki, Brad Chapman, Oliver Hofmann, Robert McEwen, Justin Johnson, Brian Dougherty, J. Carl Barrett, and Jonathan R. Dry. “__VarDict__: A Novel and Versatile Variant Caller for next-Generation Sequencing in Cancer Research.” Nucleic Acids Research 44, no. 11 (June 20, 2016): e108–e108. https://doi.org/10.1093/nar/gkw227.
- Schäffer, Alejandro A., Eneida L. Hatcher, Linda Yankie, Lara Shonkwiler, J. Rodney Brister, Ilene Karsch-Mizrachi, and Eric P. Nawrocki. “__VADR__: Validation and Annotation of Virus Sequence Submissions to GenBank.” BMC Bioinformatics 21, no. 1 (December 2020): 211. https://doi.org/10.1186/s12859-020-3537-3.
- Flageul, Alexandre, Pierrick Lucas, Edouard Hirchaud, Fabrice Touzain, Yannick Blanchard, Nicolas Eterradossi, Paul Brown, and Béatrice Grasland. “__Viral Variant Visualizer (VVV)__: A Novel Bioinformatic Tool for Rapid and Simple Visualization of Viral Genetic Diversity.” Virus Research 291 (January 2021): 198201. https://doi.org/10.1016/j.virusres.2020.198201.

