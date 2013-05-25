#!/usr/bin/env python

import os
import sys

import setuptools

from distutils.core import setup
#from ez_setup import use_setuptools
#use_setuptools()
#from setuptools import setup

# convenience is king
opj = os.path.join

# make sure org_settings.cfg is in source directory
org_settings_fn = 'org_settings.cfg'
dist_settings_path = opj(os.getcwd(),'src','chipsequtil',org_settings_fn)
if not os.path.exists(dist_settings_path) :
    sys.stderr.write('WARNING: %s could not be found \
                      in distribution root directory.  org_settings.py script may \
                      not work properly.\n'%dist_settings_path)

scripts = ['scripts/align_pipeline.py',
           'scripts/build_chipseq_infosite.py',
           'scripts/chipseq_pipeline.py',
           'scripts/combine_gerald_stats.py',
           'scripts/compare_microarray_binding.py',
           'scripts/create_pipeline_script.py',
           'scripts/create_alignment_script.py',
           'scripts/extract_promoters.py',
           'scripts/filter_bed_by_position_count.py',
           'scripts/filter_macs_peaks.py',
           'scripts/filter_gps_peaks.py',
           'scripts/filter_mapped_known_genes.py',
           'scripts/generate_stats_doc.py',
           'scripts/gerald_stats.py',
           'scripts/gerald_to_bed.py',
           'scripts/integrate_macs_ucsc.py',
           'scripts/join_mapped_known_genes.py',
           'scripts/map_intervals.py',
           'scripts/map_peaks_to_genes.py',
           'scripts/map_peaks_to_known_genes.py',
           'scripts/motif_scan.py',
           'scripts/nibFrag.py',
           'scripts/org_settings.py',
           'scripts/peaks_to_fasta.py',
           'scripts/plot_pos_vs_neg_peaks.py',
           'scripts/plot_peak_loc_dist.py',
           'scripts/probeset_to_known_gene.py',
           'scripts/rejection_sample_fasta.py',
           'scripts/sort_bed.py',
           'scripts/split_file.py',
           'scripts/split_qsub.py',
           'scripts/THEME.sh',
           'scripts/wait_for_qsub.py',
           'scripts/wait_for_jobid.py',
           'scripts/wqsub.py',
           'scripts/wqsub_drmaa.py',
           ]

# setup and install
setup(name='chipsequtil',
      version='0.6',
      author='Adam Labadorf; Shao-shan Carol Huang',
      author_email='alabadorf@gmail.com; shhuang1@gmail.com',
      package_dir={'':'src'},
      py_modules=['chipsequtil.nib','chipsequtil.util','chipsequtil.plotting',
                  'chipsequtil.sampling','chipsequtil.seq','chipsequtil.motifs'],
      packages=['chipsequtil'],
      package_data={'': ['org_settings.cfg','tool_settings.cfg','static/css/*']},
      scripts=scripts,
      #cmdclass={'uninstall': uninstall},
     )
