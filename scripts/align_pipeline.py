#!/usr/bin/env python

import os
from subprocess import Popen, PIPE
import string
import sys
import argparse

from pypeline import Pypeline, ProcessPypeStep as PPS, PythonPypeStep as PyPS, parse_steplist

usage = "%prog [options] aligner species"
description = """1st generation ChIP-Seq alignment pipeline:
"""

epilog = """Note: it is advised to leave the --*-args arguments unchanged
unless you really know what you're doing."""

parser = argparse.ArgumentParser(description=description,epilog=epilog)
subparsers = parser.add_subparsers(dest='aligner_name')

parent_parser_bowtie = argparse.ArgumentParser(add_help=False)
parent_parser_bowtie.add_argument('--auto',dest='auto',action='store_true',help='run all steps non-interactively (for batch mode, e.g.)')
parent_parser_bowtie.add_argument('--steplist',dest='steplist',default='',help='with --auto, run specific steps')
parent_parser_bowtie.add_argument('--run-name',dest='run_name',default=os.path.basename(os.getcwd()),help='name for the alignment run, used for convenience [default: current directory name]')
parent_parser_bowtie.add_argument('--reads-files',dest='reads_files',action='append',required=True,help='comma-separated list of read files')
parent_parser_bowtie.add_argument('--reads-format',dest='reads_format',choices=['.fastq.gz','.fastq'],default='.fastq.gz',help='Format of reads file')
parent_parser_bowtie.add_argument('--output-file',dest='output_file',action='append',default=[],help='output file names [default: current directory name]')
parent_parser_bowtie.add_argument('--output-format',dest='output_format',choices=['.sam'],default='.sam',help='Format of output alignment file [default: %(default)s')
parent_parser_bowtie.add_argument('--temp-dir',dest='temp_dir',default=os.getcwd(),help='directory to store the alignment output files [default: current directory name]')
parent_parser_bowtie.add_argument('--final-dir',dest='final_dir',action='append',default=[],help='final destination directory for output files')
parent_parser_bowtie.add_argument('--samtools-ind-loc',dest='samtools_ind_loc',required=True,help='path of samtools index')

parent_parser_bowtie.add_argument('--print-args',dest='print_args',action='store_true',help=argparse.SUPPRESS) # secret ninja option

parser_bowtie = subparsers.add_parser('bowtie',parents=[parent_parser_bowtie])
parser_bowtie.add_argument('--bowtie-opts',dest='bowtie_opts',default='',help='double quote wrapped options for bowtie [default: %(default)s]')
parser_bowtie.add_argument('--bowtie-ind-loc',dest='bowtie_ind_loc',required=True,help='bowtie index location')
parser_bowtie.add_argument('--bowtie-aln-nthread',dest='bowtie_aln_nthread',default=2,type=int,help='number of threads to run bowtie [default: %(default)s]')

parser_bowtie2 = subparsers.add_parser('bowtie2',parents=[parent_parser_bowtie])
parser_bowtie2.add_argument('--bowtie2-opts',dest='bowtie2_opts',default='',help='double quote wrapped options for bowtie2 [default: %(default)s]')
parser_bowtie2.add_argument('--bowtie2-ind-loc',dest='bowtie2_ind_loc',required=True,help='bowtie2 index location [default: %(default)s]')
parser_bowtie2.add_argument('--bowtie2-aln-nthread',dest='bowtie2_aln_nthread',default=2,type=int,help='number of threads to run bowtie2 [default: %(default)s]')

if __name__ == '__main__' :

    # parse command line arguments
    args = parser.parse_args(sys.argv[1:])

    if args.print_args:
        opts_strs = []
        all_opts = []
        all_opts.extend(parser._actions)
        for opt in all_opts:
            if isinstance(opt,argparse._SubParsersAction):
                opts_strs.append('    %s'%(str(getattr(args,opt.dest))))
                all_sub_opts = opt.choices[getattr(args,opt.dest)]._actions
                for sub_opt in all_sub_opts:
                    sub_opt_str = sub_opt.option_strings[0]
                    if sub_opt_str in ['-h','--help','--print-args']:
                        pass
                    elif sub_opt_str == '--steplist' and not args.auto :
                        pass
                    elif isinstance(sub_opt,argparse._StoreAction):
                        sub_arg = str(getattr(args,sub_opt.dest))
                        if sub_arg.count(' ') > 0 or sub_arg.find(' -') != -1 or sub_arg.find('(')!=-1 or sub_arg.find(')')!=-1 or sub_arg.startswith('-') or sub_arg.find('--') != -1 :
                            opts_strs.append('    %s="%s"'%(sub_opt_str,str(getattr(args,sub_opt.dest))))
                        else :
                            opts_strs.append('    %s=%s'%(sub_opt_str,str(getattr(args,sub_opt.dest))))
                    elif isinstance(sub_opt,argparse._StoreTrueAction) and getattr(args,sub_opt.dest) :
                        opts_strs.append('    %s'%sub_opt_str)
                    elif isinstance(sub_opt,argparse._AppendAction):
                        sub_arg = getattr(args,sub_opt.dest)
                        for a in sub_arg:
                            if a.count(' ') > 0 or a.find(' -') != -1 or a.find('(')!=-1 or a.find(')')!=-1 or a.startswith('-') or a.find('--') != -1 :
                                opts_strs.append('    %s="%s"'%(sub_opt_str,a))
                            else :
                                opts_strs.append('    %s=%s'%(sub_opt_str,a))
        opts_strs.append('    $@')
        sys.stdout.write(' \\\n'.join(opts_strs)+'\n')
        sys.exit(0)
    # end print-args

    # the pipeline
    #log_fn = os.path.join(opts.exp_name+'_pipeline.log')
    pipeline = Pypeline('Analysis pipeline for %s'%args.run_name)

    steps = []

    ############################################################################
    # run aligner
    ############################################################################
    aligner = args.aligner_name
    calls = []
    for (readsf,outf) in zip(args.reads_files,args.output_file):
        temp_file_fp = os.path.join(args.temp_dir,outf)
        if (aligner=='bowtie'):
            bowtie_d = {'bowtie_exec':'bowtie',
            'bowtie_opts':args.bowtie_opts,
            'bowtie_aln_nthread':args.bowtie_aln_nthread,
            'bowtie_ind_loc':args.bowtie_ind_loc,
            'bowtie_hit':'\"%s\"'%(temp_file_fp),
            'bowtie_log':'\"%s.log\"'%(temp_file_fp)}
            if (args.reads_format=='.fastq.gz'):
                bowtie_d['reads_files'] = ' '.join(['\"%s\"'%(f) for f in readsf.split(',')])
                calls.append("gzip -dc %(reads_files)s | %(bowtie_exec)s %(bowtie_opts)s -p %(bowtie_aln_nthread)s %(bowtie_ind_loc)s - %(bowtie_hit)s 2>&1 | tee %(bowtie_log)s"%(bowtie_d))
            else:
                bowtie_d['reads_files'] = readsf
                calls.append("%(bowtie_exec)s %(bowtie_opts)s -p %(bowtie_aln_nthread)s %(bowtie_ind_loc)s %(reads_files)s %(bowtie_hit)s 2>&1 | tee %(bowtie_log)s"%(bowtie_d))
        elif (aligner=='bowtie2'):
            bowtie2_d = {'bowtie2_exec':'bowtie2',
                         'reads_files':' '.join(['\"%s\"'%(f) for f in readsf.split(',')]),
                         'bowtie2_opts':args.bowtie2_opts,
                         'bowtie2_aln_nthread':args.bowtie2_aln_nthread,
                         'bowtie2_ind_loc':args.bowtie2_ind_loc,
                         'bowtie2_hit':'\"%s\"'%temp_file_fp,
                         'bowtie2_log':'\"%s.log\"'%(temp_file_fp)}
            calls.append("%(bowtie2_exec)s %(bowtie2_opts)s -p %(bowtie2_aln_nthread)s -x %(bowtie2_ind_loc)s -U %(reads_files)s -S %(bowtie2_hit)s 2>&1 | tee %(bowtie2_log)s"%(bowtie2_d))

    steps.append(PPS('Run %s'%(aligner),calls,env=os.environ))

    ###################
    # sam2bam
    ##################
    calls = []
    for outf,final_d in zip(args.output_file,args.final_dir):
        if (args.output_format=='.sam'):
            temp_sam_fp = os.path.join(args.temp_dir,outf)
            aln_bam = os.path.splitext(outf)[0]
            temp_bam_fp = os.path.join(args.temp_dir,aln_bam)
            final_bam_fp = os.path.join(final_d,aln_bam)
            sam2bam_d = {'samtools_exec':'samtools',
                         'samtools_ind_loc':args.samtools_ind_loc,
                         'temp_sam_fp':temp_sam_fp,
                         'temp_bam_fp':temp_bam_fp,
                         'final_bam_fp':final_bam_fp,
                         }
            calls.append("%(samtools_exec)s view -uSt \"%(samtools_ind_loc)s\" \"%(temp_sam_fp)s\" | samtools sort - \"%(temp_bam_fp)s\""%sam2bam_d)
            calls.append("mv \"%(temp_bam_fp)s.bam\" \"%(final_bam_fp)s.bam\""%(sam2bam_d))
    steps.append(PPS("Convert SAM to BAM and move to final destination",calls,env=os.environ))

    pipeline.add_steps(steps)
    if args.auto and args.steplist is not None :
        steplist = parse_steplist(args.steplist,pipeline)
    else :
        steplist = None
    pipeline.run(interactive=not args.auto,steplist=steplist)
