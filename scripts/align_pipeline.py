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
parent_parser_bowtie.add_argument('--reads1-files',dest='reads1_files',action='append',required=True,help='comma-separated list of read files')
parent_parser_bowtie.add_argument('--reads2-files',dest='reads2_files',action='append',default=[],help='comma-separated list of read files')
parent_parser_bowtie.add_argument('--reads-format',dest='reads_format',choices=['.fastq.gz','.fastq'],default='.fastq.gz',help='Format of reads file')
parent_parser_bowtie.add_argument('--output-file',dest='output_file',action='append',default=[],help='output file names [default: current directory name]')
parent_parser_bowtie.add_argument('--output-format',dest='output_format',choices=['.sam','.bam'],default='.sam',help='Format of output alignment file [default: %(default)s')
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

parser_tophat = subparsers.add_parser('tophat',parents=[parent_parser_bowtie])
parser_tophat.add_argument('--tophat-opts',dest='tophat_opts',default='',help='double quote wrapped options for tophat [default: %(default)s]')
parser_tophat.add_argument('--bowtie-version',dest='bowtie_version',choices=['bowtie'],required=True,help='bowtie version')
parser_tophat.add_argument('--bowtie-ind-loc',dest='bowtie_ind_loc',required=True,help='bowtie index location')
parser_tophat.add_argument('--tophat-gene-mod',dest='tophat_gene_mod',help='gene model annotation file in GTF2.2 or GFF3 format')
parser_tophat.add_argument('--tophat-tx-ind-loc',dest='tophat_tx_ind_loc',help='transcriptome index location')
parser_tophat.add_argument('--tophat-aln-nthread',dest='tophat_aln_nthread',default=2,type=int,help='number of threads to run the mapping [default: %(default)s]')

parser_tophat2 = subparsers.add_parser('tophat2',parents=[parent_parser_bowtie])
parser_tophat2.add_argument('--tophat2-opts',dest='tophat2_opts',default='',help='double quote wrapped options for tophat2 [default: %(default)s]')
parser_tophat2.add_argument('--bowtie-version',dest='bowtie_version',choices=['bowtie','bowtie2'],required=True,help='bowtie version')
parser_tophat2.add_argument('--bowtie-ind-loc',dest='bowtie_ind_loc',required=True,help='bowtie index location')
parser_tophat2.add_argument('--tophat2-gene-mod',dest='tophat2_gene_mod',help='gene model annotation file in GTF2.2 or GFF3 format')
parser_tophat2.add_argument('--tophat2-tx-ind-loc',dest='tophat2_tx_ind_loc',help='transcriptome index location')
parser_tophat2.add_argument('--tophat2-aln-nthread',dest='tophat2_aln_nthread',default=2,type=int,help='number of threads to run the mapping [default: %(default)s]')

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
    for (reads1f,reads2f,outf,final_d) in zip(args.reads1_files,args.reads2_files,args.output_file,args.final_dir):
        outf_root = os.path.splitext(outf)[0]
        temp_file_fp = os.path.join(args.temp_dir,outf)     
        temp_file_fp_root = os.path.join(args.temp_dir,outf_root)
        final_file_fp = os.path.join(final_d,outf)
        final_file_fp_root = os.path.join(final_d,outf_root)
        if (aligner=='bowtie'):
            bowtie_d = {'bowtie_exec':'bowtie',
            'bowtie_opts':args.bowtie_opts,
            'bowtie_aln_nthread':args.bowtie_aln_nthread,
            'bowtie_ind_loc':args.bowtie_ind_loc,
            'bowtie_hit':'\"%s\"'%(temp_file_fp),
            'bowtie_log':'\"%s.log\"'%(final_file_fp)}
            if (args.reads_format=='.fastq.gz'):
                if (reads2f==''): # single-end only
#                    bowtie_d['reads_files'] = '\"<(zcat %s)\"'%(' '.join(['\"%s\"'%(f) for f in reads1f.split(',')]))
                    bowtie_d['reads_files'] = '<(gzip -dc %s)'%(' '.join(['%s'%(f) for f in reads1f.split(',')]))
                else: # paired end
                    bowtie_d['reads_files'] = '-1 <(zcat %s) -2 <(zcat %s)'%(' '.join(['\"%s\"'%(f) for f in reads1f.split(',')]),
                                                                             ' '.join(['\"%s\"'%(f) for f in reads2f.split(',')]))
                #calls.append("gzip -dc %(reads_files)s | %(bowtie_exec)s %(bowtie_opts)s -p %(bowtie_aln_nthread)s %(bowtie_ind_loc)s - %(bowtie_hit)s 2>&1 | tee %(bowtie_log)s"%(bowtie_d))
            else:
                if (reads2f==''): # single-end only
                    bowtie_d['reads_files'] = reads1f
                else: # paired end
                    bowtie_d['reads_files'] = '-1 %s -2 %s'%(reads1f,reads2f)
            calls.append("bash -c \"%(bowtie_exec)s %(bowtie_opts)s -p %(bowtie_aln_nthread)s %(bowtie_ind_loc)s %(reads_files)s %(bowtie_hit)s 2>&1 | tee %(bowtie_log)s\""%(bowtie_d))
        elif (aligner=='bowtie2'):
            bowtie2_d = {'bowtie2_exec':'bowtie2',
                         'bowtie2_opts':args.bowtie2_opts,
                         'bowtie2_aln_nthread':args.bowtie2_aln_nthread,
                         'bowtie2_ind_loc':args.bowtie2_ind_loc,
                         'bowtie2_hit':'\"%s\"'%temp_file_fp,
                         'bowtie2_log':'\"%s.log\"'%(final_file_fp)}
            if (reads2f==''):
                bowtie2_d['reads_files'] = '-U %s'%(reads1f)
            else:
                bowtie2_d['reads_files'] = '-1 %s -2 %s'%(reads1f,reads2f)
            calls.append("%(bowtie2_exec)s %(bowtie2_opts)s -p %(bowtie2_aln_nthread)s -x %(bowtie2_ind_loc)s %(reads_files)s -S %(bowtie2_hit)s 2>&1 | tee %(bowtie2_log)s"%(bowtie2_d))
        elif (aligner=='tophat'):
            tophat_d = {'tophat_exec':'tophat',
                        'tophat_opts':args.tophat_opts,
                        'tophat_ebwt_base':args.bowtie_ind_loc,
                        'tophat_aln_nthread':args.tophat_aln_nthread,
                        'tophat_output_dir':'\"%s\"'%(temp_file_fp_root),
                        'tophat_log':'\"%s.log\"'%(final_file_fp_root)}
            if (reads2f==''):
                tophat_d['reads_files'] = reads1f
            else:
                tophat_d['reads_files'] = '%s %s'%(reads1f,reads2f)
            if (args.tophat_gene_mod==None):
                tophat_d['tophat_gene_mod'] = ''
            else:
                tophat_d['tophat_gene_mod'] = '-G %s'%(args.tophat_gene_mod)
            if (args.tophat_tx_ind_loc==None):
                tophat_d['tophat_tx_ind_loc'] = ''
            else:
                tophat_d['tophat_tx_ind_loc'] = '--transcriptome-index=%s'%(args.tophat_tx_ind_loc)
            calls.append("bash -c \"%(tophat_exec)s %(tophat_opts)s %(tophat_gene_mod)s %(tophat_tx_ind_loc)s -o %(tophat_output_dir)s -p %(tophat_aln_nthread)s %(tophat_ebwt_base)s %(reads_files)s 2>&1 | tee %(tophat_log)s\""%(tophat_d))
        elif (aligner=='tophat2'):
            tophat2_d = {'tophat2_exec':'tophat2',
                        'tophat2_opts':args.tophat2_opts,
                        'tophat2_ebwt_base':args.bowtie_ind_loc,
                        'tophat2_aln_nthread':args.tophat2_aln_nthread,
                        'tophat2_output_dir':'\"%s\"'%(temp_file_fp_root),
                        'tophat2_log':'\"%s.log\"'%(final_file_fp_root)}
            if (reads2f==''):
                tophat2_d['reads_files'] = reads1f
            else:
                tophat2_d['reads_files'] = '%s %s'%(reads1f,reads2f)
            if (args.bowtie_version=='bowtie'):
                tophat2_d['bowtie_version_opt'] = '--bowtie1'
            else:
                tophat2_d['bowtie_version_opt'] = ''
            if (args.tophat2_gene_mod==None):
                tophat2_d['tophat2_gene_mod'] = ''
            else:
                tophat2_d['tophat2_gene_mod'] = '-G %s'%(args.tophat2_gene_mod)
            if (args.tophat2_tx_ind_loc==None):
                tophat2_d['tophat2_tx_ind_loc'] = ''
            else:
                tophat2_d['tophat2_tx_ind_loc'] = '--transcriptome-index=%s'%(args.tophat2_tx_ind_loc)
            calls.append("bash -c \"%(tophat2_exec)s %(tophat2_opts)s %(bowtie_version_opt)s %(tophat2_gene_mod)s %(tophat2_tx_ind_loc)s -o %(tophat2_output_dir)s -p %(tophat2_aln_nthread)s %(tophat2_ebwt_base)s %(reads_files)s 2>&1 | tee %(tophat2_log)s\""%(tophat2_d))
    
    steps.append(PPS('Run %s'%(aligner),calls,env=os.environ))

    ###################
    # sam2bam
    ##################
    calls = []
    for outf,final_d in zip(args.output_file,args.final_dir):
        outf_root = os.path.splitext(outf)[0]
        temp_file_fp = os.path.join(args.temp_dir,outf)
        if (aligner in ['bowtie','bowtie2']):
            temp_bam_fp = os.path.join(args.temp_dir,outf_root)
            final_bam_fp = os.path.join(final_d,outf_root)
            if (args.output_format=='.sam'):
                temp_sam_fp = temp_file_fp
                sam2bam_d = {'samtools_exec':'samtools',
                             'samtools_ind_loc':args.samtools_ind_loc,
                             'temp_sam_fp':temp_sam_fp,
                             'temp_bam_fp':temp_bam_fp,
                             'final_bam_fp':final_bam_fp,
                             }
                calls.append("%(samtools_exec)s view -uSt \"%(samtools_ind_loc)s\" \"%(temp_sam_fp)s\" | %(samtools_exec)s sort - \"%(temp_bam_fp)s\""%sam2bam_d)
                calls.append("cp \"%(temp_bam_fp)s.bam\" \"%(final_bam_fp)s.bam\""%(sam2bam_d))
            elif (args.output_format=='.bam'):
                sam2bam_d = {'samtools_exec':'samtools',
                             'samtools_ind_loc':args.samtools_ind_loc,
                             'temp_bam_fp':temp_bam_fp,
                             'final_bam_fp':final_bam_fp,
                             }
                calls.append("%(samtools_exec)s sort \"%(temp_bam_fp)s.bam \"%(temp_bam_fp)s.sorted\""%sam2bam_d)
                calls.append("cp \"%(temp_bam_fp)s.sorted.bam\" \"%(final_bam_fp)s.bam\""%(sam2bam_d))
        elif (aligner in ['tophat','tophat2']):
            temp_bam_dir = os.path.join(args.temp_dir,outf_root)
            temp_bam_fp = os.path.join(temp_bam_dir,'accepted_hits')
            final_bam_dir = os.path.join(final_d,outf_root)
            final_bam_fp = os.path.join(final_bam_dir,'accepted_hits')
            calls.append("mkdir -p \"%s\""%(final_bam_dir))
            if (args.output_format=='.sam'):
                temp_sam_fp = os.path.join(args.temp_dir,outf_root,'accepted_hits.sam')
                sam2bam_d = {'samtools_exec':'samtools',
                             'samtools_ind_loc':args.samtools_ind_loc,
                             'temp_sam_fp':temp_sam_fp,
                             'temp_bam_dir':temp_bam_dir,
                             'temp_bam_fp':temp_bam_fp,
                             'final_bam_dir':final_bam_dir,
                             'final_bam_fp':final_bam_fp,
                             }
                calls.append("%(samtools_exec)s view -uSt \"%(samtools_ind_loc)s\" \"%(temp_sam_fp)s\" | %(samtools_exec)s sort - \"%(temp_bam_fp)s\""%sam2bam_d)
                calls.append("cd \"%(temp_bam_dir)s\"; find . -type f -not -name accepted_hits.sam -print0 | cpio -0 -pdm \"%(final_bam_dir)s\""%(sam2bam_d))
            elif (args.output_format=='.bam'):
                sam2bam_d = {'samtools_exec':'samtools',
                             'samtools_ind_loc':args.samtools_ind_loc,
                             'temp_bam_dir':temp_bam_dir,
                             'temp_bam_fp':temp_bam_fp,
                             'final_bam_dir':final_bam_dir,
                             'final_bam_fp':final_bam_fp,
                             }
                calls.append("%(samtools_exec)s sort \"%(temp_bam_fp)s.bam\" \"%(temp_bam_fp)s.sorted\""%sam2bam_d)
                calls.append("cd \"%(temp_bam_dir)s\"; find . -type f -not -name accepted_hits.bam -print0 | cpio -0 -pdm \"%(final_bam_dir)s\""%(sam2bam_d))
    steps.append(PPS("Create BAM file and move to final destination",calls,env=os.environ))

    pipeline.add_steps(steps)
    if args.auto and args.steplist is not None :
        steplist = parse_steplist(args.steplist,pipeline)
    else :
        steplist = None
    pipeline.run(interactive=not args.auto,steplist=steplist)
