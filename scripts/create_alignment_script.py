#!/usr/bin/env python

from __future__ import with_statement
import csv
import getpass
import json
import os
import textwrap

try:
    import readline
    import glob
    readline.parse_and_bind("tab: complete")
    readline.set_completer_delims('')

    comp_states = {}
    def basic_complete_file(text,state) :
        #if text.strip() == '' :
        #    text = './'
        options = dict([(i,p) for i,p in enumerate(glob.glob(text+'*'))])
        return options.get(state,None)

    readline.set_completer(basic_complete_file)

except ImportError:
    print "Module readline not available."

import re
import stat
import sys
from optparse import OptionParser
from subprocess import Popen, PIPE

import chipsequtil
from chipsequtil import get_global_settings, get_local_settings, check_org_settings, GLOBAL_SETTINGS_FN, LOCAL_SETTINGS_FN, get_global_tool_settings, get_local_tool_settings, check_tool_settings, GLOBAL_TOOL_SETTINGS_FN, LOCAL_TOOL_SETTINGS_FN
import generate_uniform_symlinks

from terminalcontroller import TERM_ESCAPE, announce, warn, error, white, bold

usage = "%prog"
description = """Script for creating a custom run script for
aligning ChIPSeq/DNAse hypersensitivity reads.  User is asked for
paths and settings required for alignment the *align_pipeline.py*
utility and produces an executable run script with helpful information on how to
run it.  Also creates a JSON formatted file containing all the parameters for
this pipeline run."""
epilog = "Note: this script only works in Unix-style environments"
parser = OptionParser(usage=usage,description=description,epilog=epilog)


script_template = """\
#!/bin/bash

# required parameters for the pipeline
ALIGNER=%(aligner)s
SPECIES=%(species)s

# align_pipeline.py is the main workhorse of this analysis
# you may change any of the arguments below from their defaults

align_pipeline.py $ALIGNER $SPECIES \\
%(def_args)s
"""

start_text = """\
This is an interactive script that creates an executable script to use for
aligning ChIP-Seq/DNase-Seq reads. When prompted for sample sheet file name, tab
completion is available a la bash or tcsh shells. Press Ctrl-C at any time to
quit.
"""

end_text = """The script %(script_fn)s has been created to run this pipeline. \
The script can now be run with:

$> ./%(script_fn)s

Have a nice day."""



def wb(st) :
    sys.stdout.write(white(bold(st)))


def input(st,default=None) :

    if default is None :
        default_str = ''
    else :
        default_str = ' [default: ' + default + ' ] '

    out = None
    while out is None :
        out = raw_input(white(bold(st))+default_str+white(bold(':'))+' \n')
        if len(out) == 0 :
            out = default

    return out


if __name__ == '__main__' :

    TERM_ESCAPE = True

    try :

        pipeline_args = []  # list of (arg_key,arg_val) tuples

        # herro
        announce('ChIP-Seq Experiment Read Alignment Pipeline Script Generator')
        print textwrap.fill(start_text)

        opts, args = parser.parse_args(sys.argv[1:])
        if len(args) > 0 :
            warn("Arguments were passed, but this script doesn't accept any arguments, rudely ignoring them...\n")

        # this dictionary will be used to generate a JSON formatted file with
        # all the relevant settings for the pipeline
        json_dict = {}
        
        ############################################################################
        # get settings
        ############################################################################
        global_tool_settings = get_global_tool_settings()
        local_tool_settings = get_local_tool_settings()
        #valid_tool_settings = global_tool_settings.keys() + local_tool_settings.keys()
        all_tool_settings = {}
        all_tool_settings.update(global_tool_settings)
        all_tool_settings.update(local_tool_settings)

        ############################################################################
        # alignment run name
        ############################################################################
        def_run_name = os.path.basename(os.getcwd())
        run_name = input('Alignment run name',def_run_name)
        run_name = run_name.replace(' ','_')
        json_dict['run name'] = run_name

        ############################################################################
        # name of the sample sheet
        ############################################################################
        samp_fn_text = """The sample sheet specifies the location of the raw reads
(non-aligned) files and other experiment parameters.  See %(example_samp_sheet)s for 
an example."""

        def_samp_fn = 'aln_sample_sheet'
        samp_fn = input('Sample sheet file name',def_samp_fn)
        samp_fn = samp_fn.replace(' ','_')

        json_dict['sample sheet file name'] = samp_fn

        ############################################################################
        # name of the alignment output path
        ############################################################################
        aln_path_text = """The alignment output files will be stored temporarily in a 
temporary directory (default to current directory).  You might want to use a local temporary 
directory (not NFS) for faster processing.  The very last step of the alignment script will convert the the aligned 
files in this directory to BAM format and move to a common location on the NFS."""

        print textwrap.fill(aln_path_text)

        def_path = os.path.basename(os.getcwd())
        aln_temp_path = input('Temporary directory',def_path)
        aln_temp_path = aln_temp_path.replace(' ','_')
        aln_final_path = input('Final alignment file destination',def_path)
        aln_final_path = aln_final_path.replace(' ','_')

        json_dict['alignment temp directory'] = aln_temp_path
        json_dict['alignment final directory'] = aln_final_path

        ############################################################################
        # Alignment settings
        ############################################################################
        announce('Alignment tool configuration')
        valid_aln_tools = global_tool_settings.get('tools',dict()).get('alignment','').split() + local_tool_settings.get('tools',dict()).get('alignment','').split()
        valid_aln_tools.sort()

        aln_text = """\
Below are the alignment tools available on this system.  The pipeline will
use the settings for the entire execution. If you do not see a set of settings that 
correspond to files you need you may add your own to %(local_tool)s.  See %(glob_tool)s for details.
"""
        print textwrap.fill(aln_text%{'local_tool':LOCAL_TOOL_SETTINGS_FN,'glob_tool':GLOBAL_TOOL_SETTINGS_FN},break_long_words=False)
        print

        wb('Available settings\n')
        aln_sets = [(k,all_tool_settings[k]) for k in valid_aln_tools]
        for aln, settings in aln_sets :
            wb(aln.ljust(8))
            print ':', settings.get('description','No description')
            #for k,v in settings.items() :
            #    print ' '*4+k+": "+str(v)
        aln = ''

        while aln not in valid_aln_tools :
            aln = input('Choose alignment tool, one of ('+','.join(valid_aln_tools)+')')
        print

        json_dict['aligner'] = aln

        ###### bowtie 1 & 2 #####
        if (aln in ['bowtie','bowtie2']):
            all_aln_settings = all_tool_settings[aln]

            announce('%s options configuration'%aln)
            valid_aln_opt_settings = all_aln_settings.get('opt_sets','').split()
            aln_opt_text = """
Below are the %(aln)s alignment options available on this system.  The pipeline will
use the settings for the entire execution. If you do not see a set of settings that 
correspond to files you need you may add your own to %(local_tool)s.  See %(glob_tool)s for details.
"""
            print textwrap.fill(aln_opt_text%{'aln':aln,'local_tool':LOCAL_TOOL_SETTINGS_FN,'glob_tool':GLOBAL_TOOL_SETTINGS_FN},break_long_words=False)
            print

            wb('Available settings\n')
            aln_opt_sets = [(k,all_aln_settings.get(k,'')) for k in all_aln_settings.get('opt_sets','').split()]
            for aln_opt, settings in aln_opt_sets :
                wb(aln_opt.ljust(8))
                print ':', settings
            aln_opt = ''
            
            while aln_opt not in valid_aln_opt_settings :
                aln_opt = input('Choose '+aln+' options, one of ('+','.join(valid_aln_opt_settings)+')')
                print
            json_dict['%s options'%(aln)] = aln_opt
            pipeline_args.append(('--%s-opts'%(aln),all_tool_settings[aln][aln_opt]))

            announce('%s index configuration'%aln)
            valid_aln_ind_list = all_aln_settings.get('indices','').split()
            aln_ind_text = """\
Below are the %(aln)s alignment indices available on this system.  The pipeline will
use the settings for the entire execution. If you do not see a set of settings that 
correspond to files you need you may add your own to %(local_tool)s.  See %(glob_tool)s for details.
"""
            print textwrap.fill(aln_opt_text%{'aln':aln,'local_tool':LOCAL_TOOL_SETTINGS_FN,'glob_tool':GLOBAL_TOOL_SETTINGS_FN},break_long_words=False)
            print

            wb('Available settings\n')
            for aln_ind in valid_aln_ind_list:
                wb(aln_ind.ljust(8))
                print ':', all_aln_settings.get(aln_ind+'.description','No description')
            aln_ind = ''

            while aln_ind not in valid_aln_ind_list:
                aln_ind = input('Choose '+aln+' index, one of ('+','.join(valid_aln_ind_list)+')')
                print

            json_dict['%s index'%(aln)] = aln_ind
            species = all_tool_settings[aln][aln_ind+'.species']
            json_dict['species'] = species
            pipeline_args.append(('--%s-ind-loc'%(aln),all_tool_settings[aln][aln_ind+'.loc']))

            nthread = ''
            while not re.search('^\d+$',nthread) :
                nthread = input('How many threads would you like to use for mapping EACH experiment?','2')
                json_dict['nthread'] = nthread
                pipeline_args.append(('--%s-aln-nthread'%(aln),nthread))

        # end if (aln in ['bowtie','bowtie2'])

        ############################################################################
        # done with input, creating script and other stuff
        ############################################################################
        # parse the sample sheet and create symlinks for the read files
        wb('Reading sample sheet and create symlink for read files if needed')
        json_dict['destination directory'] = dict()
        json_dict['reads files'] = dict()
        dest_dir_pattern = all_tool_settings['NameRule']['dest_dir']
        reads_symlink_pattern = all_tool_settings['NameRule']['reads_symlink']
        samp_ifh = open(samp_fn,'r')
        samp_reader = csv.DictReader(samp_ifh,dialect='excel-tab')
        for samp_line in samp_reader:
            local_id = samp_line['Local-ID']
            dest_dir = os.path.join(aln_final_path,dest_dir_pattern%(samp_line))
            json_dict['destination directory'][local_id] = dest_dir
            if (os.path.exists(dest_dir)):
                if (not os.path.isdir(dest_dir)):
                    raise ValueError('Alignment file destination directory name %s exists but is not a directory'%(dest_dir))
            else:
                os.makedirs(dest_dir)
            pipeline_args.append(('--dest-dir',dest_dir))

            json_dict['reads files'][local_id] = []

            samp_file_parts = samp_line['reads_file_parts'].split(',')
            samp_line_copy = samp_line.copy()
            for part_id in samp_file_parts:
                samp_line_copy['part_id'] = part_id
                symlink_src = os.path.join(samp_line['reads_dir'],'%s%s%s'%(samp_line['reads_file_prefix'],part_id,samp_line['reads_file_suffix']))
                symlink_target = os.path.join(dest_dir,reads_symlink_pattern%(samp_line_copy))
                if os.path.exists(symlink_target):
                    if (os.path.realpath(symlink_target) != os.path.abspath(symlink_src)):
                        ans = raw_input('Symlink %s in current directory points to %s but you asked for %s, overwrite symbolic link? y/[n]'%(symlink_target,os.path.realpath(symlink_target),os.path.abspath(symlink_src)))
                        if ans=='y':
                            os.remove(symlink_target)
                            os.symlink(symlink_src,symlink_target)
                else:
                    os.symlink(symlink_src,symlink_target)
                json_dict['reads files'][local_id].append(symlink_target)                        
            pipeline_args.append(('--reads-files',','.join(json_dict['reads files'][local_id])))

        # put all the command line utility args in json_dict as its own dict
        json_dict['pipeline args'] = pipeline_args
                
        # get default align_pipeline.py args
        pipeline_args = ' '.join(['%s="%s"'%(k,v) for k,v in pipeline_args])
        print 'align_pipeline.py --run-name=%s %s --print-args'%(run_name,pipeline_args)
        def_args = Popen('align_pipeline.py --run-name=%s %s --print-args'%(run_name,pipeline_args),shell=True,stdout=PIPE,stderr=PIPE).communicate()[0]

        wb('Creating script...\n')
        script_fn = '%s_aln_pipeline.sh'%run_name
        with open(script_fn,'w') as script_f :
            script_f.write(script_template%{'species':species,'aligner':aln,'def_args':def_args})
            os.chmod(script_f.name,stat.S_IRWXU|stat.S_IRWXG|stat.S_IROTH)

        print end_text%{'script_fn':script_fn}

        wb('Creating parameter file...\n')
        json_fn = '%s_aln_params.json'%run_name
        with open(json_fn,'w') as json_f :
            json.dump(json_dict,json_f,indent=4)

    except KeyboardInterrupt :
        sys.stderr.write('\n')
        error('Script creation interrupted, aborting')
