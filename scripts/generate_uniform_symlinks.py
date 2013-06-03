import csv
import logging
import math
import os
import optparse
import re
import string
import sys

import shared.LoggerFactory

DEST_DIR_PATTERN = os.path.join('%(Species)s','%(Genomic-Element-Type)s','%(Technique)s','%(File-Format)s')
READS_SYMLINK_PATTERN = '%(Target)s:Acc=%(Accession)s#TL=%(Tag-Line)s#Age=%(Age)s#LC=%(Light-Condition)s#Temp=%(Temperature)s#Treatment=%(Treatment)s:%(Technique)s:%(Replicate)s:%(Dataset-Role)s:%(Local-ID)s:%(part_id)s%(reads_file_suffix)s'

def make_symlink(samp_sheet,output_path,sample_info_keys,make_link=False,rewrite=False):
    # return list of dict of keys '
    logger = logging.getLogger()
    logger.info('Reading sample sheet %s',samp_sheet)
    ifh = open(samp_sheet,'r')
    samp_reader = csv.DictReader(ifh,dialect='excel-tab')
    samp_info = []
    for samp_line in samp_reader:
        output_path_sub = os.path.join(output_path,DEST_DIR_PATTERN%(samp_line))

        if (make_link and not os.path.isdir(output_path_sub)):
                logger.info('Making new directory for this data type: %s',output_path_sub)
                os.makedirs(output_path_sub)

        samp_file_parts = samp_line['reads_file_parts'].split(',')
        for part_id in samp_file_parts:
            symlink_src = os.path.join(samp_line['reads_dir'],'%s%s%s'%(samp_line['reads_file_prefix'],part_id,samp_line['reads_file_suffix']))
            samp_line_copy = samp_line.copy()
            samp_line_copy['part_id'] = part_id
            symlink_target = os.path.join(output_path_sub,READS_SYMLINK_PATTERN%(samp_line_copy))
            if make_link:
                if (os.path.islink(symlink_target)):
                    if (rewrite):
                        os.remove(symlink_target)
                    else:
                        raise OSError('Symlink target %s already exists and rewrite was set to False'%(symlink_target))
                os.symlink(symlink_src,symlink_target)
            samp_line_copy['Original-Filename'] = symlink_src
            samp_line_copy['Uniform-Filename'] = symlink_target
            samp_info.append(dict((k,samp_line_copy.get(k,'')) for k in sample_info_keys))
    ifh.close()
    return samp_info


def main():
    usage = 'usage: %prog [arguments]'
    parser = optparse.OptionParser(usage=usage)
    # required arguments
    optgroup_req = optparse.OptionGroup(parser,'Required arguments')
    optgroup_req.add_option('--sample-sheet',
                            help='sample sheet in tab delimited format; see an example in the examples directory',
                            action='store',dest='samp_sheet')
    optgroup_req.add_option('--output-path',
                            help='directory name for putting the symbolic links',
                            action='store',dest='output_path')
    optgroup_req.add_option('--filename-map',
                            help='output file showing mapping of original filename to uniformly named symbolic links',
                            action='store',dest='filename_map')
    parser.add_option_group(optgroup_req)

    # optional arguments
    optgroup_opt = optparse.OptionGroup(parser,'Optional arguments')
    optgroup_opt.add_option("--rewrite",
                            help='if set, rewrite any existing target symlinks for the files in the current sample sheet',
                            action="store_true", dest="rewrite")
    parser.add_option_group(optgroup_opt)
    (opts,args) = parser.parse_args()

    (opts,args) = parser.parse_args()
    mandatories = ['samp_sheet','output_path','filename_map']
    for m in mandatories:
        if getattr(opts,m) == None:
            print 'Required argument %s is missing\n'%(m)
            parser.print_help()
            sys.exit(-1)

    logger = shared.LoggerFactory.get_logger(opts.filename_map+'.log')
    shared.LoggerFactory.log_command(logger,sys.argv)

    samp_info_keys = ['Local-ID','Original-Filename','Uniform-Filename']
    samp_info = make_symlink(opts.samp_sheet,opts.output_path,samp_info_keys,make_link=True,rewrite=opts.rewrite)
    logger.info('Writing symbolic link file mapping to %s',opts.filename_map)
    ofh = open(opts.filename_map,'w')
    samp_writer = csv.DictWriter(ofh,fieldnames=samp_info_keys,dialect='excel-tab')
    samp_writer.writeheader()
    for samp_line in samp_info:
        samp_writer.writerow(samp_line)

    logger.info('Check symlinks in %s',opts.output_path)
    logger.info('Done!')

if __name__=='__main__':
    main()
