import csv
import math
import os
import optparse
import re
import string
import sys

import shared.LoggerFactory

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
                            help='output file showing mapping of original filename to uniform symbolic links',
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

    logger.info('Reading sample sheet %s',opts.samp_sheet)
    logger.info('Writing symbolic link file mapping to %s',opts.filename_map)
    ifh = open(opts.samp_sheet,'r')
    ofh = open(opts.filename_map,'w')
    samp_reader = csv.DictReader(ifh,dialect='excel-tab')
    samp_writer = csv.DictWriter(ofh,
                                 fieldnames=['Local-ID','Original-Filename','Uniform-Filename'],
                                 dialect='excel-tab')
    samp_writer.writeheader()
    for sample in samp_reader:
        output_path_sub = os.path.join(opts.output_path,
                                  sample['Species'],sample['Genomic-Element-Type'],
                                  sample['Technique'],sample['File-Format'])
        if (not os.path.isdir(output_path_sub)):
                logger.info('Making new directory for this data type: %s',output_path_sub)
                os.makedirs(output_path_sub)
        symlink_target_pref = '%(Target)s:Acc=%(Accession)s#TL=%(Tag-Line)s#Age=%(Age)s#LC=%(Light-Condition)s#Temp=%(Temperature)s#Treatment=%(Treatment)s:%(Technique)s:%(Replicate)s:%(Dataset-Role)s:%(Local-ID)s'%(sample)
        samp_file_parts = sample['reads_file_parts'].split(',')
        for part_id in samp_file_parts:
            symlink_src = os.path.join(sample['reads_dir'],'%s%s%s'%(sample['reads_file_prefix'],part_id,sample['reads_file_suffix']))
            symlink_target = os.path.join(output_path_sub,
                                          '%(pref)s:%(part)s%(suffix)s'%{'pref':symlink_target_pref,
                                                                       'part':part_id,
                                                                       'suffix':sample['reads_file_suffix']})
            if (os.path.islink(symlink_target)):
                if (opts.rewrite):
                    os.remove(symlink_target)
                else:
                    raise OSError('Symlink target %s already exists and rewrite was set to False'%(symlink_target))

            os.symlink(symlink_src,symlink_target)
            samp_writer.writerow({'Local-ID':sample['Local-ID'],
                                  'Original-Filename':symlink_src,
                                  'Uniform-Filename':symlink_target})

    ifh.close()
    ofh.close()
    logger.info('Check symlinks in %s',opts.output_path)
    logger.info('Done!')

if __name__=='__main__':
    main()
