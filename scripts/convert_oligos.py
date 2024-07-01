#! /usr/bin/env python3
import argparse
import os
import re
import sys

def convert_oligos_to_mappingfile_rj(oligofile,outdir, targetdir=None,verbose=True):
        if targetdir is None:
                targetdir = os.path.dirname(oligofile)
                print("target dir: %s" % targetdir)
        if not re.search('.oligos$',oligofile) and not re.search('.oligos.clean$',oligofile):
                raise Exception('YTError:',oligofile,' may not be a .oligos file??')
        oligobasename = os.path.basename(oligofile)
        single_primer = []
        paired_primer = []
        comment = []
        barcode = []
        with open(oligofile,'r') as ofile:
                lines = ofile.read().splitlines()
        for line in lines:
                if re.search('(?i)^(forward|reverse)',line):
                        single_primer.append(line.split('\t'))
                elif re.search('(?i)^primer',line):
                        paired_primer.append(line.split('\t'))
                elif re.search('(?i)^barcode',line):
                        barcode.append(line.split('\t'))
                elif re.search('^#',line):
                        comment.append(line)
                else:
                        if re.search('^#?\t',line):
                                raise Exception('YTError: Tabs detected!')
                        else:
                                raise Exception('YTError: Not sure how to deal with this line:\n',line)
        if len(paired_primer)==2 and len(single_primer)==0:
                #miseq run
                primers = [p[1] for p in paired_primer]
                revprimers = [p[2] for p in paired_primer]
                taxfiles = [ '1.map.txt', '2.map.txt']
                print(taxfiles)
                barcode_list = []
                for b in barcode:
                        bar = b[1]
                        samp_id = b[3]
                        barcode_list.append((bar,samp_id))
                header = ['#SampleID','BarcodeSequence','LinkerPrimerSequence','ReversePrimer','Description'] #rj added reverse primer
                samp_id_taboo_chars = '[_%+;: -]+'
                for taxfile,primer,revprimer in zip(taxfiles,primers,revprimers):
                        with open(os.path.join(outdir, taxfile),'w') as tf:
                                tf.write('\t'.join(header)+'\n')
                                for bar,samp_id in barcode_list:
                                        if re.search(samp_id_taboo_chars,samp_id):
                                                samp_id_corrected = re.sub(samp_id_taboo_chars,'.',samp_id)
                                                map_line = [samp_id_corrected,bar,primer,revprimer,samp_id]
                                                #print 'Corrected sample ID name: '+samp_id+' -> '+samp_id_corrected
                                        else:
                                                map_line = [samp_id,bar,primer,revprimer,samp_id]
                                        tf.write('\t'.join(map_line)+'\n')
        elif len(single_primer)==1 and len(paired_primer)==0:
                #454
                primers = [p[1] for p in single_primer]
                #taxfiles = [re.sub('.oligos$','.map.txt',os.path.join(targetdir,oligobasename))]
                taxfiles = [outdir + '/map.txt']
                barcode_list = [(b[1],b[2]) for b in barcode]
                header = ['#SampleID','BarcodeSequence','LinkerPrimerSequence','Description'] #rj no reverse primer here
                samp_id_taboo_chars = '[_%+;: -]+'
                for taxfile,primer in zip(taxfiles,primers):
                        with open(os.path.join(outdir, taxfile),'w') as tf:
                                tf.write('\t'.join(header)+'\n')
                                for bar,samp_id in barcode_list:
                                        if re.search(samp_id_taboo_chars,samp_id):
                                                samp_id_corrected = re.sub(samp_id_taboo_chars,'.',samp_id)
                                                map_line = [samp_id_corrected,bar,primer,samp_id]
                                                #print 'Corrected sample ID name: '+samp_id+' -> '+samp_id_corrected
                                        else:
                                                map_line = [samp_id,bar,primer,samp_id]
                                        tf.write('\t'.join(map_line)+'\n')
        else:
                raise Exception('YTError: Is this 454 or miseq??')
        print('Conversion completed: '+oligofile+' converted to:\n'+', '.join(taxfiles))
        return taxfiles

def main():
        parser = argparse.ArgumentParser(description='YT: Converts mothur oligo file to qiime mapping file.')
        parser.add_argument('oligofile',help='mothur oligo file to be converted.')
        args = parser.parse_args()
        convert_oligos_to_mappingfile_rj(args.oligofile)

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    convert_oligos_to_mappingfile_rj(snakemake.input.oligosfile, outdir=snakemake.params.outdir)
