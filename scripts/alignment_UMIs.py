#!/usr/bin/env python3
from subprocess import Popen, PIPE

## Generating alignment .sam file when -m flag activated
def generate_fasta(dictname,filename):
    
    f = open(filename + "_consensus.fasta","w")
    for umi,info in dictname.items():
        # I realized this step could be done when generating the dictionary
        if info['cons_seq'] == "":
            cons_seq = info['centroids'][1]
        else:
            cons_seq = info['cons_seq']
        f.write(">{}\n".format(umi + "_consensus_seq"))
        f.write(cons_seq + "\n")
    f.close()

def generate_fasta_all(dictname, filename):

    f = open(filename + "_consensus.fasta","w")
    for umi in dictname.keys():
        # when only one read in bin, Read_sequence = Consensus_read_sequence = Centroid
        if dictname[umi]['readids'] == {}:
            f.write(">{}\n{}\n".format(dictname[umi]['centroids'][0],dictname[umi]['centroids'][1]))
        else:
            for read in dictname[umi]['readids'].keys():
                read_seq = dictname[umi]['readids'][read]
                f.write(">{}\n{}\n"\
                        .format(read,read_seq)) 
    f.close()

# Long read alignment with minimap2

def long_read_alignment_minimap2(reference,inputfile,filename,mode):
    if mode == "longread":
        flag = " -ax splice "
    elif mode == "nonlongread":
        flag = " -a "

    if reference == None:
        print("No refence genome for alignment assigned, please input a reference with -r <reference_genome>")
    else:
        outputfile = filename + "_alignment.sam"
        logfile = filename + "_alignment.log"
        arguments = ['minimap2 ' + flag + '{} {} > {}'\
                .format(reference, inputfile, outputfile)]
        process = Popen(args = arguments,
                        shell=True,
                        stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        f = open(logfile,"w")
        f.write(stdout.decode("utf-8"))
        f.write(stderr.decode("utf-8"))
        f.close()
if __name__ == '__main__':
    pass