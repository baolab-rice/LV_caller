#!/usr/bin/env python3

from alignment_UMIs import generate_fasta,generate_fasta_all, long_read_alignment_minimap2


from subprocess import Popen, PIPE
import sys 
import argparse
import re

## Check for the inputs
# Check the input folder, should be the folder of longumi_read pipeline original output 
def check_input_folder():

    if args.directory[-1] == '/':
        args.directory = args.directory[:-1]

# Check the output file name, if it's a folder, create files named with output in that folder
def check_output_filename():

    if args.output[-1] == '/':
        args.output = args.output + "outputs"
    else:
        args.output = args.output + "/outputs"

    arguments = ['mkdir {}'.format(args.output)]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stdout,stderr)

    args.output = args.output + '/result'

# Check the names of required and optional inputs
def check_before_processing():
    check_input_folder()
    check_output_filename()

# Check if anything data during processing is empty
def check_if_empty(element,name):

    if len(element) < 1:
        raise Exception("Your {} is empty, please check your input.".format(name))

### Extract UMI ID, reads ID, centroids and consensus read sequence from longumi_read ouputs

## Extract UMI bins names
def get_umi_bin_list(_dir):

    raconx = "/raconx" + args.raconx + "/"
    arguments = ['ls {}{} | xargs'.format(_dir,raconx)]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)

    umibins_raw = stdout.decode("utf-8")
    umibins_names = [x.strip() for x in umibins_raw.split(" ") if "umi" in x]
    umibins_dirs = [_dir + raconx + x.strip() for x in umibins_raw.split(" ") if "umi" in x]

    return umibins_names,umibins_dirs

## Extract reads IDs, centroids seq and consensus seq for a specific umibin 

def read_from_overlaps(filename):

    read_list = []
    with open(filename,'r') as f:
        for line in f:
            read = line.split("\t")[0]
            read_list.append(read)
        f.close()

    return read_list

def read_from_centroids(filename):

    centroid = ["",""]
    with open(filename,'r') as f:
        for line in f:
            if ">" in line:
                centroid[0] = line.split(";")[0][1:]
            else:
                centroid[1] = centroid[1] + line.strip()
        f.close()

    return centroid

def read_from_consensus_seq(filename):

    consensus_seq = ""
    with open(filename,'r') as f:
        for line in f:
            if ">" not in line and line != "":
                consensus_seq = line.strip()
        f.close()
  
    return consensus_seq   

def read_from_centroids_counts(filename):

    counts = []
    with open(filename,'r') as f:
        for line in f:
            if ">" in line:
                counts.append(line.split("size=")[1][:-2])
        f.close()
    
    count = 0
    for number in counts:
        count += int(number) 
    return count

# meta function for read umi bin
def read_umi_bin(binname,bindir,umis):

    umi = binname.split("bins")[0]

    umis[umi] = {}

    overlaps = bindir + "/ovlp.paf"
    centroids = bindir + "/" + binname + "_centroids.fa"
    consensus_seq = bindir + "/" + binname + "_sr.fa"
    count = read_from_centroids_counts(centroids)

    readids = read_from_overlaps(overlaps)
    centroids = read_from_centroids(centroids)
    cons_seq = read_from_consensus_seq(consensus_seq)

    umis[umi]['readids'] = readids
    umis[umi]['centroids'] = centroids
    umis[umi]['cons_seq'] = cons_seq
    umis[umi]['count'] = count

    return umis

# Generate a list of polished UMIs
def read_polished_umi(filename):

    polished_umis = []
    with open(filename,'r') as f:
        for line in f:
            if ">" in line:
                umi = line[1:].split("bins")[0]
                polished_umis.append(umi)
    f.close()
    return polished_umis
 
def polishing_umis(umis,polished_list):

    umis_to_delete = []
    for umi in umis.keys():
        if umi not in polished_list:
            umis_to_delete.append(umi)
    for umi in umis_to_delete:
        del umis[umi]
    return umis

## Convert raw fasta file into a dictionary
def read_raw_fasta(filename):
    
    reads = {}
    with open(filename,'r') as f:
        for line in f:
            if ">" in line:
                title = line[1:].strip().split(" ")[0]
            else:
                reads[title] = line.strip()
        f.close()

    return reads

## map reads with read ids from umi bins 
def map_read_with_id(readids,raw_reads):

    temp_dict = {}
    for readid in readids:
        try:
            temp_dict[readid] = raw_reads[readid]
        except:pass

    return temp_dict

## Output the file with    UMI_ID#  #reads/bin consensus_read_sequence 
def output_reads(filename,umis):
    filename = filename + ".txt"
    f =  open(filename,'w')
    # Output 1: UMI_ID Read_ID Read_sequence Centroid_ID
    if args.output_style == '1' and args.reads != None:
        f.write("UMI_ID\tRead_ID\tRead_sequence\tCentroid_ID\n")
        for umi in umis.keys():
            # when only one read in bin, Read_sequence = Consensus_read_sequence = Centroid
            if umis[umi]['readids'] == {}:
                f.write("{}\t{}\t{}\t{}\n"\
                        .format(umi,umis[umi]['centroids'][0],umis[umi]['centroids'][1]\
                        ,umis[umi]['centroids'][0]))
            else:
                for read in umis[umi]['readids'].keys():
                    read_seq = umis[umi]['readids'][read]
                    consensus_seq = umis[umi]['cons_seq']
                    centroid_ID = umis[umi]['centroids'][0]
                    centroid = umis[umi]['centroids'][1]
                    f.write("{}\t{}\t{}\t{}\n"\
                            .format(umi,read,read_seq,centroid_ID))

    # Output2:
    #        UMI
    #        Consensus_seq
    #        Reads_seqs
    if args.output_style == '2' and args.reads != None:
        for umi in umis.keys():
            if umis[umi]['readids'] == {}:
                f.write(">{}\n".format(umi))
                f.write("{}\n".format(umis[umi]['centroids'][0]))
                f.write("{}\n".format(umis[umi]['centroids'][1]))
            else:
                consensus_seq = umis[umi]['cons_seq']
                centroid_ID = umis[umi]['centroids'][0]
                centroid = umis[umi]['centroids'][1]
                f.write("{}\n".format(umi))
                f.write("Consensus sequence:{}\n".format(consensus_seq))
                for read in umis[umi]['readids'].keys():
                    read_seq = umis[umi]['readids'][read]
                    f.write("{}\n".format(read+";"+read_seq))
    f.close()    

def output_stats(filename,umis):
    # stats1: UMI_ID Reads Consensus_read_sequence
    if args.stats == '1':
        f = open(filename + "_stats.txt",'w')
        f.write("UMI_ID\tRead_IDs\tConsensus_read_sequence\n")
        for umi in umis.keys():
            if umis[umi]['readids'] == {}:
                f.write("{}\t{}\t{}\n"\
                        .format(umi,umis[umi]['centroids'][0],umis[umi]['centroids'][1]))
            else:
                reads = ""
                for read in umis[umi]['readids']:
                    reads = reads + read + ","
                f.write("{}\t{}\t{}\n"\
                        .format(umi,reads[:-1],umis[umi]['cons_seq']))   
        f.close()     

    # stats2:  UMI_ID Count Consensus_read_sequence
    if args.stats == '2':
        f = open(filename + "_stats.txt",'w')
        f.write("UMI_ID\tRead_count\tConsensus_read_sequence\n")
        for umi in umis.keys():
            if umis[umi]['readids'] == {}:
                f.write("{}\t{}\t{}\n"\
                        .format(umi,1,umis[umi]['centroids'][1]))
            else:
                f.write("{}\t{}\t{}\n"\
                        .format(umi,umis[umi]['count'],umis[umi]['cons_seq']))   
        f.close() 
        
def get_umis():

    # Variable define
    umis = {}

    # Read from raw files
    print("[Gathering umibins list...]")
    umibins_names,umibins_dirs = get_umi_bin_list(args.directory)
    check_if_empty(umibins_dirs,"UMI List")

    print("[Reading from umibins...]")
    for i in range(len(umibins_names)):
        umis = read_umi_bin(umibins_names[i],umibins_dirs[i],umis)
    
    print("[Polishing umis...]")
    polished_cons_filename = args.directory + "/consensus_raconx3.fa"
    polished_umis = read_polished_umi(polished_cons_filename)
    umis = polishing_umis(umis,polished_umis)
     
    return umis

def output_UMIs(umis):

    # General output and stats
    if args.reads != None:
        print("[Reading from raw fasta file...]")
        reads = read_raw_fasta(args.reads)

        print("[Mapping reads with reads IDs...]")
        for name,value in umis.items():
            umis[name]['readids'] = map_read_with_id(value['readids'],reads)
    if args.output_style != 0:
        output_reads(args.output,umis)
    
    print("[Writing into files...]")
    output_stats(args.output,umis)

def generate_filtered_fasta(umis):
    print("[Generating fasta file for consensus seqs...]")
    if args.all == True:
        generate_fasta_all(umis,args.output)
    else:
        generate_fasta(umis,args.output)

def alignment():
    print("[Alignment using minimap2...]")
    inputfile = args.output + "_consensus.fasta"
    long_read_alignment_minimap2(args.reference,inputfile,args.output,"longread")

def HDR_mode(umis):
    print("[Running in HDR mode...]")
    from LV_caller_HDR_mode import HDR_mode_main
    HDR_mode_main(umis,args.reference,args.output,args.large_deletion_parameters)

def large_deletion():

    print("[Large deletions calling...]")
    inputfile = args.output + "_alignment.sam"
    if args.all == True:
        mode = 'all'
    else:
        mode = 'consensus'
    from large_deletions import large_deletion_calling
    LD200_size = large_deletion_calling(inputfile,args.large_deletion_parameters,mode)

    return LD200_size

def large_insertion():

    print("[Large insertions calling...]")
    inputfile = args.output + "_alignment.sam"
    from large_insertions import large_insertion_calling
    large_insertion_calling(inputfile,args.large_deletion_parameters)

def large_variants_rearrange():

    print("[Large deletions and large insertions rearranging...]")
    from LV_rearrange import rearrange
    rearrange(args.output)


def large_deletion_clustering():

    print("[Clustering large deletions (>=200bp)...]")
    inputfile = args.output + "_LD200.txt"
    from clustering import cluster_generate
    cluster_generate(inputfile, args.large_deletion_clustering_parameters)

    #SA test 20220805
    print("[Clustering intermedia deletions (50bp-200bp)...]")
    inputfile = args.output + "_LD50to200.txt"
    cluster_generate(inputfile, args.large_deletion_clustering_parameters)

    print("[Clustering all large deletions (>=50bp)...]")
    inputfile = args.output + "_LD.txt"
    cluster_generate(inputfile, args.large_deletion_clustering_parameters)    

def distribute_LD():
    #Gernate distribution figure
    print("[Generating distribution figure for large deletions (>=200bp)...]")
    from distribution import distribution_generate
    distribution_generate(args.output + "_LD200.txt", int(args.large_deletion_parameters))

    #Gernate distribution figure
    print("[Generating distribution figure for clustered large deletions (>=200bp)...]")
    distribution_generate(args.output + "_LD200_cluster.txt", int(args.large_deletion_parameters))

    #SA test 20220805
    #Gernate distribution figure
    print("[Generating distribution figure for intermediate deletions (50bp-200bp)...]")
    distribution_generate(args.output + "_LD50to200.txt", int(args.large_deletion_parameters))

    #Gernate distribution figure
    print("[Generating distribution figure for clustered intermediate deletions (50bp-200bp)...]")
    distribution_generate(args.output + "_LD50to200_cluster.txt", int(args.large_deletion_parameters))

    #Gernate distribution figure
    print("[Generating distribution figure for all large deletions (>=50bp)...]")
    distribution_generate(args.output + "_LD.txt", int(args.large_deletion_parameters))

    #Gernate distribution figure
    print("[Generating distribution figure for clustered all large deletions (>=50bp)...]")
    distribution_generate(args.output + "_LD_cluster.txt", int(args.large_deletion_parameters))   

def map_LI():
    #Mapping large insertions to reference genome
    print("[Mapping large insertions (>=50bp) to reference genome...]")
    from large_insertion_mapping import LI_mapping
    LI_mapping (args.output, args.genome)

def remove_temp():

    print("[Removing temp files...]")
    arguments = ['rm ' + args.output[:-6] + '*temp*']
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)

def generate_stats():

    print("[Generating stats...]")
    arguments = ['bash stats.sh -i {} > {}'.format(args.output[:-15],args.output + "_output_stats.txt")]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)


def main():

    # Check everything before data processing
    check_before_processing()

    umis = get_umis()
    output_UMIs(umis)
    generate_filtered_fasta(umis)

    if args.mapping == True and args.reference != None and args.HDR_mode == False:
        alignment()
        if args.large_deletion == True:
            LD200_size = large_deletion()
            large_insertion()
            large_variants_rearrange()
            if args.large_deletion_clustering == True and LD200_size != 0:
                large_deletion_clustering()
                distribute_LD()
                remove_temp()
                map_LI()
            generate_stats()  
    elif args.reference != None and args.HDR_mode == True:
        HDR_mode(umis)
    

    print("Program finished.")

if __name__ == "__main__":

    desc = "longumi_read downstraem analysis: Version 3\
            this script is for extracting data from longumi_read pipeline output.\
            You may use 'seqtk seq -a <fastq> > <fasta> to convert file if raw reads need to be processed.'"
    parser = argparse.ArgumentParser(description=desc)
    # Required 
    parser.add_argument('-d','--directory',required=True,help="Directory of outout folder with longumi_read pipeline, \
                                                            the output folder should contain a raconx subfolder.")
    parser.add_argument('-o','--output',required=True, help="Output directory. Output the completly organized file.")

    parser.add_argument('-st','--stats',choices=['1','2'],default='2', help="Output the stats filem (default=2):\n\
                         -st 1: UMI_ID Read_IDs Consensus_read_sequence;\n \
                         -st 2: UMI_ID Read_count Consensus_read_sequence")


    # Optional for reads related output
    parser.add_argument('-os','--output_style', choices=['1','2'],default='2',help="If also involve raw reads, can also produce a file contaning all reads.\n \
                         -os 1: UMI_ID Read_ID Read_sequence Centroid_ID.\n\
                         -os 2: UMI\
                                Consensus_seq\
                                Reads_seqs")
    parser.add_argument('-r','--reads',help="Input PicBio ccs fasta file. (converted with seqtk)")


    # Alignment
    parser.add_argument('-m','--mapping',default=False,action='store_true',\
                        help="Mapping all filetered read to reference amp using minimap2.\
                            For the large deletion analysis option, could ONLY use minimap2.")
    parser.add_argument('-g','--reference',help="Alignment reference amp. If use HDR mode, \
        the expected HDR amplicon sequence is required.")

    # Large deletions calling
    parser.add_argument('-ld','--large_deletion',default=False,action='store_true',\
                        help="Large deletion calling, devide LDs as small INDELs or unmodified, \
                            50bp-200bp, and >=200bp.")
    parser.add_argument('-ld_ps','--large_deletion_parameters',\
                        help="Large deletion analysis parameters:")

    # Large insertions calling 
    parser.add_argument('-li','--large_insertion',default=False,action='store_true',\
                        help="Large insertion (>=50bp) calling.")
    parser.add_argument('-G','--genome',default=False,help="Reference Genome.")

    # Large deletion clustering
    parser.add_argument('-ld_c','--large_deletion_clustering',default=False,action='store_true',\
                        help="Large deletion clustering based on deletion size and deletion start site.")
    parser.add_argument('-ld_cps','--large_deletion_clustering_parameters',default='10d10l',\
                        help="Large deletion clustering parameters: \
                        FORMAT: deletion_size_tolenrance+d+deletion_position_tolerance+l \
                        Default: 10d10l")
    # HDR mode
    parser.add_argument('-hdr','--HDR_mode',default=False,action='store_true',\
                        help="Running in HDR mode. The expected HDR amplicon sequence is required as a reference in -g.")
 
    # For all reads
    parser.add_argument('-a','--all',default=False,action='store_true',\
                        help="Process all reads. ")
    
    parser.add_argument('-rxn','--raconx',default="3",\
                        help="raconx number used in longread_umi pipeline (default:3).")

    args = parser.parse_args()

    main()
