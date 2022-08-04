#!/usr/bin/env python3
# The HDR mode is for large integration (>=200bp) with donor.
from os import rename
from subprocess import Popen, PIPE
from alignment_UMIs import long_read_alignment_minimap2
from large_insertions import large_insertion_calling
from large_deletions import large_deletion_calling
from LV_rearrange import rearrange

def rename_file(filename1,filename2):
    arguments = ['mv ' + filename1 + ' > ' + filename2]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)   

def remove_file(filename):
    arguments = ['rm ' + filename]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)     

def partial_ki_identification(infilename,outfilename,code,symbol):

    cut_pos = int(code)
    f = open(outfilename,'w')
    f.write("UMI\tStart\tEnd\n")
    f.close()
    arguments = ['awk ' + "'{if ($2 " + symbol + cut_pos + ") print $1,$2,$3} '"  + infilename + ' >> ' + outfilename]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)       

def HDR_mode_main(umis,reference,filename,ldcode):
    # Initial alignment 
    print("[Initial alignment...]")
    filename_1 = filename + "_initial"
    file_intial_sam = filename_1 + "_alignment.sam"
    long_read_alignment_minimap2(reference,filename_1)

    print("[Large gene modification calling...]")
    large_deletion_calling(file_intial_sam,ldcode,"consensus")
    large_insertion_calling(file_intial_sam)
    rearrange(filename_1)

    ## w/o LI
    # Unmod or premature
    print("[Identifying unmodied and premature reads...]")
    file_unmod = filename + "_HDR_unmod.txt" 
    file_premature = filename + "_HDR_premature.txt"
    file_LD200 = filename_1 + "_LD200.txt"
    partial_ki_identification(file_LD200,file_unmod,ldcode,"<=")
    partial_ki_identification(file_LD200,file_premature,ldcode,">")

    # Perfect TI 
    print("[Identifying perfect TI...]")
    file_small_and_unmod = filename_1 + "_small_INDELs_and_unmod_temp.txt"
    file_perfect_TI = filename + "_perfect_TI.txt"
    rename_file(file_small_and_unmod,file_perfect_TI)

    # remove the file not need
    remove_file(file_LD200)

    ## w/ LI - remap
    print("[Remapping large insertions detected...]")
    file_LI_LD200 = filename + "_LI_with_LD200.fasta"
    file_LI_LD50 = filename + "_LI_with_LD50.fasta"
    file_LI_other = filename + "_LI_other.fasta"

    filename_2_200 = filename + "_remap_LI_with_LD200"
    filename_2_50 = filename + "_remap_LI_with_LD50"
    filename_2_other = filename + "_remap_LI_other"
    long_read_alignment_minimap2(reference,filename_2_200)
    long_read_alignment_minimap2(reference,filename_2_50)
    long_read_alignment_minimap2(reference,filename_2_other)




### 5'TI + 3'HI