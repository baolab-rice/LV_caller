import sys
from SMRTseq_data_processing import read_raw_fasta

def read_inputs(filename):
    dictname = {}
    with open(filename,'r') as f:
        lines = f.readlines()
        lines = lines[1:]
        for line in lines:
            umi = line.split('\t')[0]
            dictname[umi] = line.strip()
    f.close()
    return dictname

def rearrange(filename):

    file_LI = filename + "_LI_temp.fasta"
    file_200 = filename + "_LD200_temp.txt"
    file_50 = filename + "_LD50to200_temp.txt"
    file_small_and_unmod = filename + "_small_INDELs_and_unmod_temp.txt"

    umis_LI = read_raw_fasta(file_LI)
    umis_LD200 = read_inputs(file_200)
    umis_LD50 = read_inputs(file_50)
    umis_small = read_inputs(file_small_and_unmod)

    LD200_and_LI = {}
    LD50_and_LI = {}
    LI_others = {}
    for umi in umis_LI.keys():
        if umi in umis_LD200.keys():
            LD200_and_LI[umi] = umis_LI[umi]
            del umis_LD200[umi]
        if umi in umis_LD50.keys():
            LD50_and_LI[umi] = umis_LI[umi]
            del umis_LD50[umi]
        if umi in umis_small.keys():
            LI_others[umi] = umis_LI[umi]
            del umis_small[umi]

    file_LI_LD200 = filename + "_LI_with_LD200.fasta"
    file_LI_LD50 = filename + "_LI_with_LD50.fasta"
    file_LI_other = filename + "_LI_other.fasta"
    file_200_ = filename + "_LD200.txt"
    file_50_ = filename + "_LD50to200.txt"
    file_small_ = filename + "_small_INDELs_and_unmod.txt"

    f = open(file_200_,'w')
    f.write("UMI\tStart\tEnd\tDeletion_length\tIf_cover_cutsite\tRead\tAlternative_deletion\n")
    for umi,line in umis_LD200.items():
        f.write(line)
        f.write('\n')
    f.close()

    f = open(file_50_,'w')
    f.write("UMI\tStart\tEnd\tDeletion_length\tIf_cover_cutsite\tRead\tAlternative_deletion\n")
    for umi,line in umis_LD50.items():
        f.write(line)
        f.write('\n') 
    f.close()

    f = open(file_small_,'w')
    f.write("UMI\tRead\n")
    for umi,line in umis_small.items():
        f.write(line)
        f.write('\n') 
    f.close()
  

    f = open(file_LI_LD200,'w')
    for umi in LD200_and_LI.keys():
        f.write(">{}\n".format(umi))
        f.write(LD200_and_LI[umi])
        f.write('\n')
    f.close()

    f = open(file_LI_LD50,'w')
    for umi in LD50_and_LI.keys():
        f.write(">{}\n".format(umi))
        f.write(LD50_and_LI[umi])
        f.write('\n')
    f.close()

    f = open(file_LI_other,'w')
    for umi in LI_others.keys():
        f.write(">{}\n".format(umi))
        f.write(LI_others[umi])
        f.write('\n')
    f.close()

    #SA test 20220805
    file_LD = filename + "_LD.txt"
    f = open(file_LD,'w')
    f.write("UMI\tStart\tEnd\tDeletion_length\tIf_cover_cutsite\tRead\tAlternative_deletion\n")
    for umi,line in umis_LD200.items():
        f.write(line)
        f.write('\n')
    for umi,line in umis_LD50.items():
        f.write(line)
        f.write('\n') 
    f.close()

def main():
    rearrange(sys.argv[1][:-9])

if __name__ == '__main__':
    main()