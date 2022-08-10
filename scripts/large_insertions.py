from large_deletions import read_from_sam, seperate_cigar
import sys
import re

def large_insertions(dictname):

    umis_to_be_deleted = []
    for umi in dictname.keys():
        marker = 0 
        for pos in dictname[umi]['pos']:
            if pos[1] == 'I':
                if int(pos[0]) >= 50:
                    marker = 1     
        if marker == 0:
            umis_to_be_deleted.append(umi)
    
    for umi in umis_to_be_deleted:
        del dictname[umi]

    return dictname

def calculate_insertion_pos(dictname,output,code):

    cutsite = int(code)
    f = open(output,'w')
    #SA test 20220805
    file_umi_cut_pos = output.split("_LI_temp.fasta")[0] + "_LI_around_cutsite.txt"
    f2 = open(file_umi_cut_pos,'w')
    for umi in dictname.keys():

        start = 0
        match =  re.findall('I',dictname[umi]['cigar'])

        for pos in dictname[umi]['pos']:
            if pos[1] == 'I' and int(pos[0]) >= 50:
                end = int(pos[0]) + start
                f.write(">{}\n".format(umi))
                f.write(dictname[umi]['read'][start:end])
                f.write('\n')
                if cutsite - 100 <= start <= cutsite + 100:
                    f2.write(umi)
                    f2.write('\n')
                start = end
            elif re.match('[MS=X]',pos[1]):
                start = int(pos[0]) + start
    f.close()
    f2.close()
            
    return dictname

def select_large_insertions(dictname,output,code):

    dictname = seperate_cigar(dictname)
    dictname = large_insertions(dictname)
    dictname = calculate_insertion_pos(dictname,output,code)
    
    return dictname

def large_insertion_calling(filename,code):

    output = filename.split("_alignment")[0] + "_LI_temp.fasta"
    umis = read_from_sam(filename)
    umis = select_large_insertions(umis,output,code)

def main():
    
    umis = read_from_sam(sys.argv[1])
    umis = select_large_insertions(umis)

if __name__ == '__main__':
    main()

