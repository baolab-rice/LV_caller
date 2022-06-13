import sys
import re

def read_from_sam(filename):

    umis = {}
    with open(filename,'r') as f:
        for line in f:
            if not line.lstrip().startswith('@'):
                columns = line.split('\t')
                umi = columns[0].split("_consensus")[0]
                umis[umi] = {}
                umis[umi]['cigar'] = columns[5]
                umis[umi]['start_pos'] = columns[3]
                umis[umi]['read'] = columns[9]
    f.close()
    
    return umis

## The functions for large deletion selection
### ERROR!!
def seperate_cigar(dictname):

    p = re.compile("[A-Z]") # actually only a few characters
    pos = []
    for umi in dictname.keys():
        cigar = dictname[umi]['cigar']
        dictname[umi]['pos'] = []
        pos = dictname[umi]['pos']
        start = 0
        for match in p.finditer(cigar):
            pos.append((cigar[start:match.start()],match.group()))
            start = match.start() + 1
    
    return dictname

def invalidated_seqs(dictname):

    new_dict = dictname.copy()
    umis_to_be_deleted = []
    for umi in new_dict.keys():
        marker = 0 
        if dictname[umi]['cigar'] == '*':
            marker = 1          
        if marker == 0:
            umis_to_be_deleted.append(umi)
    
    for umi in umis_to_be_deleted:
        del new_dict[umi]

    return new_dict

def large_deletions_200(dictname):

    new_dict = dictname.copy()
    umis_to_be_deleted = []
    for umi in new_dict.keys():
        marker = 0 
        for pos in new_dict[umi]['pos']:
            if pos[1] == 'N':
                if int(pos[0]) >= 200:
                    marker = 1          
        if marker == 0:
            umis_to_be_deleted.append(umi)
    
    for umi in umis_to_be_deleted:
        del new_dict[umi]

    return new_dict

def large_deletion_50to200(dictname):
    new_dict = dictname.copy()
    umis_to_be_deleted = []
    for umi in new_dict.keys():
        marker = 0 
        for pos in new_dict[umi]['pos']:
            if pos[1] == 'N' and marker != 2:
                if 50<= int(pos[0]) < 200:
                    marker = 1   
                if int(pos[0]) >= 200:
                    marker = 2                 
        if marker != 1:
            umis_to_be_deleted.append(umi)
    
    for umi in umis_to_be_deleted:
        del new_dict[umi]

    return new_dict

def large_deletion_no_LD(dictname1,dictname2,dictname3,dictname4):

    new_dict = dictname1.copy()
    umis_to_be_deleted = []
    for umi in dictname2.keys():
        umis_to_be_deleted.append(umi)
    for umi in dictname3.keys():
        umis_to_be_deleted.append(umi)
    for umi in dictname4.keys():
        umis_to_be_deleted.append(umi)
    for umi in umis_to_be_deleted:
        del new_dict[umi]
    return new_dict   

def calculate_deletion_pos(dictname,cut_pos):

    for umi in dictname.keys():

        start = int(dictname[umi]['start_pos'])
        match =  re.findall('N',dictname[umi]['cigar'])

        if len(match) < 2:
            for pos in dictname[umi]['pos']:
                if re.match('[MD=X]',pos[1]):
                    start = int(pos[0]) + start
                elif pos[1] == 'N':
                    end = int(pos[0]) + start
                    length = int(pos[0]) 
                    break
            dictname[umi]['deletion_start'] = start
            dictname[umi]['deletion_end'] = end
            dictname[umi]['length'] = length
            dictname[umi]['alt'] = 'None'
            if start <= cut_pos and end >= cut_pos:
                dictname[umi]['if_cover_cutsite'] = "yes"
            else:
                dictname[umi]['if_cover_cutsite'] = "no"
        else:
            # if no deletion can cover the cut site, then take the start of the first deletion and the end of the last deletion
            n = len(match)
            start_alt = [None]*n
            end_alt = [None]*n
            i = 0
            for pos in dictname[umi]['pos']:
                if re.match('[MD=X]',pos[1]):
                    start = int(pos[0]) + start
                elif pos[1] == 'N':
                    end = int(pos[0]) + start
                    end_alt[i] = end
                    start_alt[i] = start

                    start = int(pos[0]) + start
                    i += 1
            dictname[umi]['if_cover_cutsite'] = "no"
            dictname[umi]['deletion_start'] = start_alt[0]
            dictname[umi]['deletion_end'] = end_alt[n-1]
            dictname[umi]['length'] = end_alt[n-1] - start_alt[0]
            for u in range(n):
                # if there's a deletion can cover the cut site, then output that deletion's start and end
                if start_alt[u] <= cut_pos and end_alt[u] >= cut_pos:
                    dictname[umi]['if_cover_cutsite'] = "yes"
                    dictname[umi]['deletion_start'] = start_alt[u]
                    dictname[umi]['deletion_end'] = end_alt[u]
                    dictname[umi]['length'] = end_alt[u] - start_alt[u]
                    start_del = start_alt[u]
                    end_del = end_alt[u]
            try:
                del start_del,end_del
            except:pass
            mappings = map(lambda x,y: (x,y), start_alt, end_alt)
            dictname[umi]['alt'] = [x for x in mappings]
            
    return dictname

def small_INDEL(dictname, cutsite):
    umis_un = dictname.copy()
    umis_sINDEL = {}
    for umi in dictname.keys():
        start = int(dictname[umi]['start_pos'])
        marker = 0
        for pos in dictname[umi]['pos']:
            if re.match('[MD=NX]',pos[1]):
                start = int(pos[0]) + start
                if start > cutsite and re.match('[D]',pos[1])::
                        umis_sINDEL[umi] = dictname[umi]
                    break
                # Assum insertion only happends at cutsite
                elif start = cutsite and re.match('[ID]',pos[1]):
                    umis_sINDEL[umi] = dictname[umi]
                    break
    
    for umi in umis_sINDEL.keys():
        del umis_un[umi]

    return umis_sINDEL, umis_un

def select_large_deletions(dictname,code):

    cut_pos = int(code)

    dictname = seperate_cigar(dictname)

    dict_non = invalidated_seqs(dictname)

    dict_200 = large_deletions_200(dictname)
    dict_200 = calculate_deletion_pos(dict_200,cut_pos)
    
    dict_50 = large_deletion_50to200(dictname)
    dict_50 = calculate_deletion_pos(dict_50,cut_pos)

    dict_no_LD = large_deletion_no_LD(dictname,dict_200,dict_50,dict_non)

    dict_sINDEL,dict_un = small_INDEL(dict_no_LD, cut_pos)

    
    return dict_200,dict_50,dict_sINDEL,dict_un,dict_non

#Output the final results

def output(dictname):
    print("UMI\tStart\tEnd\tDeletion_length\tIf_cover_cutsite\tRead\tAlternative_deletion")
    for umi in dictname.keys():
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(umi.split("bin")[0],dictname[umi]['deletion_start'],\
            dictname[umi]['deletion_end'],dictname[umi]['length'],dictname[umi]['if_cover_cutsite'],\
            dictname[umi]['read'],dictname[umi]['alt'][:]))

def write_output_large(dictname,filename):

    f = open(filename,'w')
    f.write("UMI\tStart\tEnd\tDeletion_length\tIf_cover_cutsite\tRead\tAlternative_deletion\n")
    for umi in dictname.keys():
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(umi.split("bin")[0],dictname[umi]['deletion_start'],\
            dictname[umi]['deletion_end'],dictname[umi]['length'],dictname[umi]['if_cover_cutsite'],\
            dictname[umi]['read'],dictname[umi]['alt'][:]))
    f.close()

def write_output_small(dictname,filename):

    f =  open(filename,'w')
    f.write("UMI\tRead\n")
    for umi in dictname.keys():
        f.write("{}\t{}\n".format(umi.split("bin")[0],dictname[umi]['read']))
    f.close()

def write_output_invalidated_umis(dictname,filename):

    f =  open(filename,'w')
    f.write("UMI\n")
    for umi in dictname.keys():
        f.write("{}\n".format(umi))
    f.close()

def large_deletion_calling(filename,code,mode):

    umis = read_from_sam(filename)
    dict_200,dict_50,dict_small,dict_un,dict_non = select_large_deletions(umis,code)
    file_200 = filename.split("_alignment")[0] + "_LD200_temp.txt"
    file_50 = filename.split("_alignment")[0] + "_LD50to200_temp.txt"
    file_small = filename.split("_alignment")[0] + "_small_INDELs_temp.txt"
    file_un = filename.split("_alignment")[0] + "_unmodified_temp.txt"
    file_discard = filename.split("_alignment")[0] + "_invalidated_umis.txt"

    write_output_large(dict_200,file_200)
    write_output_large(dict_50,file_50)
    write_output_small(dict_small,file_small)
    write_output_small(dict_un,file_un)
    write_output_invalidated_umis(dict_non,file_discard)

    return len(dict_200)

def main():
    
    large_deletion_calling(sys.argv[1],"2816c10t","all")
    # umis = read_from_sam(sys.argv[1])
    # umis = select_large_deletions(umis)
    # output(umis)

if __name__ == '__main__':
    main()
