import sys

def read_fq(filename):
    f = open(filename,'r')
    lines = f.readlines()
    dict_reads = {}
    print(len(lines))
    for i in range(int(len(lines)/4)):
        dict_reads[lines[i*4][1:].strip()] = (lines[i*4+1].strip(),lines[i*4+2].strip(),lines[i*4+3].strip())
    f.close()

    return dict_reads

def read_ids(filename):
    list_reads = []
    with open(filename,'r') as f:
        for line in f:
            list_reads.append(line.strip())
    f.close()
    print(len(list_reads))
    return list_reads

def output(dictname,listname):
    f = open(sys.argv[3],'w')
    #for read in dictname.keys():
    for read in listname:
        title = "@" + read
        f.write(title)
        f.write('\n')
        f.write(dictname[read][0])
        f.write('\n')
        f.write(dictname[read][1])
        f.write('\n')
        f.write(dictname[read][2])
        f.write('\n')
    f.close() 

def main():
    dict_reads = read_fq(sys.argv[1])
    list_reads = read_ids(sys.argv[2])

    # for read in list_reads:
    #     del dict_reads[read]
    
    output(dict_reads,list_reads)

if __name__ == '__main__':
    main()

