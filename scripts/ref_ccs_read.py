import sys

def read_umi(filename):
    list_reads = []
    with open(filename,'r') as f:
        for line in f:
            list_reads.append(line.strip())
    f.close()
    return list_reads

def read_stats2(filename,listname):
    f = open(filename,'r')
    lines = f.readlines()
    dict_reads = {}
    counts = 0
    for line in lines[1:]:
        UMI_ID = line.split('\t')[0]
        if UMI_ID not in listname:
            count = int(line.split('\t')[1])
            counts = counts + count
    f.close()
    print(counts)

def main():
    umis_del = read_umi(sys.argv[1])
    read_stats2(sys.argv[2],umis_del)

if __name__ == '__main__':
    main()