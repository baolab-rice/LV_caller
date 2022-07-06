import sys 

def read_plots(filename,num,filetype):
    if filetype == "SMRT":
        size = 0
    elif filetype == "LongAmp":
        size = 1
    filelist = []
    templist = []
    with open(filename,'r') as f:
        next(f)
        for line in f:
            if int(line.split('\t')[0]) > size:
                newline = line.strip().split('\t')
                filelist.append([int(newline[0]),int(newline[1]),int(newline[2]),int(newline[2])-int(newline[1])])
    f.close()

    filelist.sort(reverse=True)
    if len(filelist) >= 500:
        for line in filelist[:500]:
            templist.append([int(line[3]),float(int(line[1])+int(line[3])/2-num)])
    else:
        for line in filelist:
            templist.append([int(line[3]),float(int(line[1])+int(line[3])/2-num)])
    return templist

    

def clustering(set_s,set_l):
    newplot = []
    set_s_sorted = sorted(set_s,key=lambda x:(x[0],x[1]))
    set_l_sorted = sorted(set_l,key=lambda x:(x[0],x[1]))

    ds = 20
    ls = 20

    for i in range(len(set_s_sorted)):
        for j in range(len(set_l_sorted)):
            if set_s_sorted[i][0]-ds <= set_l_sorted[j][0] <= set_s_sorted[i][0]+ds:
                if set_s_sorted[i][1]-ls <= set_l_sorted[j][1] <= set_s_sorted[i][1]+ls :
                    newplot.append(((set_s_sorted[i][0]+set_l_sorted[j][0])/2,(set_l_sorted[j][1]+set_s_sorted[i][1])/2))
                    set_l_sorted[j][0] = 0
                    set_l_sorted[j][1] = 0
                    set_s_sorted[i][0] = 0
                    set_s_sorted[i][1] = 0

    return newplot,set_s_sorted,set_l_sorted

def main():
    dataset1 = read_plots(sys.argv[1],int(sys.argv[4]),"SMRT")
    dataset2 = read_plots(sys.argv[2],int(sys.argv[4]),"LongAmp")

    # if len(dataset1) < len(dataset2):
    #     dataset2 = dataset2[:len(dataset1)]
    # elif len(dataset1) > len(dataset2):
    #     dataset1 = dataset1[:len(dataset2)]

    newplot,s,l = clustering(dataset1,dataset2)
    fsname = sys.argv[3]+"_SMRT.txt"
    flname = sys.argv[3]+"_LongAmp.txt"
    fs = open(fsname,'w')
    fl = open(flname,'w')
    f = open(sys.argv[3] + '.txt','w')
    f.write("deletion_middle_point_location,Deletion_length\n")
    oc = 0
    sc = 0
    lc = 0
    for plot in newplot:
        if plot[0] != 0:
            f.write("{},{}\n".format(plot[1],plot[0]))
            oc += 1
    fs.write("deletion_middle_point_location,Deletion_length\n")
    for plot in s:
        if plot[0] != 0:
            fs.write("{},{}\n".format(plot[1],plot[0]))
            sc += 1
    fl.write("deletion_middle_point_location,Deletion_length\n")
    for plot in l:
        if plot[0] != 0:
            fl.write("{},{}\n".format(plot[1],plot[0]))    
            lc += 1
    
    print(sc,lc,oc)

if __name__ == '__main__':
    main()

