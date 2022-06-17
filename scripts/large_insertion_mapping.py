from subprocess import Popen, PIPE
import sys

def blat_exe(filename,ref):

    output = filename.replace('.fasta','.psl')
    arguments = ['blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 {} {} {}'.format(ref,filename,output)]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)

def formatting_psl(filename):

    inputfile = filename.replace('.fasta','.psl')
    temp = inputfile.replace('.psl','_temp.psl')
    output = filename.replace('.fasta','_besthit.txt')
    #'NR==FNR{a[$10]++;next}{print a[$10],$0}'

    a1 = "awk 'NR==FNR{a[$10]++;next}{print a[$10],$0}' " + inputfile + " " +inputfile
    arguments = [a1 + "| sort -Vr | sed 's/[0-9]* //' > {}".format(temp)]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)    

    a1 = "awk '!a[$10]++' " + temp
    a2 = "| awk -v OFS='\t' 'BEGIN{print \"UMI\tChr\tMatch\tMismatch\tStrand\tStart\tEnd\"} ; $10 ~ /^umi*/ {print$10,$14,$1,$2,$9,$16,$17}' >" + output
    arguments = [a1+a2]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)     

    arguments = ["rm {}".format(temp)]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr) 

def LI_mapping(filename, ref):

    file_LI_LD200 = filename + "_LI_with_LD200.fasta"
    file_LI_LD50 = filename + "_LI_with_LD50.fasta"
    file_LI_other = filename + "_LI_other.fasta"

    blat_exe(file_LI_LD200,ref)
    blat_exe(file_LI_LD50,ref)
    blat_exe(file_LI_other,ref)

    formatting_psl(file_LI_LD200)
    formatting_psl(file_LI_LD50)
    formatting_psl(file_LI_other)

def main():
    read_mapping(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
    main()

    


