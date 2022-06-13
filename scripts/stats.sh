#!/bin/bash

while getopts i: flag
do
    case "${flag}" in
        i) pathname=${OPTARG};;
    esac
done

echo "Ref_CCS_reads:"
valid_con="$(python3 '/home/mc131/Desktop/Projects/SMRTseq/LV_calling-main/scripts/ref_ccs_read.py' "${pathname}/outputs/result_invalidated_umis.txt" "${pathname}/outputs/result_stats.txt")"
echo $valid_con


echo "UMI_consensus:"
UMI_con_ori=$(wc -l ${pathname}/outputs/result_consensus.fasta | awk '{print $1}')
echo $UMI_con_ori/2 | bc 

echo "Invalidated_UMI_consensus:"
UMI_con_in=$(wc -l ${pathname}/outputs/result_invalidated_umis.txt | awk '{print $1}')
echo "$UMI_con_in-1" | bc

echo "Ref_UMI_consensus:"
echo "${UMI_con_ori}/2 - $UMI_con_in + 1" |bc

echo "LD200:"
LD200=$(wc -l ${pathname}/outputs/result_LD200.txt | awk '{print $1}')
echo "$LD200 -1" | bc

echo "LD50-200:"
LD50=$(wc -l ${pathname}/outputs/result_LD50to200.txt | awk '{print $1}')
echo "$LD50 -1" | bc

echo "samllINDEL:"
INDEL=$(wc -l ${pathname}/outputs/result_small_INDELs.txt | awk '{print $1}')
echo "$INDEL -1" | bc

echo "unmodified:"
UN=$(wc -l ${pathname}/outputs/result_unmodified.txt | awk '{print $1}')
echo "$UN -1" | bc

echo "LI_with_LD200:"
LI1=$(wc -l ${pathname}/outputs/result_LI_with_LD200.fasta | awk '{print $1}')
echo "$LI1/2" |bc

echo "LI_with_LD50:"
LI2=$(wc -l ${pathname}/outputs/result_LI_with_LD50.fasta | awk '{print $1}')
echo "$LI2/2" |bc

echo "LI_other:"
LI3=$(wc -l ${pathname}/outputs/result_LI_other.fasta | awk '{print $1}')
echo "$LI3/2" |bc

echo "clusters in LD200:"
clusters=$(wc -l ${pathname}/outputs/result_LD200_cluster.txt | awk '{print $1}')
echo "$clusters -1" | bc
