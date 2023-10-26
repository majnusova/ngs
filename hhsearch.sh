#!/bin/sh

 
cd /home/users/mvladka/projects/microsporidia/data_hhsuite/msa

source activate /home/users/mvladka/miniconda3/envs/hhsuite

ffindex_build -s ../mydatabase.ff{data,index} .
cd ..

ffindex_apply mydatabase.ff{data,index} -i mydatabase_a3m_wo_ss.ffindex -d mydatabase_a3m_wo_ss.ffdata -- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0

rm mydatabase.ff{data,index}

mv mydatabase_a3m_wo_ss.ffindex mydatabase_a3m.ffindex
mv mydatabase_a3m_wo_ss.ffdata mydatabase_a3m.ffdata

ffindex_apply mydatabase_a3m.ff{data,index} -i mydatabase_hhm.ffindex -d mydatabase_hhm.ffdata -- hhmake -i stdin -o stdout -v 0 

cstranslate -f -x 0.3 -c 4 -I a3m -i mydatabase_a3m -o mydatabase_cs219 

sort -k3 -n -r mydatabase_cs219.ffindex | cut -f1 > sorting.dat
    
ffindex_order sorting.dat mydatabase_hhm.ff{data,index} mydatabase_hhm_ordered.ff{data,index}
mv mydatabase_hhm_ordered.ffindex mydatabase_hhm.ffindex
mv mydatabase_hhm_ordered.ffdata mydatabase_hhm.ffdata
    
ffindex_order sorting.dat mydatabase_a3m.ff{data,index} mydatabase_a3m_ordered.ff{data,index}
mv mydatabase_a3m_ordered.ffindex mydatabase_a3m.ffindex
mv mydatabase_a3m_ordered.ffdata mydatabase_a3m.ffdata


for query in $(ls *.fasta); do perl /home/users/mvladka/miniconda3/envs/hhsuite/scripts/reformat.pl fas a3m $query $(basename $query .fasta).a3m -M 50 &> $(basename $query .fasta).reformatlog; done

for a3m_query in $(ls *.a3m); do hhmake -i $a3m_query -M a3m -add_cons &> $(basename a3m_query .a3m).hhmakelog; done
for hhm in $(ls *.hhm); do hhsearch -i $hhm -d mydatabase -M a3m &> $(basename $hhm .hhm).hhsearchlog; done

