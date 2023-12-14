for i in $(ls *.fa); do makeblastdb -in $i -dbtype nucl -parse_seqids -out /home/majnusova/all/eustig_genomes/blast_nucl/$i; done
