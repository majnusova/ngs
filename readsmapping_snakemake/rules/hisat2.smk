rule hisat:
    input:
        i1 = f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}_FP.fq.gz",
        i2 = f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}_RP.fq.gz",
        #i4 = f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}_FU.fq.gz",
        #i3 = f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}_RU.fq.gz",
        
    output:
        f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}.bam"
    params:
        index_file = "/home/users/mvladka/projects/jotnarlogs/intermediate_files/rhihy_index"
    log:
        f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/logs/{{sample}}_hisat.log"
    threads:
        40
    shell:
        "hisat2 -q -x {params.index_file} -1 {input.i1} -2 {input.i2} | samtools sort -o {output} &> {log}" 
        #"hisat2 -x {params.index_file} -1 {input.i1} -2 {input.i2} -U {input.i3},{input.i4} | samtools sort -o {output} &> {log}" = command co zalignuje i unpaired readys