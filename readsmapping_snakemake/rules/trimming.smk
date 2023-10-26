rule trimmomatic_pe:
    input:
        reads = get_fastq_files
    output:
        # forward reads
        forward_paired = f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}_FP.fq.gz",
        forward_unpaired = f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}_FU.fq.gz",
        # reverse reads
        reverse_paired = f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}_RP.fq.gz",
        reverse_unpaired = f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{{sample}}/{{sample}}_RU.fq.gz"
    params:
        adapters = "/home/users/mvladka/miniconda3/envs/snakemake/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    log:
        f"/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/logs/{{sample}}_trimmomatic.log"
    
    shell:
        "trimmomatic PE -threads 20 {input.reads} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} ILLUMINACLIP:{params.adapters}:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 &> {log}"