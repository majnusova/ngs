rule featurecounts:
    input:
        bam_files=expand("/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{sample}/{sample}.bam", sample=SAMPLES.index)
        
    output:
        result="/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/featurecounts/aurantiochytrium_combined_featurecounts.txt"

    params:
        annotation = "/home/users/mvladka/projects/jotnarlogs/data/aurantiochytrium/Aurli1_GeneCatalog_genes_20120618.gff"

    log:
        "/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/logs/aurantiochytrium_featurecounts.log"
    threads:
        80
    shell:
        """
        featureCounts -a {params.annotation} -s 2 -g name -o {output.result} {input.bam_files} &> {log}
        """