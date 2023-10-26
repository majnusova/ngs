rule featurecounts:
    input:
        bam_files=expand("/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/{sample}/{sample}.bam", sample=SAMPLES.index)
        
    output:
        result="/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/processed/featurecounts/combined_featurecounts.txt"

    params:
        annotation = "/home/users/mvladka/projects/jotnarlogs/data/rhihy/Rhihy1_GeneCatalog_genes_20151219.gtf"

    log:
        "/home/users/mvladka/projects/jotnarlogs/cilia_snake/data/logs/all_featurecounts.log"
    threads:
        80
    shell:
        """
        featureCounts -a {params.annotation} -o {output.result} -p {input.bam_files} &> {log}
        """
