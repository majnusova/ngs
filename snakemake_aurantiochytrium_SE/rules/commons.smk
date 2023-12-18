def get_fastq_files(wildcards):
    sample_df = SAMPLES.loc[wildcards.sample]
    fastq_list = [
        sample_df["fq1"]
        #sample_df["fq2"],
    ]
    return fastq_list




#s touto funkci pipeline fungovala, ale myslim, ze neni treba os.path.join. kdyz mam v dataframe full path k fq1 i fq2:

#def get_fastq_files(wildcards):
    #sample_df = SAMPLES.loc[wildcards.sample]
    #fastq_list = [
        #os.path.join("data", sample_df["fq1"]),
        #os.path.join("data", sample_df["fq2"]),
    #]
    #return fastq_list