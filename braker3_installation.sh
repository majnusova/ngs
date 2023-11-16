mamba create --name myenvname braker3


#check writability of augustus (jestlize je braker3 installed via conda/mamba, je to writable):
cd /home/users/mvladka/miniconda3/envs/braker3/config/ 
ls -l #w = writable

git clone https://github.com/gatech-genemark/GeneMark-ETP.git

echo 'export PATH=/home/users/mvladka/braker3_tools/GeneMark-ETP/tools:$PATH' >> ~/.bashrc
source ~/.bashrc
echo $PATH #for verification 

#check if you have python and perl (python3 --version) and download tools for cpan perl 

mamba -c bioconda perl-local-lib # local::lib needed for downloading from CPAN perl lib without sudo


#perl -MCPAN -Mlocal::lib -e 'CPAN::install(LWP)' - did not work

mamba install -c bioconda perl-app-cpanminus 

cpanm Cwd
cpanm Data::Dumper
cpanm File::Path
cpanm File::Spec
cpanm File::Temp
cpanm FindBin
cpanm Getopt::Long
cpanm Hash::Merge
cpanm List::Util
cpanm MCE::Mutex
cpanm Math::Utils
cpanm Parallel::ForkManager
cpanm Statistics::LineFit
cpanm Storable #if err: mamba install -c conda-forge gcc)
cpanm Thread::Queue
cpanm YAML
cpanm YAML::XS
cpanm threads


needed to add path to bashrc: (find ~/ -name gmetp.pl 2>/dev/null)
export PATH="/home/users/mvladka/braker3_tools/GeneMark-ETP/bin:$PATH"
source ~/.bashrc

#(base) [mvladka@head mvladka]$ gmetp.pl - now it works


prothint
git clone https://github.com/gatech-genemark/ProtHint.git
export PATH="/home/users/mvladka/ProtHint/bin:$PATH"

export PATH="/home/users/mvladka/braker3_tools/GeneMark-ETP/bin/gmes:$PATH"
echo 'export PATH="/home/users/mvladka/braker3_tools/GeneMark-ETP/bin/gmes:$PATH"' >> ~/.bashrc
source ~/.bashrc


                                                                                                                         
#RNA-seq (no prot available):

#!/bin/sh
#PBS -N braker3
#PBS -q batch
#PBS -l nodes=1:ppn=80
#PBS -l walltime=999:00:00
#PBS -j oe

source activate braker3

braker.pl --species=pseudellipsoidon --genome=/home/users/mvladka/projects/pseudellipsoidon/Pseudellipsoidion_assembly.fasta \
   --rnaseq_sets_ids=pseuedell-edaph \
   --rnaseq_sets_dirs=/home/users/mvladka/projects/pseudellipsoidon &> braker.log








