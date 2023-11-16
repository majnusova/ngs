# BRAKER3
mamba install braker3

"""
The config/ directory from AUGUSTUS can be accessed with the variable AUGUSTUS_CONFIG_PATH.
BRAKER3 requires this directory to be in a writable location, so if that is not the case, copy this directory to a writable location, e.g.:
cp -r /home/users/mvladka/miniconda3/envs/braker3/config/ /absolute_path_to_user_writable_directory/
export AUGUSTUS_CONFIG_PATH=/absolute_path_to_user_writable_directory/config

Due to license and distribution restrictions, GeneMark-ETP and ProtHint should be additionally installed for BRAKER3 to fully work.
These packages can be either installed as part of the BRAKER3 environment, or the PATH variable should be configured to point to them.
The GeneMark key should be located in /home/users//mvladka/.gm_key and GENEMARK_PATH should include the path to the GeneMark executables gmes_petap.pl or gmetp.pl.
"""

#check writability of augustus (jestlize je braker3 installed via conda/mamba, je to writable):
cd /home/users/mvladka/miniconda3/envs/braker3/config/ 
ls -l #w = writable

#GENEMARK-ETP do stejneho env: https://github.com/gatech-genemark/GeneMark-ETP/blob/main/INSTALL

git clone https://github.com/gatech-genemark/GeneMark-ETP.git

echo 'export PATH=/home/users/mvladka/braker3_tools/GeneMark-ETP/tools:$PATH' >> ~/.bashrc
source ~/.bashrc
echo $PATH #for verification 

python3 --version 
perl --version

mamba -c bioconda perl-local-lib # local::lib needed for downloading from CPAN perl lib without sudo (yum -y install perl-App-cpanminus.noarch requires to be root)


#perl -MCPAN -Mlocal::lib -e 'CPAN::install(LWP)' - did not work

mamba install -c bioconda perl-app-cpanminus # mozna je zbytecne mit local__lib i cpanminus?

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
cpanm Storable #if err: mamba install -c conda-forge gcc
cpanm Thread::Queue
cpanm YAML
cpanm YAML::XS
cpanm threads


# add path to bashrc if needed
find ~/ -name gmetp.pl 2>/dev/null
export PATH="/home/users/mvladka/braker3_tools/GeneMark-ETP/bin:$PATH"
source ~/.bashrc

#(base) [mvladka@head mvladka]$ gmetp.pl - now it works


#prothint: https://github.com/gatech-genemark/ProtHint#installation

cpanm MCE::Mutex threads YAML Thread::Queue Math::Utils #but you should have it already


git clone https://github.com/gatech-genemark/ProtHint.git
export PATH="/home/users/mvladka/ProtHint/bin:$PATH"

export PATH="/home/users/mvladka/braker3_tools/GeneMark-ETP/bin/gmes:$PATH"
echo 'export PATH="/home/users/mvladka/braker3_tools/GeneMark-ETP/bin/gmes:$PATH"' >> ~/.bashrc
source ~/.bashrc


                                                                                                                         
#RNA-seq (if no prots available):

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








