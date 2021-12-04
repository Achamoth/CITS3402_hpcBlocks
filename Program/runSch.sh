#PBS -l nodes=4:ppn=1
source /etc/bash.bashrc
mpirun Program/blocks && sleep 20;kill $!
#PBS -m abe
#PBS -M [my email address]
