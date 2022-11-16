##PBS -N matlab1
#PBS -q route
#PBS -l nodes=1,walltime=80:00:00,mem=16000mb,gres=matlab
#PBS -M luisgo@umich.edu
#PBS -m abe
#PBS -V
#
echo "I ran on:"
cat $PBS_NODEFILE
#
#cd to your execution directory first
cd /home/luisgo/Finitedifferenceintegralwall/
#
matlab < hernandezcoils2.m
#