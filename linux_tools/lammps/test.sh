#PBS -S /bin/bash
#PBS -l nodes=1,walltime=5000:00:00
 
cd /home/liulc/src/mycode/lammps
echo $HOSTNAME > host.log
mkdir /state/partition1/liulc/test
cp lammps_input /state/partition1/liulc/test/
cp data_lammps /state/partition1/liulc/test/
cp ffield /state/partition1/liulc/test/
cp lammps_restart /state/partition1/liulc/test/
cd /state/partition1/liulc/test/
lammps < lammps_input 
 
cp * /home/liulc/src/mycode/lammps
if [ $status != 0 ]
then
        echo "Copying results from `hostname`:/state/partition1/liulc/test/ back to $PBS_O_WORKDIR failed." 
        echo "After fixing the problem be sure to remove the directory `hostname`:/state/partition1/liulc/test/"
else
        cd  /home/liulc/src/mycode/lammps
        rm -rf /state/partition1/liulc/test/
fi
 
