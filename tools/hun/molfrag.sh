
# argument
if [ $# -ne 1 ]
then
  echo "Usage : molfrag.sh [control.file]"
  exit
fi

Logfile="log.molfrag"

echo "Molecular Fragment Analysis Toolkit"
echo "-----------------------------------"
echo ""
echo "Programmed by Hyungjun Kim (linus) at 2007"
echo ""

#echo "Molecular Fragment Analysis Toolkit" > $Logfile
#echo "-----------------------------------" >> $Logfile
#echo "" >> $Logfile
#echo "Programmed by Hyungjun Kim (linus) at 2007" >> $Logfile
#echo "" >> $Logfile

echo " *Creating Bond Change Logs."
#echo " *Creating Bond Change Logs." >> $Logfile
bondlog $1
#cp dataTATB/Bond.log dataTATB/Bond.log.bak
#cp dataTATB/Nbondchg.log dataTATB/Nbondchg.log.bak

echo ""
echo " *Reducing Bond Change Logs."
echo " (Deleting Redundant reactions)"

#echo "" >> $Logfile
#echo " *Reducing Bond Change Logs." >> $Logfile
#echo " (Deleting Redundant reactions)" >> $Logfile
netbondlog $1

echo ""
echo " *Molecular Fragment Analysis"

#echo "" >> $Logfile
#echo " *Molecular Fragment Analysis" >> $Logfile
molfrag $1
molstat $1

MolStatList=$(grep MolStatTotal $1)
MolStatList=${MolStatList#MolStatTotal}

fragtable $1 $MolStatList

echo ""
echo " *Rxn Analysis"

#echo "" >> $Logfile
#echo " *Molecular Fragment Analysis" >> $Logfile
rxnanal $1

echo "DONE"
#rm bondlog netbondlog molfrag rxnanal molstat fragtable
