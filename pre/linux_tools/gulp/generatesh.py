import os
import time

list = ["112", "121", "211", "221", "212", "122", "222", "444"]

for i in list:
    o = open("a%s.sh"%i, 'w')
    o.write("""#!/bin/bash
#PBS -l nodes=1:ppn=1
cd $PBS_O_WORKDIR
/home/chengtao/gulp/Src/gulp < TATB_%s.gin > a%s.log
"""%(i,i))
    o.close()

for i in list:
    time.sleep(3)
    os.system("qsub a%s.sh"%i)

