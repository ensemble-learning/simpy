import paramiko
import os

def connect_server(servername, path, cmd):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(servername, 22, timeout=4)
    print "processing in node: %s"%servername
    print "processing in folder: %s"%path
    stdin, stdout, stderr = client.exec_command('cd %s; %s'%(path, cmd))
    for std in stdout.readlines():
       #print stdout,
       print std,
    client.close()

#CMD = "python /home/chengtao/src/amd_analysis/AMD-analysis-tools/findreaction_fun.py"
#CMD = "head H2O1_reax_table.txt"
#CMD = "tail ./out.csv"
CMD = "ls -l"
#CMD = "grep step10000000 water.mol -A 4"
#CMD = "grep ^10000000 out.csv"
for i in range(20):
    folder = "r%02d"%i
    os.chdir(folder)
    f = open('host.log', 'r')
    fc = f.readlines()
    node = fc[0].strip()
    rf = fc[1].strip()
    connect_server(node, rf, CMD)
    os.chdir('..')

