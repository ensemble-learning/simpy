import paramiko
import os

nodes = [
"node-4-1",
"node-4-2",
"node-4-3",
"node-4-4",
"node-4-5",
"node-4-6",
"node-4-7",
"node-4-9",
"node-4-10",
]

def connect_server(servername, cmd):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(servername, 22, timeout=4)
    print "processing in node: %s"%servername
    #stdin, stdout, stderr = client.exec_command('cd %s; %s'%(path, cmd))
    try:
        stdin, stdout, stderr = client.exec_command('%s'%cmd)
    except SshException:
        pass
    except paramiko.NoValidConnectionsError:
        pass
    except NoValidConnectionsError:
        pass
    except NoValidConnectionsError(errors):
        pass
    except paramiko.ssh_exception.NoValidConnectionsError:
        pass
    except paramiko.ssh_exception.NoValidConnectionsError, exception:
        pass
    lines = stdout.readlines()
    #print lines
    client.close()

CMD = "nvidia-smi"
for i in range(len(nodes)):
    node = nodes[i]
    connect_server(node, CMD)

