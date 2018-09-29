import os
import random

nodes_permit = [
#"compute-0-0",
#"compute-0-1",
"compute-0-10",
#"compute-0-2",
#"compute-0-3",
"compute-0-4",
"compute-0-6",
"compute-0-7",
"compute-0-8",
"compute-0-9",
"compute-1-0",
"compute-1-1",
"compute-1-11",
"compute-1-12",
"compute-1-14",
"compute-1-16",
"compute-1-17",
"compute-1-18",
"compute-1-2",
"compute-1-3",
"compute-1-5",
"compute-1-6",
"compute-1-7",
"compute-1-8",
"compute-1-9",
]

def nodes_free( nodes = 8 ):
# This function is used to find available nodes in permitted nodes list

    nodes_list =[]

    nodes_infor = os.popen("freenodes", 'r')
    for i in nodes_infor:
        a = i.strip().split(":")
        if int(a[-1]) >= nodes and a[0] in nodes_permit:
            nodes_list.append(a[0])

    if len(nodes_list) > 0:
        return nodes_list[0]
    else:
        temp = random.randint(0,(len(nodes_permit)-1))
        return nodes_permit[temp]
        #return nodes_free( nodes - 1 )

if __name__ == "__main__":
    print nodes_free()
