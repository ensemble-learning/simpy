import os

def file_input(x, y):
    os.system("mkdir a%03d_%03d"%(x,y))
    os.system("cp * ./a%03d_%03d"%(x,y))

def change_itp(x, y):
    os.chdir("./a%03d_%03d"%(x,y))
    f = open("or.itp", 'r')
    o = open("ethanol.itp", 'w')
    
    a  = 1

    for i in f:
        if a == 8:
            sigma = float(i[45:52])
            epsilin = float(i[53:64])
            o.write(i[:45]+'%7.4f'%(sigma*x/100)+'%12.4f'%(epsilin*y/100)+'\n')
        #if a == 11:
        #    sigma = float(i[45:52])
        #    epsilin = float(i[53:64])
        #    o.write(i[:45]+'%7.4f'%(sigma*x/100)+'%12.4f'%(epsilin*y/100)+'\n')
        elif a == 12:
            sigma = float(i[45:52])
            epsilin = float(i[53:64])
            o.write(i[:45]+'%7.4f'%(sigma*x/100)+'%12.4f'%(epsilin*y/100)+'\n')
        else:
            o.write(i)
        a += 1
    os.chdir("../")

for a in range(95,106):
    for b in range(95,106):
        file_input( a , b )
        change_itp(a , b)
