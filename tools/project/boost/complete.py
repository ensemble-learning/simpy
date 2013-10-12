import os

def get_sim_time(n = 66):

    f = open("water.mol", "r")
    for i in f :
        if "step" in i:
             step = int(i.strip().split()[0][4:])
        if "H2O1" in i:
            tokens = i.strip().split()
            nwater = int(tokens[0])
            if nwater >= n:
                break
    f.close()
    
    return step

def get_real_time(nstep):
    
    if not os.path.exists("water.bboost"):
        realtime = nstep
    else:
        counter = 0
        realtime = 0.0
        f = open("water.bboost", "r")
        for i in f:
            if counter > 0:
                tokens = i.strip().split()
                step = int(tokens[0])
                if step > nstep:
                    break
                realtime += float(tokens[7])
            counter += 1
        f.close()
    return realtime

def get_half_time():
    nstep = get_sim_time(33)
    realtime = get_real_time(nstep)
    print nstep, realtime

if __name__ == "__main__":
    get_half_time()
