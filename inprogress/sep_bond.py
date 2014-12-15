import sys

def get_slice_single(t0, t1, dt):
    
    l0 = t0/dt * 2211849
    l1 = t1/dt * 2211849
    
    #-------------------parse bonds -----------------
    if 1:
        print "get the bonds slice"
        f = open("rdxshock.bonds", "r")
        for i in f:
            if "Timestep %d"%t0 in i:
                break
    
        o = open("step%06d.bonds"%t0, "w")
        for i in f:
            if "Timestep %d"%t1 in i:
                break
            else:
                o.write(i)
    
        o.close()
        f.close()
    
    #-------------------parse coords----------------
    if 0:
        print "get the coords slice"
        f = open("dump.shock", "r")
        o = open("dump_%06d.lammpstrj"%t0, "w")
    
        print l0, l1
        counter = 0
        for i in f:
            if counter < l0 :
                pass
            elif counter >= l1:
                break
            else:
                o.write(i)
            counter += 1
    
        o.close()
        f.close()


def main():
    dt = 4000
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    for i in range(a, b, 4000):
        t0 = i
        t1 = i + 4000
        print "processing %d out of 460000"%i
        get_slice_single(t0, t1, dt)

if __name__ == "__main__":
    main()
