def get_slice(nstep):
    f = open("dump.%d.cfg"%nstep, "r")
    o = open("slice.%08d.cfg"%nstep, "w")
    
    counter = 0
    for i in f:
        o.write(i)
        if counter >= 15:
            break
        counter += 1
    
    n = 1
    counter = 0
    tmp = []
    for i in f:
        tmp.append(i)
        if n%3 == 0:
            tokens = tmp[2].strip().split()
            z = float(tokens[2])
            if z > 0.1 and z <= 0.12:
                for j in tmp:
                    o.write(j)
                counter += 1
            n = 0
            tmp = []
        n += 1
    o.seek(0)
    o.write("Number of particles = %7d"%counter)
    o.close()
    f.close()

def main():
    t0 = 0
    t1 = 100000
    dt = 4000
    
    for i in range(t0, t1, dt):
        get_slice(i)

if __name__ == "__main__":
    main()
