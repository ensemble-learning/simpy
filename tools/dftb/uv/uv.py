import numpy as np
PRE = 1.3062974e8
EV2NM = 1239.84193
DELTA = EV2NM/0.03

def gen_uv(x, f, l1):

    f1 = PRE*f/(1e7*DELTA)
    f2 = np.power((1/x - 1/l1)/(1/DELTA), 2)
    f3 = f1*np.exp(-f2)
    return f3

def get_dftb_output():
    params = []
    f = open("EXC.DAT", "r")
    n_skip = 5
    n = 0
    for i in f:
        if n >= 5:
            tokens = i.strip().split()
            params.append(tokens)
        n += 1
    print(params)
    return params

def output(x, y):
    y_max = np.max(y)
    y = y/y_max
    o = open("uv-vis.dat", "w")
    for i in range(len(x)):
        o.write("%12.4f"%x[i])
        o.write("%26.8f"%y[i])
        o.write("\n")
    o.close()

def main():
    x = np.linspace(100, 600, 200)
    y = np.zeros(len(x))

    params = get_dftb_output()
    for i in range(len(params)):
        l0 = float(params[i][0])
        l1 = EV2NM/l0
        f = float(params[i][1])
        s = float(params[i][5])
        #f = f*s*s
        u = gen_uv(x, f, l1)
        y = y + u

    output(x, y)

if __name__ == "__main__":
    main()



