steps = []

f = open("H2.csv", "r")

for i in f:
    tokens = i.strip().split(",")
    steps.append(int(tokens[0]))

f.close()

start = steps[0] 
final = steps[-1]

counter = 0
real_step = 0
real_steps = []
skip = 1000
f = open("water.bboost", "r")

for i in f:
    if counter == 0:
        pass
    else:
        tokens = i.strip().split()
        if len(tokens) == 10:
            nstep = int(tokens[0])
            real_step += float(tokens[7])
            if nstep >= start and nstep <= final:
                if (counter-1)%skip == 0:
                    real_steps.append([nstep, real_step])
    counter += 1

f.close()

o = open("real_time.csv", "w")
for i in real_steps:
    o.write("%-15d,%30.3f\n"%(int(i[0]), i[1]))
o.close()
