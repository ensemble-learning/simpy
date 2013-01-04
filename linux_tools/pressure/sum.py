f = open("out.log" ,'r')


matrix_p = []
counter = -1
title = "test"

for i in f:
    if i.strip()[:3] == title:
        matrix_p[counter].append(float(i.split()[1]))
        matrix_p[counter].append(float(i.split()[2]))
        matrix_p[counter].append(float(i.split()[3]))
    else:
        counter += 1
        matrix_p.append([])
        matrix_p[counter].append(float(i.split()[1]))
        matrix_p[counter].append(float(i.split()[2]))
        matrix_p[counter].append(float(i.split()[3]))
        title = i.strip()[:3]
f.close()

sum =[[],[],[]]
for i in range(len(matrix_p)):
    for j in range(len(matrix_p[i])):
        if j % 3 == 0:
            matrix_p[i][j] = matrix_p[i][j] / matrix_p[i][90] 
            sum[0].append(matrix_p[i][j])
        elif j % 3 == 1:
            matrix_p[i][j] = matrix_p[i][j] / matrix_p[i][91]
            sum[1].append(matrix_p[i][j])
        elif j % 3 == 2:
            matrix_p[i][j] = matrix_p[i][j] / matrix_p[i][92]
            sum[2].append(matrix_p[i][j])

pre = [[], [], []]
pre[0] = 31*[0]
pre[1] = 31*[0]
pre[2] = 31*[0]

for i in range(len(sum)):
    for j in range(len(sum[i])):
        pre[i][j%31] = pre[i][j%31] + sum[i][j]

o = open("pre_out.log" ,'w')
for i in range(31):
    for j in range(3):
        o.write("%10.4f"%pre[j][i])
    o.write("\n")
o.close()
