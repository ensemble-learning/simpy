INF = 99999

def orientation(p, q, r):
    val = (q[1] - p[1])*(r[0] - q[0]) -(q[0] - p[0])*(r[1] - q[1])
    if val == 0:
        return 0
    if val > 0:
        return 1
    if val < 0:  
        return 2
    
def onSegement(p, q, r):
    if (q[0] <= max(p[0], r[0])) and (q[0] >= min(p[0], r[0])): 
        if (q[1] <= max(p[1], r[1])) and (q[1] >= min(p[1], r[1])): 
            return True

def doIntersect(p1, q1, p2, q2):
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)
    print(o1, o2, o3, o4)

    if (o1 != o2) and (o3 != o4):
        return True
    if (o1 == 0) and onSegement(p1, p2, q1):
        return True
    if (o2 == 0) and onSegement(p1, q2, q1):
        return True
    if (o3 == 0) and onSegement(p2, p1, q2):
        return True
    if (o4 == 0) and onSegement(p2, q1, q2):
        return True

    return False

def isInside(polygon, n, p):
    extreme = [INF, p[1]]
    counter = 0
    i = 0
    flag = True
    while (flag):
        nextn = (i+1)%n
        if doIntersect(polygon[i], polygon[nextn], p, extreme):
            if orientation(polygon[i], polygon[nextn], p) == 0:
                return onSegement(polygon[i], p, polygon[nextn])
            counter += 1
        i = nextn
        print(i, counter)
        if i == 0:
            flag = False

    if counter%2 == 0:
        return False
    else:
        return True

def main():
    polygon = [[0,0], [10,0], [10,10], [5,15], [0,10]]
    n = len(polygon)
    p = [5,5]
    flag = isInside(polygon, n, p)
    if flag:
        print('Yes')
    else:
        print('No')

main()
