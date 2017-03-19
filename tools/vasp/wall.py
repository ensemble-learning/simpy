o = open("HILLSPOT", "w")

cv1 = 5.0
cv2 = 5.0

hd = 0.1
hh = 0.5

x = 0
step = 0.05

while(x-cv1 <= step):
    o.write("%8.4f%8.4f%8.4f%8.4f\n"%(x, cv2, hd, hh))
    x += step

x = 0
while(x-cv2 <= step):
    o.write("%8.4f%8.4f%8.4f%8.4f\n"%(cv1, x, hd, hh))
    x += step

o.close()

