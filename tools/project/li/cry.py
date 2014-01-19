a = 5.76
b = 3.51
for i in range(1,6):
     for j in range(1,12):
         val = abs((a*i - b*j)/(a*i))
         if val < 0.05:
             print i, j, a*i - b*j, a*i, val
         
    
