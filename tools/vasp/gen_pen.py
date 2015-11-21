import sys

def main():
    H = 2.0 
    W = 0.2

    a0 = float(sys.argv[1])
    a1 = float(sys.argv[2])
    a2 = a1 + W
    b0 = float(sys.argv[3])
    b1 = float(sys.argv[4])
    b2 = b1 + W

    x = a0
    y = b1

    o = open("PENALTYPOT", "w")

    while( x < a2 ):
        o.write("%9.4f%9.4f%9.2f%9.2f\n"%(
                x, b1, H, W))
        x += W
        
    x = a1
    y = b0

    while( y < b2 ):
        o.write("%9.4f%9.4f%9.2f%9.2f\n"%(
            a1, y, H, W))
        y += W

    o.write("%9.4f%9.4f%9.2f%9.2f\n"%(
            a2, b2, H, W))
    o.close()

if __name__ == "__main__":
    print "pen.py a0 a1 b0 b1"
    if len(sys.argv) >= 5:
        main()
