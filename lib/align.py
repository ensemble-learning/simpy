"""
Calculated Rotation Matrix to align Vector A to Vector B in 3d using
Rodrigues' rotation formula.

octave code:
v1 = [1 1 1]';
v2 = [0 0 1]';
w = cross(v1,v2);
w = w/norm(w);
w_hat = [0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
cos_tht = v1'*v2/norm(v1)/norm(v2);
tht = acos(cos_tht);
R = eye(size(v1,1)) + w_hat*sin(tht) + w_hat^2*(1-cos(tht));
t1 = R*v1;
t1 = t1/norm(t1)*norm(v2);

"""
import math
import numpy as np

def align_axis(v1, v2):
    w = np.cross(v1,v2)
    w = w/np.linalg.norm(w)
    w_hat = np.matrix([[    0, -w[2],  w[1]], 
                       [ w[2],     0, -w[0]],
                       [-w[1],  w[0],     0]])
    cos_tht = np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v2)
    tht = np.arccos(cos_tht)
    R = np.identity(3) + w_hat*np.sin(tht) + w_hat*w_hat*(1-np.cos(tht))
    v_new = np.dot(R, v1)
    print v_new
    return R

def test():
    v1 = np.array([1.0, 1.0, 1.0])
    v2 = np.array([0.0, 0.0, 1.0])
    R = align_axis(v1, v2)
    print R
if __name__ == "__main__":
    test()
