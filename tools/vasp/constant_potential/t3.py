F_pzc = -1065.7647177216472301 # hartree
N_pzc = 191 # N
Mu_pzc = -0.289507276 # hartree
C = 9.2535 # e/hartree

def cal_g(F_pzc, N_pzc, Mu_pzc, C, Mu):

    n2 = N_pzc
    dn = (Mu - Mu_pzc)*C 
    n1 = n2 + dn
    mu1 = Mu
    mu2 = Mu_pzc

    a = n1 - n2
    b = (n1 - n2)*mu2
    c = mu1 - mu2
    d = 0.5*(n1 -n2)*(mu1-mu2)
    e = n1*mu1
    F = F_pzc + b + d
    G = F - e

    G = F_pzc - 0.5*C*c*c - n2*mu1
    return mu1, G

Mus = [
-0.134453,
-0.148153,
-0.162963,
-0.177258,
-0.191481,
]

G_sim = [
-1040.2019953160327077,
-1037.5643743294467640,
-1034.7153598387719740,
-1031.9678769762815591,
-1029.2367137438520786,
]

n = 0
for i in Mus:
    Mu = i
    mu, G = cal_g(F_pzc, N_pzc, Mu_pzc, C, Mu)
    print mu, G, "%6.2f"%((G-G_sim[n])*627/23)
    n += 1
