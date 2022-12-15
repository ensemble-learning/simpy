def get_charge(fname):
    o = open('q.dat', 'w')
    f = open(fname, 'r')
    for i in f:
        if i.strip().startswith('m_atom'):
            break
    for i in f:
        if i.strip().startswith(':::'):
            break
    for i in f:
        if i.strip().startswith(':::'):
            break
        else:
            tokens = i.strip().split()
            q = tokens[7]
            o.write(q+'\n')
    f.close()
    o.close()

if __name__ == '__main__':
    fname = 'mmod_mini_7.mae'
    get_charge(fname)
