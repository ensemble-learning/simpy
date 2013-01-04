def parse_top(filename):
    top = {}
    f = open(filename, 'r')
    counter = True
    while counter:
        counter = False
        for i in f:
            counter = True
            if i.startswith('#'):
                top['itp'] = i
            elif i.strip().startswith('['):
                key = i.split()[1]
                break
        top[key] = []
        for i in f:
            if len(i.strip()) == 0:
                break
            elif i.startswith(';'):
                pass
            else:
                top[key].append(i)
    return top['atoms']

def main():
    filename = 'out.top'
    parse_top(filename)


if __name__ == "__main__":
    main()


