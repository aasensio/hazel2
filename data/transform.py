

ff = open('LINEAS', 'r')
flines = ff.readlines()
ff.close()

ff = open('LINEAS.csv', 'w')
for l in flines:
    tmp = l.split()
    idx, element = tmp[0].split('=')

    ff.write('{0},{1},{2},{3}\n'.format(idx, element,tmp[1],tmp[2]))

ff.close()
