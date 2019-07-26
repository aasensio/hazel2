import numpy as np
import glob

def transform(filein, fileout):

    tmp = np.loadtxt(filein, skiprows=1)
    n = len(tmp[:,0])
    f = open(fileout, 'w')
    f.write('ff\n')
    f.write('1.0\n\n')
    f.write('  logtau     T        Pe           vmic        v            Bx           By         Bz\n')
    for i in range(n):
        f.write('{0:6.3f}   {1:7.1f}   {2:e}   {3:7.2f}    {4:7.2f}    {5:7.2f}     {6:7.2f}    {7:7.2f}\n'.format(tmp[i,0], tmp[i,1], tmp[i,2], tmp[i,3]*1e-5, tmp[i,5], 0.0, 0.0, 0.0))
    f.close()

files = glob.glob('*.mod')

for f in files:
    print(f)
    tmp = f.split('.')
    transform(f, '{0}.1d'.format(tmp[0]))



