import hazel
import glob
import os

def test_file_generators():

    tmp = hazel.tools.File_observation(mode='single')
    tmp.set_size(n_lambda=128, n_pixel=1)
    tmp.save('test')

    tmp = hazel.tools.File_observation(mode='multi')
    tmp.set_size(n_lambda=128, n_pixel=10)
    tmp.save('test2')

    tmp = hazel.tools.File_photosphere(mode='single')
    tmp.set_default(n_pixel=1)
    tmp.save('photosphere')

    tmp = hazel.tools.File_photosphere(mode='multi')
    tmp.set_default(n_pixel=10)
    tmp.save('photosphere2')

    tmp = hazel.tools.File_chromosphere(mode='single')
    tmp.set_default(n_pixel=1)
    tmp.save('chromosphere')

    tmp = hazel.tools.File_chromosphere(mode='multi')
    tmp.set_default(n_pixel=10)
    tmp.save('chromosphere2')

    try:
        for f in glob.glob('test*.*'):
            os.remove(f)
    except:
        pass

    try:
        for f in glob.glob('photosphere*.*'):
            os.remove(f)
    except:
        pass

    try:
        for f in glob.glob('chromosphere*.*'):
            os.remove(f)
    except:
        pass