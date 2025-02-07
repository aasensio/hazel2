import numpy as np
from hazel.model import Model

def synth_model(config_file, f, active=True, verbose=0):
    """
    Synthesize a model using a configuration file and the results of the inversion as
    defined in the output file.

    Args:
        config_file (str): Filename of the configuration file.
        f (struct): Structure with the results of the inversion obtained directly from the HDF5 file.
        active (bool, optional): list of atmospheres to keep from the output. Defaults to True (include all). For example, active=['ph1', 'ch1', 'te1'].
        verbose (int, optional): Verbosity level. Defaults to 0.

    Returns:
        wl: Wavelength axis of the synthesized spectrum.
        stokes: Stokes parameters of the synthesized spectrum.
    """

    conf = f.attrs['configuration']
    for key in conf:
        if 'Topology =' in key:            
            topology = key.split('=')[1][0:-1]
            break
    
    atms = topology.split('->')

    if active is True:
        active = atms

    mod2 = Model(config_file, working_mode='synthesis', verbose=verbose)

    topology = ''
    first = True

    for atm in atms:
        atm = atm.strip()
        if atm in active:

            if verbose > 3:
                print(f'{atm} - active')

            if first:
                topology += f'{atm}'
                first = False
            else:
                topology += f' -> {atm}'

            if 'ph' in atm:                    
                n = len(f[atm]['T'][0, 0, 0, :])
                model_ph = np.zeros((n, 8))
                model_ph[:, 0] = f[atm]['log_tau'][:]
                model_ph[:, 1] = f[atm]['T'][0, 0, 0, :]
                model_ph[:, 2] = -1.0 * np.ones(model_ph.shape[0])
                model_ph[:, 3] = f[atm]['vmic'][0, 0, 0, :]
                model_ph[:, 4] = f[atm]['v'][0, 0, 0, :]
                model_ph[:, 5] = f[atm]['Bx'][0, 0, 0, :]
                model_ph[:, 6] = f[atm]['By'][0, 0, 0, :]
                model_ph[:, 7] = f[atm]['Bz'][0, 0, 0, :]
                
                mod2.atmospheres[atm].set_parameters(model_ph, 1.0, 0.0)
            
            if 'ch' in atm:                        
                Bx1 = f[atm]['Bx'][0, 0, 0, 0]
                By1 = f[atm]['By'][0, 0, 0, 0]
                Bz1 = f[atm]['Bz'][0, 0, 0, 0]
                v1 = f[atm]['v'][0, 0, 0, 0]
                tau1 = f[atm]['tau'][0, 0, 0, 0]
                a1 = f[atm]['a'][0, 0, 0, 0]
                beta1 = f[atm]['beta'][0, 0, 0, 0]
                deltav1 = f[atm]['deltav'][0, 0, 0, 0]

                if verbose >= 3:
                    labels = ['Bx', 'By', 'Bz', 'tau', 'v', 'deltav', 'beta', 'a']                    
                    print(f'{atm} - (values for all cycles)')
                    for l in labels:
                        print(f"  * {l} : {f[atm][l][:].flatten()}")

                model_atm = [Bx1, By1, Bz1, tau1, v1, deltav1, beta1, a1]
                mod2.atmospheres[atm].set_parameters(model_atm, 1.0)

            if 'te' in atm:
                lambda0 = f[atm]['lambda0'][0, 0, 0, 0]
                sigma = f[atm]['sigma'][0, 0, 0, 0]
                depth = f[atm]['depth'][0, 0, 0, 0]
                a = f[atm]['a'][0, 0, 0, 0]

                if verbose >= 3:
                    labels = ['lambda0', 'sigma', 'depth', 'a']
                    print(f'{atm} - (values for all cycles)')
                    for l in labels:
                        print(f"  * {l} : {f[atm][l][:].flatten()}")

                mod2.atmospheres[atm].set_parameters([lambda0, sigma, depth, a], 1.0)

        else:
            mod2.remove_atmosphere(atm)
            if verbose > 3:
                print(f'{atm} - inactive')

    
    mod2.set_topologies([topology])
    mod2.setup()
    mod2.synthesize()

    wl = mod2.spectrum['spec1'].wavelength_axis
    stokes = mod2.spectrum['spec1'].stokes

    del mod2

    return wl, stokes