import hazel
import numpy as np

def test_coordinates():
    """
    Test for magnetic fields given in spherical/cartesian geometries 
    and vertical/line-of-sight reference frames
    """

    # Define LOS theta angle
    theta_los = 25.0

    # Atmosphere with vertical/cartesian
    mod = hazel.Model(working_mode='synthesis')
    mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
        'LOS': [theta_los,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
    mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 
        'Wavelength': [10826, 10833], 'Reference frame' : 'vertical', 'Coordinates for magnetic field vector': 'cartesian'})
    mod.setup()

    # Atmosphere with vertical/spherical
    mod2 = hazel.Model(working_mode='synthesis')
    mod2.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
        'LOS': [theta_los,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
    mod2.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 
        'Wavelength': [10826, 10833], 'Reference frame' : 'vertical', 'Coordinates for magnetic field vector': 'spherical'})
    mod2.setup()

    # Atmosphere with line-of-sight/cartesian
    mod3 = hazel.Model(working_mode='synthesis')
    mod3.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
        'LOS': [theta_los,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
    mod3.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 
        'Wavelength': [10826, 10833], 'Reference frame' : 'line-of-sight', 'Coordinates for magnetic field vector': 'cartesian'})
    mod3.setup()

    # Atmosphere with line-of-sight/spherical
    mod4 = hazel.Model(working_mode='synthesis')
    mod4.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
        'LOS': [theta_los,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
    mod4.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 
        'Wavelength': [10826, 10833], 'Reference frame' : 'line-of-sight', 'Coordinates for magnetic field vector': 'spherical'})
    mod4.setup()

    # Magnetic field points along the LOS
    B = 0.5
    thB = 25.0
    phiB = 0.0

    mod2.atmospheres['ch1'].set_parameters([B, thB, phiB,1.0,0.0,8.0,1.0,0.0],1.0)
    mod2.synthesize()

    # Transform to cartesian and check that they are the same
    Bx = B * np.sin(thB * np.pi/180) * np.cos(phiB * np.pi/180)
    By = B * np.sin(thB * np.pi/180) * np.sin(phiB * np.pi/180)
    Bz = B * np.cos(thB * np.pi/180)
    mod.atmospheres['ch1'].set_parameters([Bx, By, Bz,1.0,0.0,8.0,1.0,0.0],1.0)
    mod.synthesize()

    # Compute now in LOS/spherical and check that they are the same
    # We use thetaB=0 in this case because LOS is pointing along the field.
    mod4.atmospheres['ch1'].set_parameters([B, 0.0, phiB,1.0,0.0,8.0,1.0,0.0],1.0)
    mod4.synthesize()

    np.testing.assert_allclose(mod.spectrum['spec1'].stokes, mod2.spectrum['spec1'].stokes, rtol=1e-5, atol=0)
    np.testing.assert_allclose(mod.spectrum['spec1'].stokes, mod4.spectrum['spec1'].stokes, rtol=1e-5, atol=0)
