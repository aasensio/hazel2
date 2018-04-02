# cdiazbas@iac.es
import numpy as np


# Return the angles in the plane of the sky given angles with respect
# to the vertical for observations on the limb (in degrees!)
def absolute_to_sky(thetaB, chiB):
    thetaB = np.deg2rad(thetaB)
    chiB = np.deg2rad(chiB)

    t1 = np.sin(thetaB) * np.sin(chiB)
    t2 = -np.cos(thetaB)
    t3 = np.sin(thetaB) * np.cos(chiB)

    thetaSky = np.arccos(t3)
    sinthSky = np.sqrt(1.e0 - t3**2)

    sinChiSky = t1 / sinthSky
    cosChiSky = t2 / sinthSky

# Test for the quadrant
    chiSky_preliminary = np.arccos(cosChiSky)
    if (np.sign(sinChiSky) > 0.e0):
        chiSky = chiSky_preliminary
    else:
        chiSky = -chiSky_preliminary

    return [np.rad2deg(thetaSky), np.rad2deg(chiSky)]


# Return the angles in the vertical system given angles in the
# plane of the sky for observations on the limb (in degrees!)
def sky_to_absolute(thetaSky, chiSky):
    thetaSky = np.deg2rad(thetaSky)
    chiSky = np.deg2rad(chiSky)

    t1 = np.sin(thetaSky) * np.sin(chiSky)
    t2 = np.cos(thetaSky)
    t3 = -np.sin(thetaSky) * np.cos(chiSky)

    thetaB = np.arccos(t3)
    sinthB = np.sqrt(1.e0 - t3**2)

    sinChiB = t1 / sinthB
    cosChiB = t2 / sinthB

# Test for the quadrant
    chiB_preliminary = np.arccos(cosChiB)
    if (np.sign(sinChiB) > 0.e0):
        chiB = chiB_preliminary
    else:
        chiB = -chiB_preliminary

    return [np.rad2deg(thetaB), np.rad2deg(chiB)]


# Return the angles in the plane of the sky given angles with respect
# to the vertical for observations at angle theta (in degrees!)
def absolute_to_sky_general(theta, thetaB, chiB):
    theta = np.deg2rad(theta)
    thetaB = np.deg2rad(thetaB)
    chiB = np.deg2rad(chiB)

    cosThetaSky = np.cos(theta) * np.cos(thetaB) + \
        np.sin(theta) * np.sin(thetaB) * np.cos(chiB)
    sinThetaSky = np.sqrt(1.e0 - cosThetaSky**2)

    thetaSky = np.arccos(cosThetaSky)

    cosChiSky = (np.cos(theta) * np.sin(thetaB) * np.cos(chiB) -
                 np.cos(thetaB) * np.sin(theta)) / sinThetaSky
    sinChiSky = (np.sin(thetaB) * np.sin(chiB)) / sinThetaSky

# Test for the quadrant
    chiSky_preliminary = np.arccos(cosChiSky)
    if (np.sign(sinChiSky) > 0.e0):
        chiSky = chiSky_preliminary
    else:
        chiSky = -chiSky_preliminary

    return [np.rad2deg(thetaSky), np.rad2deg(chiSky)]


# Return the angles in the plane of the sky given angles with respect
# to the vertical for observations at angle theta (in degrees!)
def sky_to_absolute_general(theta, thetaSky, chiSky):
    theta = np.deg2rad(theta)
    thetaSky = np.deg2rad(thetaSky)
    chiSky = np.deg2rad(chiSky)

    cosThetaB = np.cos(theta) * np.cos(thetaSky) - \
        np.sin(theta) * np.sin(thetaSky) * np.cos(chiSky)
    sinThetaB = np.sqrt(1.e0 - cosThetaB**2)

    thetaB = np.arccos(cosThetaB)

    cosChiB = (np.cos(theta) * np.sin(thetaSky) * np.cos(chiSky) +
               np.cos(thetaSky) * np.sin(theta)) / sinThetaB
    sinChiB = (np.sin(thetaSky) * np.sin(chiSky)) / sinThetaB

# Test for the quadrant
    chiB_preliminary = np.arccos(cosChiB)
    if (np.sign(sinChiB) > 0.e0):
        chiB = chiB_preliminary
    else:
        chiB = -chiB_preliminary

    return [np.rad2deg(thetaB), np.rad2deg(chiB)]


if __name__ == '__main__':

    pass
