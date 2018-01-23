from __future__ import division
import numpy as np
import pickle
import h5py
import astropy.units as u
import eagleSqlTools as sql
import matplotlib.pyplot as plt
from matplotlib import colors
from pyquaternion import Quaternion
from mpl_toolkits.mplot3d import Axes3D
from astropy import constants as const
import math
from BulgeDiskSeparator4 import *

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def GalaxyRotate(GNs, SGNs, Centres, selection=0):
    '''for a set of galaxies specified by group and subgroup numbers, computes
        histograms of the circularities of their particles, plots them for each galaxy,
        and puts the plots together'''
    plt.figure(figsize=(8, 14))
    n = len(GNs)
    for i in xrange(n):
        # loading galaxy info
        gas = read_galaxy(0, GNs[i], SGNs[i], Centres[i])
        dm = read_galaxy(1, GNs[i], SGNs[i], Centres[i])
        stars = read_galaxy(4, GNs[i], SGNs[i], Centres[i])
        bh = read_galaxy(5, GNs[i], SGNs[i], Centres[i])
        print 'Progressii:', i / n

        # separate data
        r = stars['coords'] - Centres[i]
        bd_temp = stars['bd']
        v = stars['velocity']
        m = stars['mass']
        # print np.sum(m), StellarMass[i]
        N_particles = len(m)

        # Finding axis of galaxy rotation
        L = compute_L(r, m, v, N_particles)
        Galaxy_RotationAxis = compute_RotationAxis(L)
        r_new = np.zeros((N_particles, 3))


        temp1 = np.cross(Galaxy_RotationAxis, (0,0,1))
        temp2 = np.linalg.norm(temp1)
        axis = temp1/temp2
        theta = np.arcsin(temp2)

        rotation_matrixi = rotation_matrix(axis, theta)

        for j in xrange(N_particles):
            r_new[j] = np.dot(rotation_matrixi, r[j])

        print r_new

        x = r_new[:,0]
        y = r_new[:,1]
        z = r_new[:,2]
        print len(z)
        x = [value for value in x if not math.isnan(value)]
        y = [value for value in y if not math.isnan(value)]
        z = [value for value in z if not math.isnan(value)]
        print np.sum(x), np.sum(y)

        print 'completed loop'

        plt.subplot(n, 3, (3 * i) + 1, aspect='equal')
        plt.hist2d(x, y, (300,300), cmap='Greys', norm=colors.LogNorm())
        plt.xlim(-0.1, 0.1)
        plt.ylim(-0.1, 0.1)
        plt.ylabel('GN:'+str(int(GNs[i]))+' '+'SGN:'+str(int(SGNs[i])))#, rotation=80)
        if i+1 == n:
            plt.xlabel('xy plane /mpc')

        plt.subplot(n, 3, (3 * i) + 2, aspect='equal')
        plt.hist2d(y, z, (300, 300), cmap='Greys', norm=colors.LogNorm())
        plt.xlim(-0.1, 0.1)
        plt.ylim(-0.1, 0.1)
        if i + 1 == n:
            plt.xlabel('yz plane / mpc')

        plt.subplot(n, 3, (3 * i) + 3, aspect='equal')
        plt.hist2d(z, x, (300, 300), cmap='Greys', norm=colors.LogNorm())
        plt.xlim(-0.1, 0.1)
        plt.ylim(-0.1, 0.1)
        if i + 1 == n:
            plt.xlabel('zx plane / mpc')

    #plt.tight_layout()
    plt.savefig('GalaxyRotator.png')

    return 0



myQuery = 'SELECT \
                        SH.MassType_Star, \
                        SH.GalaxyID, \
                        SH.GroupNumber, \
                        SH.SubGroupNumber, \
                        SH.CentreOfPotential_x, \
                        SH.CentreOfPotential_y, \
                        SH.CentreOfPotential_z \
                    FROM \
                        RefL0012N0188_SubHalo as SH \
                    WHERE \
                        SH.SnapNum = 28 \
                        and SH.MassType_star > 1e10 \
                    '
myData = read_galaxy_database(myQuery, 'BulgeDiskSeparator4.pickle')

GroupNum, SubGroupNum, GalaxyCentres = extract_galaxydata(myData)

GalaxyRotate(GroupNum, SubGroupNum, GalaxyCentres)

