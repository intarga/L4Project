from __future__ import division
import numpy as np
import h5py
import astropy.units as u
import eagleSqlTools as sql
import matplotlib.pyplot as plt
from matplotlib import colors

def read_dataset(itype, att, nfiles=16):
    """ Read a selected dataset, itype is the PartType and att is the attribute name. """

    # Output array.
    data = []

    # Loop over each file and extract the data.
    for i in range(nfiles):
        f = h5py.File('./data/snap_028_z000p000.%i.hdf5'%i, 'r')
        tmp = f['PartType%i/%s'%(itype, att)][...]
        data.append(tmp)

        # Get conversion factors.
        cgs     = f['PartType%i/%s'%(itype, att)].attrs.get('CGSConversionFactor')
        aexp    = f['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
        hexp    = f['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')

        # Get expansion factor and Hubble parameter from the header.
        a       = f['Header'].attrs.get('Time')
        h       = f['Header'].attrs.get('HubbleParam')

        f.close()

    # Combine to a single array.
    if len(tmp.shape) > 1:
        data = np.vstack(data)
    else:
        data = np.concatenate(data)

    # Convert to physical.
    if data.dtype != np.int32 and data.dtype != np.int64:
        data = np.multiply(data, cgs * a**aexp * h**hexp, dtype='f8')

    return data

def read_header():
    """ Read various attributes from the header group. """
    f       = h5py.File('./data/snap_028_z000p000.0.hdf5', 'r')
    a       = f['Header'].attrs.get('Time')         # Scale factor.
    h       = f['Header'].attrs.get('HubbleParam')  # h.
    boxsize = f['Header'].attrs.get('BoxSize')      # L [Mph/h].
    f.close()

    return a, h, boxsize

def read_galaxy(itype, gn, sgn, centre):
    """ For a given galaxy (defined by its GroupNumber and SubGroupNumber)
    extract the coordinates and mass of all particles of a selected type.
    Coordinates are then wrapped around the centre to account for periodicity. """

    data = {}

    # Load data, then mask to selected GroupNumber and SubGroupNumber.
    gns = read_dataset(itype, 'GroupNumber')
    sgns = read_dataset(itype, 'SubGroupNumber')
    mask = np.logical_and(gns == gn, sgns == sgn)
    if itype == 1:
        # data['mass'] = read_dataset_dm_mass()[mask] * u.g.to(u.Msun)
        print "dark matter requested!"
    else:
        data['mass'] = read_dataset(itype, 'Mass')[mask] * u.g.to(u.Msun)
    if itype == 4:
        data['bd'] = read_dataset(itype, 'StellarFormationTime')[mask]
    data['coords'] = read_dataset(itype, 'Coordinates')[mask] * u.cm.to(u.Mpc)

    # Periodic wrap coordinates around centre.
    a, h, init_boxsize = read_header()
    boxsize = init_boxsize / h
    data['coords'] = np.mod(data['coords'] - centre + 0.5 * boxsize, boxsize) + centre - 0.5 * boxsize

    return data


mySims = np.array([('RefL0012N0188', 12.)])
con = sql.connect('aabraham', password='LM277HBz')

for sim_name, sim_size in mySims:
    print sim_name

    myQuery = 'SELECT \
	            AP.Mass_Star, \
                AP.SFR, \
                SH.MassType_Star, \
                SH.GalaxyID, \
                SH.GroupNumber, \
                SH.SubGroupNumber, \
                SH.CentreOfPotential_x, \
                SH.CentreOfPotential_y, \
                SH.CentreOfPotential_z \
            FROM \
                RefL0012N0188_SubHalo as SH, \
                RefL0012N0188_Aperture as AP \
            WHERE \
	            AP.GalaxyID = SH.GalaxyID \
                and SH.SnapNum = 28 \
	            and AP.ApertureSize = 3 \
            '

    myData = sql.execute_query(con, myQuery)

    #print myData

BulgeMass = np.zeros(len(myData))
MassFraction = np.zeros(len(myData))
sSFR = np.zeros(len(myData))
Index = np.zeros(len(myData))
GroupNum = np.zeros(len(myData))
SubGroupNum = np.zeros(len(myData))
GalaxyCentres = np.zeros((len(myData), 3))

for i in xrange(len(myData)):
    #print(myData[i][3])
    if np.log10(myData[i][2]) > 9.5 and np.log10(myData[i][0]) > 0:
        BulgeMass[i] = np.log10(myData[i][0])
        MassFraction[i] = myData[i][0] / myData[i][2]
        Index[i] = myData[i][3]
        sSFR[i] = np.log10(myData[i][1] / myData[i][0])
        GroupNum[i] = myData[i][4]
        SubGroupNum[i] = myData[i][5]
        GalaxyCentres[i][0] = myData[i][6]
        GalaxyCentres[i][1] = myData[i][7]
        GalaxyCentres[i][2] = myData[i][8]

delarr = np.where(BulgeMass == 0)
BulgeMass = np.delete(BulgeMass, delarr)
MassFraction = np.delete(MassFraction, delarr)
Index = np.delete(Index, delarr)
sSFR = np.delete(sSFR, delarr)
GroupNum = np.delete(GroupNum, delarr)
SubGroupNum = np.delete(SubGroupNum, delarr)
xCentre = GalaxyCentres[:,0]
yCentre = GalaxyCentres[:,1]
zCentre = GalaxyCentres[:,2]
xCentre = np.delete(xCentre, delarr)
yCentre = np.delete(yCentre, delarr)
zCentre = np.delete(zCentre, delarr)
GalaxyCentres = np.column_stack((xCentre,yCentre,zCentre))

bd1 = np.empty(0)
bd2 = np.empty(0)
bd = np.empty(0)
r_plot = np.empty(0)

n = len(GroupNum)

for i in xrange(len(GroupNum)):
    x = read_galaxy(4, GroupNum[i], SubGroupNum[i], GalaxyCentres[i])
    #print i/n
    #print x['coords'], x['mass']
    #print x['coords'][0][1]
    #print x['coords'], 'gc', GalaxyCentres[i]
    r_temp = np.linalg.norm(x['coords'] - GalaxyCentres[i], axis=1)
    bd_temp = x['bd']
    #mask = np.argsort(r)
    #r = r[mask]
    #print r
    #b = x['bd'][mask]
    delarr1 = np.where(r_temp < 0.003)
    delarr2 = np.where(r_temp > 0.003)
    r1_temp = np.delete(r_temp, delarr1)
    r2_temp = np.delete(r_temp, delarr2)
    bd1_temp = np.delete(bd_temp, delarr1)
    bd2_temp = np.delete(bd_temp, delarr2)
    #plt.hist(np.log10(r1_temp), bins=50, histtype='step', color='red', normed=1)
    #plt.hist(np.log10(r2_temp), bins=50, histtype='step', color='blue', normed=1)
    #plt.hist(np.log10(bd1_temp), bins=50, histtype='step', color='blue', normed=1)
    #plt.hist(np.log10(bd2_temp), bins=50, histtype='step', color='red', normed=1)
    #plt.plot(r1_temp, bd1_temp)
    #plt.plot(r2_temp, bd2_temp)
    #plt.xlabel('Log10[Gas Density /gcm^-3]')
    #plt.ylabel('pdf')
    #print 'plotted'
    #plt.savefig('BulgeBDH1.png')
    bd1 = np.concatenate((bd1, bd1_temp))
    bd2 = np.concatenate((bd2, bd2_temp))
    bd = np.concatenate((bd, bd_temp))
    r_plot = np.concatenate((r_plot, r_temp))

plotchoice = 0

if plotchoice == 0:
    plt.hist(bd1, bins=50, histtype='step', color='blue', normed=1)
    plt.hist(bd2, bins=50, histtype='step', color='red', normed=1)
    plt.hist(bd, bins=50, histtype='step', color='black', normed=1)
    plt.xlabel('Birth Expansion Factor')
    plt.ylabel('Number Density of Particles')
    #print 'plotted'
    plt.savefig('BulgeDiskBEH.png')

if plotchoice == 1:
    plt.hist2d(np.log10(r_plot), np.log10(bd), (500,500), norm=colors.LogNorm(), cmap='Greys')
    plt.colorbar()
    plt.axvline(x=np.log10(0.003), linewidth=1, color='red')
    plt.xlabel('Log10[Radius /Mpc')
    plt.ylabel('Log10[Birth Gas Density /gcm^-3]')
    # print 'plotted'
    plt.savefig('BulgeDiskBEH.png')