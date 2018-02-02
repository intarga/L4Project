from __future__ import division
import numpy as np
import pickle
import h5py
import astropy.units as u
import eagleSqlTools as sql
import matplotlib.pyplot as plt
from matplotlib import colors

def read_dataset(itype, att, nfiles=16, snapnum=28):
    """ Read a selected dataset, itype is the PartType and att is the attribute name. """
    if snapnum == 28:
        p = '000'
    elif snapnum == 23:
        p = '503'
    elif snapnum == 27:
        p = '101'

    # Output array.
    data = []

    # Loop over each file and extract the data.
    for i in range(nfiles):
        f = h5py.File('RefL0012N0188/snap_0%s_z000p%s.%i.hdf5'%(snapnum, p, i), 'r')
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
    f       = h5py.File('RefL0012N0188/snap_028_z000p000.0.hdf5', 'r')
    a       = f['Header'].attrs.get('Time')         # Scale factor.
    h       = f['Header'].attrs.get('HubbleParam')  # h.
    boxsize = f['Header'].attrs.get('BoxSize')      # L [Mph/h].
    f.close()

    return a, h, boxsize

def read_galaxy(itype, gn, sgn, centre, snapnum=28):
    """ For a given galaxy (defined by its GroupNumber and SubGroupNumber)
    extract the coordinates and mass of all particles of a selected type.
    Coordinates are then wrapped around the centre to account for periodicity. """

    data = {}

    # Load data, then mask to selected GroupNumber and SubGroupNumber.
    gns = read_dataset(itype, 'GroupNumber', snapnum=snapnum)
    sgns = read_dataset(itype, 'SubGroupNumber', snapnum=snapnum)
    mask = np.logical_and(gns == gn, sgns == sgn)
    if itype == 1:
        # data['mass'] = read_dataset_dm_mass()[mask] * u.g.to(u.Msun)
        print "dark matter requested!"
    else:
        data['mass'] = read_dataset(itype, 'Mass', snapnum=snapnum)[mask] * u.g.to(u.Msun)
    if itype == 4:
        data['bd'] = read_dataset(itype, 'BirthDensity', snapnum=snapnum)[mask]
        data['particleID'] = read_dataset(itype, 'ParticleIDs', snapnum=snapnum)[mask]
    data['coords'] = read_dataset(itype, 'Coordinates', snapnum=snapnum)[mask] * u.cm.to(u.Mpc)

    # Periodic wrap coordinates around centre.
    a, h, init_boxsize = read_header()
    boxsize = init_boxsize / h
    data['coords'] = np.mod(data['coords'] - centre + 0.5 * boxsize, boxsize) + centre - 0.5 * boxsize

    return data

data_request = 0

if data_request == 1:

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

        pickle_out = open('request1.pickle','wb')
        pickle.dump(myData, pickle_out)
        pickle_out.close()

elif data_request == 0:
    pickle_in = open('request1.pickle', 'rb')
    myData = pickle.load(pickle_in)

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


pid1 = np.empty(0)
pid2 = np.empty(0)
r_plot = np.empty(0)

MasterPid = read_dataset(4, 'ParticleIDs', snapnum=27)
#print MasterPid
MasterGN = read_dataset(4, 'GroupNumber', snapnum=27)
MasterSGN = read_dataset(4, 'SubGroupNumber', snapnum=27)
mask = np.argsort(MasterPid)
MasterPid = MasterPid[mask]
MasterGN = MasterGN[mask]
MasterSGN = MasterSGN[mask]

n = len(GroupNum)

for i in xrange(1):
    x = read_galaxy(4, GroupNum[i], SubGroupNum[i], GalaxyCentres[i])
    #print i/n
    print 'gn,', GroupNum[i]
    r_temp = np.linalg.norm(x['coords'] - GalaxyCentres[i], axis=1)
    pid_temp = x['particleID']
    delarr1 = np.where(r_temp < 0.003)
    delarr2 = np.where(r_temp > 0.003)
    r1_temp = np.delete(r_temp, delarr1)
    r2_temp = np.delete(r_temp, delarr2)
    pid1_temp = np.delete(pid_temp, delarr1)
    pid2_temp = np.delete(pid_temp, delarr2)
    pid1 = np.concatenate((pid1, pid1_temp))
    pid2 = np.concatenate((pid2, pid2_temp))
    r_plot = np.concatenate((r_plot, r_temp))

    parent_gn_sgn1 = np.zeros((len(pid1),2))
    parent_gn_sgn2 = np.zeros((len(pid2),2))
    #parent_sgn1 = np.zeros(len(pid1))
    #parent_sgn2 = np.zeros(len(pid2))

    for j in xrange(len(pid1)):
        k = np.searchsorted(MasterPid, pid1[j])
        parent_gn_sgn1[j] = [MasterGN[k],MasterSGN[k]]
        #parent_sgn1[j] = MasterSGN[k]

    for j in xrange(len(pid2)):
        k = np.searchsorted(MasterPid, pid2[j])
        parent_gn_sgn2[j] = [MasterGN[k],MasterSGN[k]]
        #parent_sgn2[j] = MasterSGN[k]

    #plt.scatter(parent_gn1, pid1, color='blue')
    #plt.scatter(parent_gn2, pid2, color='red')

    uniques = np.unique(np.vstack((parent_gn_sgn1, parent_gn_sgn2)), axis=0)
    print uniques
    N_bars = len(uniques)
    print np.where(parent_gn_sgn1 == uniques[3])
    plotdata1 = np.zeros(N_bars)
    plotdata2 = np.zeros(N_bars)

    for j in xrange(N_bars):
        plotdata1[j] = np.log10(len(np.where(parent_gn_sgn1 == uniques[j])[0]))
        plotdata2[j] = np.log10(len(np.where(parent_gn_sgn2 == uniques[j])[0]))

    uniques1, plotdata1_temp = np.unique(parent_gn_sgn1, return_counts=True, axis=0)
    uniques2, plotdata2_temp = np.unique(parent_gn_sgn2, return_counts=True, axis=0)

    #for j in xrange(N_bars):
    #    if (uniques[j] != uniques1[j]).all():
    #        uniques1 = np.insert(uniques1, j, [0,0])
    #        plotdata1_temp = np.insert(plotdata1_temp, j, 0)
    #        print 'inserted1'
    #    if (uniques[j] != uniques2[j]).all():
    #        uniques2 = np.insert(uniques2, j, [0,0])
    #        plotdata2_temp = np.insert(plotdata2_temp, j, 0)
    #        print 'inserted2'
    #if len(plotdata1_temp) == N_bars:
    #    print 'yay1'
    #if len(plotdata2_temp) == N_bars:
    #    print 'yay2'

    print 'uniques!', uniques, uniques1, uniques2





ind = np.arange(N_bars)  # the x locations for the groups
width = 0.8  # the width of the bars: can also be len(x) sequence
print plotdata1_temp

p1 = plt.bar(ind, plotdata1, width, color='blue')
p2 = plt.bar(ind, plotdata2, width, color='red', bottom=plotdata1)

plt.xlabel('Parent Galaxy /GN.SGN')
plt.ylabel('Log10[N]')
plt.xticks(ind, uniques, rotation=50, fontsize = 6)
#plt.yticks(np.arange(0, 81, 10))
plt.legend((p1[0], p2[0]), ('Disk', 'Bulge'))

plt.savefig('ParticleTracer.png')