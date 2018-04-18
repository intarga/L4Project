from __future__ import division
import numpy as np
import pickle
import h5py
import astropy.units as u
import eagleSqlTools as sql
import matplotlib.pyplot as plt
from astropy import constants as const
import pdb

# basedir = '/opt/local/l4astro/tt/L0025N0375/REFERENCE/snapshot_028_z000p000/'
basedir = "RefL0012N0188/"

def read_galaxy_database(query, filename, username, password, simname, data_request=0):
    if data_request == 1:

        mySims = np.array([(simname, 12.)])
        con = sql.connect(username, password=password)

        for sim_name, sim_size in mySims:
            # print sim_name

            myData = sql.execute_query(con, query)

            # print myData

            pickle_out = open(filename, 'wb')
            pickle.dump(myData, pickle_out)
            pickle_out.close()

    elif data_request == 0:
        pickle_in = open(filename, 'rb')
        myData = pickle.load(pickle_in)
        pickle_in.close()

    return myData

def extract_galaxydata(myData, snapnum):
    # Initialising arrays for Galaxy Data
    a, h, init_boxsize, aexp, hexp = read_header(snapnum)
    print myData.shape
    if myData.shape != ():

        myDataLength = len(myData)
        BulgeMass = np.zeros(myDataLength)
        StellarMass = np.zeros(myDataLength)
        sSFR = np.zeros(myDataLength)
        Index = np.zeros(myDataLength)
        GroupNum = np.zeros(myDataLength)
        SubGroupNum = np.zeros(myDataLength)
        GalaxyCentres = np.zeros((myDataLength, 3))
        GalaxyVelocities = np.zeros((myDataLength, 3))
        # ImageURLs = np.zeros(myDataLength, dtype='object')

        # Extracting Galaxy Data from SQL results
        for i in xrange(myDataLength):
            StellarMass[i] = myData[i][0]
            Index[i] = myData[i][1]
            GroupNum[i] = myData[i][2]
            SubGroupNum[i] = myData[i][3]
            GalaxyCentres[i][0] = myData[i][4] * a
            GalaxyCentres[i][1] = myData[i][5] * a
            GalaxyCentres[i][2] = myData[i][6] * a
            GalaxyVelocities[i][0] = myData[i][7] * (u.km/u.s).to(u.Mpc/u.s)
            GalaxyVelocities[i][1] = myData[i][8] * (u.km/u.s).to(u.Mpc/u.s)
            GalaxyVelocities[i][2] = myData[i][9] * (u.km/u.s).to(u.Mpc/u.s)
            # linkend = len(myData[i][9]) - 3
            # ImageURLs[i] = myData[i][9][11:linkend]

    else:
        myData = str(myData)
        myData = myData[2:-1]
        myData = myData.split(",")
        myData = map(str.strip, myData)
        myData = list(map(float, myData))
        print myData
        myDataLength = 1
        BulgeMass = np.zeros(myDataLength)
        StellarMass = np.zeros(myDataLength)
        sSFR = np.zeros(myDataLength)
        Index = np.zeros(myDataLength)
        GroupNum = np.zeros(myDataLength)
        SubGroupNum = np.zeros(myDataLength)
        GalaxyCentres = np.zeros((myDataLength, 3))
        GalaxyVelocities = np.zeros((myDataLength, 3))

        StellarMass[0] = myData[0]
        Index[0] = myData[1]
        GroupNum[0] = myData[2]
        SubGroupNum[0] = myData[3]
        GalaxyCentres[0][0] = myData[4] * a
        GalaxyCentres[0][1] = myData[5] * a
        GalaxyCentres[0][2] = myData[6] * a
        GalaxyVelocities[0][0] = myData[7] * (u.km / u.s).to(u.Mpc / u.s)
        GalaxyVelocities[0][1] = myData[8] * (u.km / u.s).to(u.Mpc / u.s)
        GalaxyVelocities[0][2] = myData[9] * (u.km / u.s).to(u.Mpc / u.s)

    return GroupNum, SubGroupNum, GalaxyCentres, GalaxyVelocities

def read_dataset(itype, att, nfiles=16, snapnum=28):
    """ Read a selected dataset, itype is the PartType and att is the attribute name. """

    z_ref = ["020p000", "015p132", "009p993", "008p988", "008p075",
             "005p050", "005p971", "005p487", "005p037", "004p485",
             "003p984", "003p528", "003p017", "002p478", "002p237",
             "002p012", "001p737", "001p487", "001p259", "001p004",
             "000p865", "000p736", "000p615", "000p503", "000p366",
             "000p271", "000p183", "000p101", "000p000"]

    # Output array.
    data = []

    # Loop over each file and extract the data.
    for i in range(nfiles):
        # print itype, att
        fname = basedir + 'snap_0{}_z{}.{}.hdf5'.format(snapnum, z_ref[int(snapnum)], i)
        f = h5py.File(fname, 'r')
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

def read_dataset_dm_mass(snapnum):
    """ Special case for the mass of dark matter particles. """

    z_ref = ["020p000", "015p132", "009p993", "008p988", "008p075",
             "005p050", "005p971", "005p487", "005p037", "004p485",
             "003p984", "003p528", "003p017", "002p478", "002p237",
             "002p012", "001p737", "001p487", "001p259", "001p004",
             "000p865", "000p736", "000p615", "000p503", "000p366",
             "000p271", "000p183", "000p101", "000p000"]

    fname       = basedir + 'snap_0{}_z{}.0.hdf5'.format(snapnum, z_ref[int(snapnum)])
    f           = h5py.File(fname, 'r')
    h           = f['Header'].attrs.get('HubbleParam')
    a           = f['Header'].attrs.get('Time')
    dm_mass     = f['Header'].attrs.get('MassTable')[1]
    n_particles = f['Header'].attrs.get('NumPart_Total')[1]

    # Create an array of length n_particles each set to dm_mass.
    m = np.ones(n_particles, dtype='f8') * dm_mass

    # Use the conversion factors from the mass entry in the gas particles.
    cgs  = f['PartType0/Mass'].attrs.get('CGSConversionFactor')
    aexp = f['PartType0/Mass'].attrs.get('aexp-scale-exponent')
    hexp = f['PartType0/Mass'].attrs.get('h-scale-exponent')
    f.close()

    # Convert to physical.
    m = np.multiply(m, cgs * a**aexp * h**hexp, dtype='f8')

    return m

def read_header(snapnum):
    """ Read various attributes from the header group. """

    z_ref = ["020p000", "015p132", "009p993", "008p988", "008p075",
             "005p050", "005p971", "005p487", "005p037", "004p485",
             "003p984", "003p528", "003p017", "002p478", "002p237",
             "002p012", "001p737", "001p487", "001p259", "001p004",
             "000p865", "000p736", "000p615", "000p503", "000p366",
             "000p271", "000p183", "000p101", "000p000"]

    fname   = basedir + 'snap_0{}_z{}.0.hdf5'.format(snapnum, z_ref[int(snapnum)])
    f       = h5py.File(fname, 'r')
    a       = f['Header'].attrs.get('Time')         # Scale factor.
    h       = f['Header'].attrs.get('HubbleParam')  # h.
    boxsize = f['Header'].attrs.get('BoxSize')      # L [Mph/h].
    hexp = f['PartType0/Mass'].attrs.get('h-scale-exponent')
    aexp = f['PartType0/Mass'].attrs.get('aexp-scale-exponent')
    f.close()

    return a, h, boxsize, aexp, hexp

def read_galaxy(itype, gn, sgn, centre, snapnum, extras=[]):
    """ For a given galaxy (defined by its GroupNumber and SubGroupNumber)
    extract the coordinates and mass of all particles of a selected type.
    Coordinates are then wrapped around the centre to account for periodicity. """

    data = {}

    # Load data, then mask to selected GroupNumber and SubGroupNumber.
    gns = read_dataset(itype, 'GroupNumber', snapnum=snapnum)
    sgns = read_dataset(itype, 'SubGroupNumber', snapnum=snapnum)
    mask = np.logical_and(gns == gn, sgns == sgn)
    if itype == 1:
        data['mass'] = read_dataset_dm_mass(snapnum)[mask] * u.g.to(u.Msun)
    else:
        data['mass'] = read_dataset(itype, 'Mass', snapnum=snapnum)[mask] * u.g.to(u.Msun)
    if itype == 4:
        #data['bd'] = read_dataset(itype, 'BirthDensity', snapnum=snapnum)[mask]
        data['velocity'] = read_dataset(itype, 'Velocity', snapnum=snapnum)[mask] * (u.cm/u.s).to(u.Mpc/u.s)
        data['ParticleIDs'] = read_dataset(itype, 'ParticleIDs', snapnum=snapnum)[mask]
        for att in extras:
            data[att] = read_dataset(itype, att, snapnum=snapnum)[mask]
    if itype == 5:
        data["BH_Mass"] = read_dataset(itype, "BH_Mass", snapnum=snapnum)[mask]
    data['coords'] = read_dataset(itype, 'Coordinates', snapnum=snapnum)[mask] * u.cm.to(u.Mpc)

    # Periodic wrap coordinates around centre.
    a, h, init_boxsize, aexp, hexp = read_header(snapnum)
    print "init_boxsize", init_boxsize
    boxsize = init_boxsize * a / h
    print boxsize
    data['coords'] = np.mod(data['coords'] - centre + 0.5 * boxsize, boxsize) + centre - 0.5 * boxsize
    return data

def compute_L(r, m, v, N_particles):
    '''for a set of particles with r relative to a common centre, compute
    their angular momentum around that centre.'''
    L = np.zeros((N_particles, 3))
    for j in xrange(N_particles):
        L[j] = np.cross(r[j], m[j]*v[j])

    return L

def compute_RotationAxis(L):
    '''for a set of particles, compute the axis of their collective rotation.'''
    Galaxy_L = np.sum(L, axis=0)
    Galaxy_RotationAxis = Galaxy_L / np.linalg.norm(Galaxy_L)
    #print "rotation axis:", Galaxy_RotationAxis, 'L:', np.linalg.norm(Galaxy_L)

    return Galaxy_RotationAxis

def compute_Lperp(L, Galaxy_RotationAxis, N_particles):
    '''Given particle angular momenta of the particles and the rotation axis
    of the galaxy, compute the angular momenta of the particles in the
    direction of galaxy rotation.'''
    Lperp = np.zeros(N_particles)
    for j in xrange(N_particles):
        Lperp[j] = np.dot(L[j], Galaxy_RotationAxis)

    return Lperp

def compute_EnclosedMass(gas, dm, stars, bh, GalaxyCentre):
    '''Takes all particle data for a galaxy, and the coordinates of the centre
    and for the radius of every particles, computes the enclosed mass within that
    radius. Returns arrays of radii and enclosed mass'''

    # combining datasets
    combined = {}
    combined['mass'] = np.concatenate((gas['mass'], dm['mass'], stars['mass'], bh['mass']))
    combined['coords'] = np.vstack((gas['coords'], dm['coords'], stars['coords'], bh['coords']))
    r_combined = np.linalg.norm(combined['coords'] - GalaxyCentre, axis=1)
    r_combined_vect = combined['coords']
    # print 'lengths comparison', len(r_combined), len(r_combined_vect)

    # sorting
    mask = np.argsort(r_combined)
    r_combined = r_combined[mask]
    m_combined = combined['mass'][mask]
    r_combined_vect = r_combined_vect[mask] #check for error

    # summing
    EnclosedMass = np.cumsum(m_combined)

    return EnclosedMass, r_combined, r_combined_vect, m_combined

def compute_L_c(r, m, G, EnclosedMass, r_combined):
    '''for a given particle in a galaxy, computes its L if it were on a circular
    orbit at a specified radius r'''
    v_c_1 = np.sqrt(G * EnclosedMass[np.searchsorted(r_combined, r) - 1] / r)
    L_c_1 = m * r * v_c_1

    return L_c_1

def sr_compute_All_L_c(r, m, G, EnclosedMass, r_combined, N_particles):
    '''runs compute_L_c for all star particles in a galaxy'''
    L_c = np.zeros(N_particles)
    v_c = np.zeros(N_particles)
    for j in xrange(N_particles):
        mod_r = np.linalg.norm(r[j])
        L_c[j] = compute_L_c(mod_r, m[j], G, EnclosedMass, r_combined)

    return L_c

def get_perp_vector(Galaxy_RotationAxis):
    perp_vector = np.cross(Galaxy_RotationAxis, [0, 0, 1])
    perp_vector = perp_vector / np.linalg.norm(perp_vector)  # Normalising
    return perp_vector

def initialise_interpolation(N_radii, N_All_Particles, Galaxy_RotationAxis, GC, EnclosedMass, m_combined, r_combined_vect, r_combined, G, eps=0.0001):
    perp_vector = get_perp_vector(Galaxy_RotationAxis)
    Pot_temp = np.zeros(N_All_Particles)
    InterpRadii = np.logspace(-4, 0, N_radii)
    InterpPotential = np.zeros(N_radii)
    InterpEnergyPUM = np.zeros(N_radii)
    eps = 0.0001 #replace with eagle paper value

    for j in xrange(N_radii):
        print 'radcount:', j
        r_temp = (InterpRadii[j] * perp_vector) + GC
        for k in xrange(N_All_Particles):
            # print 'j:',j,'kfrac:',k/N_All_Particles
            dist = np.sqrt(np.linalg.norm(r_combined_vect[k] - r_temp) ** 2 + eps ** 2)
            #print dist
            Pot_temp[k] = -(G * m_combined[k]) / (dist)
            # print 'Pot_temp',Pot_temp[k],'k:',k
        InterpPotential[j] = np.sum(Pot_temp)
        # print InterpPotential[j]
        InterpEnergyPUM[j] = (0.5 * G * EnclosedMass[np.searchsorted(r_combined, InterpRadii[j]) - 1] / InterpRadii[j]) + InterpPotential[j]
        # print InterpEnergyPUM[j]

    return InterpRadii, InterpPotential, InterpEnergyPUM

def initialise_interpolation_sample(N_radii, N_All_Particles, Galaxy_RotationAxis, GC, EnclosedMass, m_combined, r_combined_vect, r_combined, G, eps=0.0001, n_sample=10):
    perp_vector = get_perp_vector(Galaxy_RotationAxis)
    Pot_temp = np.zeros(N_All_Particles)
    InterpRadii = np.logspace(-4, 0, N_radii)
    InterpPotential = np.zeros(N_radii)
    InterpEnergyPUM = np.zeros(N_radii)
    eps = 0.0001

    for j in xrange(N_radii):
        print 'radcount:', j
        r_temp = (InterpRadii[j] * perp_vector) + GC
        for k in xrange(0, N_All_Particles, n_sample):
            # print 'j:',j,'kfrac:',k/N_All_Particles
            dist = np.sqrt(np.linalg.norm(r_combined_vect[k] - r_temp) ** 2 + eps ** 2)
            #print dist
            Pot_temp[k] = -(G * m_combined[k]) / (dist)
            # print 'Pot_temp',Pot_temp[k],'k:',k
        InterpPotential[j] = np.sum(Pot_temp)*n_sample
        # print InterpPotential[j]
        InterpEnergyPUM[j] = (0.5 * G * EnclosedMass[np.searchsorted(r_combined, InterpRadii[j]) - 1] / InterpRadii[j]) + InterpPotential[j]
        # print InterpEnergyPUM[j]

    return InterpRadii, InterpPotential, InterpEnergyPUM

def interp_compute_All_L_c(N_particles, r, m, v, InterpRadii, InterpPotential, InterpEnergyPUM, EnclosedMass, r_combined, G):
    L_c = np.zeros(N_particles)
    for j in xrange(N_particles):
        mod_r = np.linalg.norm(r[j])
        Pot_current = np.interp(mod_r, InterpRadii, InterpPotential)
        E_current = (0.5 * (np.linalg.norm(v[j]) ** 2)) + Pot_current
        r_new = np.interp(E_current, InterpEnergyPUM, InterpRadii)

        #print 'r_new:', r_new, 'Ecurr:', E_current, 'Eint:', InterpEnergyPUM
        L_c[j] = compute_L_c(r_new, m[j], G, EnclosedMass, r_combined)

    return L_c

def compute_circularities(simname, snapshot, GNs, SGNs, Centres, Velocities, OrbitFindingMethod=2, selection=0):
    """for a set of galaxies specified by group and subgroup numbers, computes
    the circularities of their particles as well as bulge mass particles for each
    and for the total, and exports these as a pickle."""
    n = len(GNs)
    #n=1
    AllBulgeIDs = []
    AllDiscIDs = []
    simsnap = simname + "/" + snapshot + "/"
    #plt.figure(figsize=(2, 2*n))
    for i in xrange(n):
        print 'Progress:', i / n
        # loading galaxy info
        try:
            stars = read_galaxy(4, GNs[i], SGNs[i], Centres[i], snapshot)
            bh = read_galaxy(5, GNs[i], SGNs[i], Centres[i], snapshot)
        except:
            "no stars/bh"
            continue
        gas = read_galaxy(0, GNs[i], SGNs[i], Centres[i], snapshot)
        dm = read_galaxy(1, GNs[i], SGNs[i], Centres[i], snapshot)

        # separate data
        r = stars['coords'] - Centres[i]

        # plt.figure(figsize=[5,5])
        # plt.plot(r[:,0], r[:,1], "b.", ms=0.5)
        # # pdb.set_trace()
        # plt.show()

        #bd_temp = stars['bd']
        pID = stars['ParticleIDs']
        v = stars['velocity'] - Velocities[i]
        m = stars['mass']
        #print np.sum(m), StellarMass[i]
        N_particles = len(m)

        # Finding axis of galaxy rotation
        L = compute_L(r, m, v, N_particles)
        Galaxy_RotationAxis = compute_RotationAxis(L)

        # Finding L around rotation axis
        Lperp = compute_Lperp(L, Galaxy_RotationAxis, N_particles)

        # computing M(<R) for all R
        EnclosedMass, r_combined, r_combined_vect, m_combined = compute_EnclosedMass(gas, dm,
                                                                            stars, bh, Centres[i])#galaxycentres
        # plt.scatter(np.log10(r_combined), np.log10(m_combined))

        N_All_Particles = len(r_combined)

        # Finding L of Circular orbit with the same energy
        # 0 - Aperture
        # 1 - Assume Same Radius
        # 2 - Linear interpolation
        # 3 - Log interpolation WIP
        # 4 - Hernquist profile WIP

        myG = const.G.to(u.Mpc ** 3 * u.Msun ** -1 * u.s ** -2).value

        if OrbitFindingMethod == 1:
            L_c = sr_compute_All_L_c(r, m, myG, EnclosedMass, r_combined, N_particles)

        if OrbitFindingMethod == 2:

            InterpRadii, InterpPotential, InterpEnergyPUM = initialise_interpolation_sample(10, N_All_Particles,
                                                    Galaxy_RotationAxis, Centres[i], EnclosedMass, m_combined,
                                                    r_combined_vect, r_combined, myG, eps=0.00001, n_sample=10)

            # plt.figure()
            # plt.plot(InterpRadii, InterpPotential, 'b.', label="nsample=10")
            # plt.xlabel("Log radius")
            # plt.ylabel("Log potential")
            #
            # InterpRadii, InterpPotential, InterpEnergyPUM = initialise_interpolation_sample(10, N_All_Particles,
            # Galaxy_RotationAxis, Centres[i], EnclosedMass, m_combined, r_combined_vect, r_combined, myG, eps=0.00001,
            # n_sample=1)

            # plt.plot(InterpRadii, InterpPotential, 'r.', label="nsample=1")
            # plt.legend()
            # plt.savefig("randSampComparison.png")

            L_c = interp_compute_All_L_c(N_particles, r, m, v, InterpRadii, InterpPotential, InterpEnergyPUM,
                                         EnclosedMass, r_combined, myG)

        Circularity = Lperp / L_c

        #plt.subplot(n, 1, i+1)
        # plt.hist(Circularity, bins=40, range=(-2,2))
        # plt.xlim(-2,2)
        # plt.show()

        print Circularity
        centrestring = str(Centres[i][0])+'-'+str(Centres[i][1])+'-'+str(Centres[i][2])
        filename = 'galdata/' + simsnap + str(GNs[i])[:-2] + '-' + str(SGNs[i])[:-2] + '-' + centrestring +'.pickle'
        circ_dict = dict(zip(pID, Circularity))
        pickle_out_circ = open(filename, 'wb')
        pickle.dump(circ_dict, pickle_out_circ)
        pickle_out_circ.close()

    return 0

def get_credentials():
    file = open('LoginDetails.txt')
    filelines = file.readlines()
    file.close()
    username = filelines[0].rstrip()
    password = filelines[1].rstrip()
    return username, password

def generate_circdata(simname, snapshot):
    # snapshot = "28"
    username, password = get_credentials()
    print username, 'string', password

    myQuery = 'SELECT \
                                SH.MassType_Star, \
                                SH.GalaxyID, \
                                SH.GroupNumber, \
                                SH.SubGroupNumber, \
                                SH.CentreOfPotential_x, \
                                SH.CentreOfPotential_y, \
                                SH.CentreOfPotential_z, \
                                SH.Velocity_x, \
                                SH.Velocity_y, \
                                SH.Velocity_z \
                            FROM \
                                {}_SubHalo as SH \
                            WHERE \
                                SH.SnapNum = {} \
                                and SH.MassType_star > 1e10\
            '.format(simname, snapshot)



    myData = read_galaxy_database(myQuery, 'SQLpickles/CircularitiesGeneratorSQL/{}_{}.pickle'.format(simname,snapshot),
                                  username, password, simname, data_request=1)

    print myData

    GroupNum, SubGroupNum, GalaxyCentres, GalaxyVelocities = extract_galaxydata(myData, snapshot)

    print "GalaxyCentres:", len(GroupNum)

    compute_circularities(simname, snapshot, GroupNum, SubGroupNum, GalaxyCentres, GalaxyVelocities, 2)

def main():
    simname = "RefL0012N0188"
    # simname = "RefL0025N0376"
    for i in [28]:
        snum = str(i)
        if i < 10:
            snum = "0" + snum
        print snum
        generate_circdata(simname, snum)


if __name__ == '__main__':
    main()
