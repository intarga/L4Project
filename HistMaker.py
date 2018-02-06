from BulgeDiskSeparator4 import *

if __name__ == '__main__':
    simsnap = 'RefL0012N0188snap28'

    pickle_in_bulge = open('galdata/' + simsnap + '-all-bulge.pickle', 'rb')
    BulgeIDs = pickle.load(pickle_in_bulge)
    pickle_in_bulge.close()

    pickle_in_disc = open('galdata/' + simsnap + 'all-disc.pickle', 'rb')
    DiscIDs = pickle.load(pickle_in_disc)
    pickle_in_disc.close()

    pIDs = read_dataset(4, 'ParticleIDs')

    maskB = pIDs == BulgeIDs
    maskD = pIDs == DiscIDs

    birth_densities_Bulge = read_dataset(4, 'BirthDensity')[maskB]
    birth_densities_Disk = read_dataset(4, 'BirthDensity')[maskD]
    birth_densities_All = np.concatenate((birth_densities_Bulge, birth_densities_Disk))

    #plt.figure()
    plt.hist(np.log10(birth_densities_All), bins=50, histtype='step', color='black', normed=1)
    plt.xlabel('Log10[Birth Gas Density /gcm^-3]')
    plt.ylabel('Number Density of Particles')
    # print 'plotted'
    plt.savefig('HistMakerBD.png')