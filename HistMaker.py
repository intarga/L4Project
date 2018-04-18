from __future__ import division
from CircularitiesGenerator import *
import os
from matplotlib import colors
import matplotlib
import matplotlib.pyplot as plt
# import pickle

del matplotlib.font_manager.weight_dict['roman']
matplotlib.font_manager._rebuild()

plt.rcParams['text.latex.preamble']=[r"\usepackage{times}"]
#Options
params = {#'text.usetex' : True,
          'font.size' : 12,
          'font.family' : 'Times New Roman',
          #'text.latex.unicode': True,
          }
plt.rcParams.update(params)
#
# fig = plt.figure()

if __name__ == '__main__':
    # basedir = '/opt/local/l4astro/tt/L0025N0375/REFERENCE/snapshot_028_z000p000/'
    basedir = "RefL0012N0188/"
    simname = "RefL0012N0188"
    # simname = "RefL0025N0376"
    snapshot = "28"
    simsnap = simname + "/" + snapshot + "/"
    atts = ["Coordinates", "Mass", "BirthDensity", 'ElementAbundance/Oxygen', "ElementAbundance/Magnesium", 'ElementAbundance/Iron',
            "ElementAbundance/Hydrogen", 'StellarFormationTime', "Feedback_EnergyFraction", "ParticleIDs", "Circularity"]

    bulgeDataAll = {}
    discDataAll = {}
    allDataAll = {}
    for att in atts:
        bulgeDataAll[att] = []
        discDataAll[att] = []
        allDataAll[att] = []

    radius = []
    radcirc = []

    for filename in os.listdir("galdata10e10/"+simsnap):
        if filename[-7:] == ".pickle":
            f = open("galdata/"+simsnap+filename, 'rb')
            circularities = pickle.load(f)
            f.close()

            print circularities
            #circularities = circularities.values()
            circlist = []

            split = filename[0:-7].split('-')
            gn = float(split[0])
            sgn = float(split[1])
            centre = [float(split[2]), float(split[3]), float(split[4])]

            # print centre

            stars = read_galaxy(4, gn, sgn, centre, 28, atts[:-1])
            stars["Coordinates"] = stars["Coordinates"] * u.cm.to(u.Mpc)
            stars["Coordinates"] = stars["Coordinates"] - centre

            radius_temp = np.zeros(len(stars["Coordinates"]))
            for i in xrange(len(radius_temp)):
                radius_temp[i] = np.linalg.norm(stars["Coordinates"][i])

            print len(stars["BirthDensity"]), len(stars["Coordinates"])

            bulgeData = {}
            discData = {}
            allData = {}
            for att in atts:
                bulgeData[att] = []
                discData[att] = []
                allData[att] = []

            for i in xrange(len(circularities)):
                circlist.append(circularities[stars["ParticleIDs"][i]])
                if circularities[stars["ParticleIDs"][i]] < 0:
                    for att in atts:
                        if att == "Circularity":
                            bulgeData[att].append(circularities[stars["ParticleIDs"][i]])
                        else:
                            bulgeData[att].append(stars[att][i])

                elif circularities[stars["ParticleIDs"][i]] > 0.85:
                    for att in atts:
                        if att == "Circularity":
                            discData[att].append(circularities[stars["ParticleIDs"][i]])
                        else:
                            discData[att].append(stars[att][i])

            bulgeMassFrac = 2*np.sum(bulgeData["Mass"])/np.sum(stars["Mass"])
            print "Bulge Mass Fraction:", bulgeMassFrac
            print "Stellar Mass:", np.sum(stars["Mass"]) * u.g.to(u.Msun)/1E10
            bh = read_galaxy(5, gn, sgn, centre, 28, "BH_Mass")
            print "Black Hole Mass:", np.sum(bh["BH_Mass"]) * u.g.to(u.Msun)/1E10
            if bulgeMassFrac < 0.5:
                for att in atts:
                    allData[att] = bulgeData[att] + discData[att]
                    bulgeDataAll[att] += bulgeData[att]
                    discDataAll[att] += discData[att]
                    allDataAll[att] += allData[att]
                    radius += list(radius_temp)
                    radcirc += circlist

            print filename, "done"
        else:
            print filename, "not accepted"

    plt.figure(figsize=(3.15, 2.8))
    plt.hist(np.log10(allDataAll["BirthDensity"]), bins=50, histtype='step', color='black', normed=1)
    plt.hist(np.log10(bulgeDataAll["BirthDensity"]), bins=50, histtype='step', color='red', normed=1)
    plt.hist(np.log10(discDataAll["BirthDensity"]), bins=50, histtype='step', color='blue', normed=1)
    plt.xlabel('Log$_{10}$[Birth Gas Density /gcm$^{-3}$]')
    plt.ylabel('Number Density of Particles')
    # print 'plotted'
    plt.tight_layout()
    plt.savefig('HistMakerBD.png')

    plt.figure(figsize=(3.15, 2.8))
    plt.hist(allDataAll["StellarFormationTime"], bins=50, histtype='step', color='black', normed=1)
    plt.hist(bulgeDataAll["StellarFormationTime"], bins=50, histtype='step', color='red', normed=1)
    plt.hist(discDataAll["StellarFormationTime"], bins=50, histtype='step', color='blue', normed=1)
    plt.xlabel('Birth Expansion Factor')
    plt.ylabel('Number Density of Particles')
    # print 'plotted'
    plt.tight_layout()
    plt.savefig('HistMakerBEF.png')

    plt.figure(figsize=(3.15, 2.8))
    plt.hist(allDataAll["Feedback_EnergyFraction"], bins=50, histtype='step', color='black', normed=1)
    plt.hist(bulgeDataAll["Feedback_EnergyFraction"], bins=50, histtype='step', color='red', normed=1)
    plt.hist(discDataAll["Feedback_EnergyFraction"], bins=50, histtype='step', color='blue', normed=1)
    plt.xlabel('Feedback Energy Fraction')
    plt.ylabel('Number Density of Particles')
    # print 'plotted'
    plt.tight_layout()
    plt.savefig('HistMakerFEF.png')

    bulgeOFrac = np.log10(np.asarray(bulgeDataAll["ElementAbundance/Magnesium"])
                          / np.asarray(bulgeDataAll["ElementAbundance/Iron"]))
    discOFrac = np.log10(np.asarray(discDataAll["ElementAbundance/Magnesium"])
                         / np.asarray(discDataAll["ElementAbundance/Iron"]))
    allOFrac = np.log10(np.asarray(allDataAll["ElementAbundance/Magnesium"])
                        / np.asarray(allDataAll["ElementAbundance/Iron"]))

    bulgeHFrac = np.log10(np.asarray(bulgeDataAll["ElementAbundance/Iron"])
                          / np.asarray(bulgeDataAll["ElementAbundance/Hydrogen"]))
    discHFrac = np.log10(np.asarray(discDataAll["ElementAbundance/Iron"])
                         / np.asarray(discDataAll["ElementAbundance/Hydrogen"]))
    allHFrac = np.log10(np.asarray(allDataAll["ElementAbundance/Iron"])
                        / np.asarray(allDataAll["ElementAbundance/Hydrogen"]))

    # bulgeOFrac = bulgeOFrac + np.log10(55.85/16) - np.log10(4.9/0.282)
    # discOFrac = discOFrac + np.log10(55.85/16) - np.log10(4.9/0.282)
    # allOFrac = allOFrac + np.log10(55.85/16) - np.log10(4.9/0.282)

    bulgeOFrac = bulgeOFrac + np.log10(55.85 / 24.3) - np.log10(0.347 / 0.282)
    discOFrac = discOFrac + np.log10(55.85 / 24.3) - np.log10(0.347 / 0.282)
    allOFrac = allOFrac + np.log10(55.85 / 24.3) - np.log10(0.347 / 0.282)

    bulgeHFrac = bulgeHFrac + np.log10(1.008/55.85) - np.log10(2.82E-5)
    discHFrac = discHFrac + np.log10(1.008/55.85) - np.log10(2.82E-5)
    allHFrac = allHFrac + np.log10(1.008/55.85) - np.log10(2.82E-5)

    # bulgeOFrac = bulgeOFrac[np.logical_not(np.isnan(bulgeOFrac))]
    # discOFrac = discOFrac[np.logical_not(np.isnan(discOFrac))]
    # allOFrac = allOFrac[np.logical_not(np.isnan(allOFrac))]
    # bulgeOFrac = bulgeOFrac[bulgeOFrac < 2]
    # discOFrac = discOFrac[discOFrac < 2]
    # allOFrac = allOFrac[allOFrac < 2]
    #
    # bulgeHFrac = bulgeHFrac[np.logical_not(np.isnan(bulgeHFrac))]
    # discHFrac = discHFrac[np.logical_not(np.isnan(discHFrac))]
    # allHFrac = allHFrac[np.logical_not(np.isnan(allHFrac))]
    # bulgeHFrac = bulgeHFrac[bulgeHFrac < 2]
    # discHFrac = discHFrac[discHFrac < 2]
    # allHFrac = allHFrac[allHFrac < 2]

    mask = np.logical_not(np.isnan(allOFrac))
    allOFrac = allOFrac[mask]
    allHFrac = allHFrac[mask]

    mask = np.logical_not(np.isnan(allHFrac))
    allOFrac = allOFrac[mask]
    allHFrac = allHFrac[mask]

    mask = np.logical_not(np.isnan(bulgeOFrac))
    bulgeOFrac = bulgeOFrac[mask]
    bulgeHFrac = bulgeHFrac[mask]

    mask = np.logical_not(np.isnan(bulgeHFrac))
    bulgeOFrac = bulgeOFrac[mask]
    bulgeHFrac = bulgeHFrac[mask]

    mask = np.logical_not(np.isnan(discOFrac))
    discOFrac = discOFrac[mask]
    discHFrac = discHFrac[mask]

    mask = np.logical_not(np.isnan(discHFrac))
    discOFrac = discOFrac[mask]
    discHFrac = discHFrac[mask]

    # allMask = np.ones(len(allOFrac), dtype=bool)
    # allONan = np.isnan(allOFrac)
    # allHNan = np.isnan(allHFrac)
    # for i in allMask:
    #     if np.any(allONan[i], allHNan[i]):
    #         allMask[i] = False
    # allOFrac = allOFrac[allMask]
    # allHFrac = allHFrac[allMask]

    # print np.sort(allOFrac)

    plt.figure(figsize=(3.15, 2.8))
    plt.hist(allOFrac, bins=50, histtype='step', color='black', normed=1, range=(-1, 2))
    plt.hist(bulgeOFrac, bins=50, histtype='step', color='red', normed=1, range=(-1, 2))
    plt.hist(discOFrac, bins=50, histtype='step', color='blue', normed=1, range=(-1, 2))
    # plt.xlabel('Log10[Birth Gas Density /gcm^-3]')
    plt.ylabel('Number Density of Particles')
    plt.xlabel("[Mg/Fe]")
    # print 'plotted'
    plt.tight_layout()
    plt.savefig('HistMakerAlpha.png')

    plt.figure(figsize=(6,4))
    plt.subplot(211)
    #plt.hist2d(bulgeHFrac, bulgeOFrac, (100,100), norm=colors.LogNorm(), cmap='Greys', range=[[-3, 1], [-0.5, 1]])
    plt.plot(bulgeHFrac, bulgeOFrac, 'r.', ms=0.25)
    plt.xlim(-3, 1)
    plt.ylim(0.5, 2)
    plt.subplot(212)
    #plt.hist2d(discHFrac, discOFrac, (100,100), norm=colors.LogNorm(), cmap='Greys', range=[[-3, 1], [-0.5, 1]])
    plt.plot(discHFrac, discOFrac, 'b.', ms=0.25)
    plt.xlim(-3, 1)
    plt.ylim(0.5, 2)
    plt.xlabel('[Fe/H]')
    plt.ylabel('[O/Fe]')
    # print 'plotted'
    plt.tight_layout()
    plt.savefig('HistMakerAlpha2dOxygen.png')

    print len(radius), len(allDataAll["Circularity"])
    plt.figure(figsize=(6, 4.5))
    plt.hist2d(radius, radcirc, (100, 100), [[0, 0.05], [-1, 1]], norm=colors.LogNorm(), cmap='Greys')
    plt.xlabel("Radius from galactic centre /Mpc")
    plt.ylabel("Circularity")
    plt.axhline(0, color="red")
    plt.axhline(0.85, color="blue")
    plt.axhline(-0.85, color="green")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig('HistMaker2DCR')
