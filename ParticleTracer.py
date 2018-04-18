from __future__ import division
import os
from CircularitiesGenerator import *
import matplotlib
import random

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

def extract_snapshot_circdata(simname, snapshot):
    simsnap = simname + "/" + snapshot + "/"
    snapdata = {}

    for filename in os.listdir("galdata/"+simsnap):
        if filename[-7:] == ".pickle":
            f = open("galdata/"+simsnap+filename, 'rb')
            tempdata = pickle.load(f)
            f.close()
            snapdata.update(tempdata)

    return snapdata

def extract_all_circdata(simname):
    data = []

    for i in xrange(29):

        print "eac", i

        getdata = 1

        if getdata == 1:

            if i < 10:
                snapshot = "0" + str(i)
            else:
                snapshot = str(i)
            snapdata = extract_snapshot_circdata(simname, snapshot)

            # print "begin pickle {}".format(i)
            # pickle_out = open("ParticleTracer{}.pickle".format(i), 'wb')
            # pickle.dump(snapdata, pickle_out)
            # pickle_out.close()
            # print "end pickle"

        else:
            # print "begin unpickle {}".format(i)
            # f = open("ParticleTracer{}.pickle".format(i), 'rb')
            # snapdata = pickle.load(f)
            # f.close()
            # print "end unpickle"
            pass
        data.append(snapdata)

    return data

def get_bulgeIDs(simname):
    print "GbIDs"
    snapshot = "28"
    simsnap = simname + "/" + snapshot + "/"
    bulgeIDs = []

    for filename in os.listdir("galdata10e10/" + simsnap):
        if filename[-9:] == "29.pickle":
            f = open("galdata/" + simsnap + filename, 'rb')
            tempdata = pickle.load(f)
            f.close()
            for key, value in tempdata.iteritems():
                if value > -2:
                    bulgeIDs.append(key)

    return bulgeIDs

if __name__ == '__main__':
    # basedir = '/opt/local/l4astro/tt/L0025N0375/REFERENCE/snapshot_028_z000p000/'
    basedir = "RefL0012N0188/"
    # simname = "RefL0025N0376"
    simname = "RefL0012N0188"

    data = extract_all_circdata(simname)
    # print data

    bulgeIDs = get_bulgeIDs(simname)

    plt.figure(figsize=(6, 4.5))

    random.shuffle(bulgeIDs)

    bulgeIDs = bulgeIDs[0:1000]

    bidlength = len(bulgeIDs)
    # print bidlength
    k = 1

    end_circs = []
    start_circs = []

    for id in bulgeIDs:
        print k/bidlength
        snaps = []#np.zeros(29)
        particle_circs = []#np.zeros(29)
        k+=1
        for i in xrange(29):
            if id in data[i]:
                particle_circs.append(data[i][id])#[i] = data[i][id]
                #print particle_circs[i]
                snaps.append(i)#[i] = i

        plt.plot(snaps, particle_circs, 'r-', alpha = 0.05)
        start_circs.append(particle_circs[0])
        end_circs.append(particle_circs[-1])
        #plt.show()

    plt.ylim(-1.1, 1.1)
    plt.xlabel("Snapshot")
    plt.ylabel("Circularity")
    plt.tight_layout()
    plt.savefig("ParticleTracer.png")

    plt.figure(figsize=(6,5))
    plt.plot(start_circs, end_circs, 'b.', ms=0.5)
    plt.plot([-1, 1], [-1,1], 'r-')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.xlabel("Birth Circularity")
    plt.ylabel("Present Circularity")
    plt.tight_layout()
    plt.savefig("ParticleTracer2.png")

    #print data


