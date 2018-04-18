import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib

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

def read_dataset(itype, att, nfiles=16):
    """ Read a selected dataset, itype is the PartType and att is the attribute name. """

    # Output array.
    data = []

    # Loop over each file and extract the data.
    for i in range(nfiles):
        f = h5py.File('RefL0012N0188/snap_028_z000p000.%i.hdf5'%i, 'r')
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

# Get values
x = read_dataset(0, 'Density')
y = read_dataset(0, 'Temperature')
print x, y

# Create plot
plt.figure(figsize=(6, 4.5))
plt.hist2d(np.log10(x), np.log10(y), (500,500), norm=colors.LogNorm(), cmap='Greys')#, 'b.', ms=2)
plt.colorbar()
plt.axvline(x=-29, color='red', label='Critical Density')
plt.axvline(x=-29+np.log10(60*0.04), color='blue', label='Mean Galaxy Density')
plt.axvline(x=-29+np.log10(0.04), color='black', label='Mean Baryon Density')
plt.axvline(x=-25.5, color='green', label='Star Formation Minimum Density')
#plt.axvline(x=np.log10(6.8)-23, color='yellow', label='1 M$_{\odot}$ pc$^{-3}$')
plt.axvline(x=-24, color='teal', label='Mean Disk Formation Density')
plt.axvline(x=-21.3, color='orange', label='Mean Bulge Formation Density')
plt.legend(fontsize='x-small', labelspacing=0.2, borderpad=0.2)
plt.xlabel('Log$_{10}$[Gas Density /gcm$^{-3}$]')
plt.ylabel('Log$_{10}$[Gas Temperature /K]')
print 'plotted'
plt.tight_layout()
plt.savefig('TemperatureDensityGas.png')
#plt.show()