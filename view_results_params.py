"""
Lauren R. Madden 2023
"""
import pickle
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import signal
import glob
import os
import subprocess

#############################################################################################################
######################### DISPLAY EXTRACELLULAR VOLTAGE RECORDING RESULTS ###################################
#############################################################################################################

# get most recent output file
os.chdir('Results')
allfiles = glob.glob(os.getcwd().replace('\\','/')+'/*')
latest_c_file = max(allfiles, key=os.path.getctime)
latest_m_file = max(allfiles, key=os.path.getmtime)
if latest_c_file != latest_m_file:
    latest_file = latest_m_file
else:
    latest_file = latest_c_file
filename = latest_file.replace("\\","/")

if filename.split('_')[-1] == 'membrane':
    variable = 'i_membrane'
elif filename.split('_')[-1] == 'v':
    variable = 'v'
elif filename.split('_')[-1] == 'withAEC':
    if filename.split('_')[-2] == 'membrane':
        variable = 'i_membrane_withAEC'
    else:
        variable = 'v_withAEC'

if 'i_membrane' in variable:
    output = open(filename,'rb')
    data = pickle.load(output)
    output.close()

    t = np.asarray(data['t'])
    wav_ds = np.asarray(data['waveform_ds'])
    secs = data['section_rows']
    sections = list()
    for sec in secs:
        sections.append(sec)

    if 1/(t[1] - t[0]) == 200:
        fs = 200e3
    else:
        fs = 500e3

    def downsample(sig, og_fs, samp_fs):
        samp_bin = range(0,len(sig),int(len(sig)/np.ceil(len(sig)*(samp_fs/og_fs))))
        new_sig = sig[samp_bin]
        return new_sig

    t_ds = downsample(t,fs,30e3)

    # Make Filter
    b,a = signal.butter(4,250,'highpass',output='ba',fs=30e3) # Make Filter

    # Filter Signal
    wav_ds_filt = signal.lfilter(b,a,wav_ds)

    plt.figure()
    plt.plot(t_ds,wav_ds_filt*1e3)
    plt.ylabel("uV")
    plt.xlabel("ms")
    simdetails = filename.split('/resultfile_')[-1].split('_i_membrane')[0]
    plt.title('Filtered Downsampled Waveform:\n'+simdetails)

else:
    output = open(filename,'rb')
    data = pickle.load(output)
    output.close()

    t = np.asarray(data['t'])
    v = np.asarray(data['v_array'])
    secs = data['section_rows']

    iseg_lastidx = 0
    for name in secs.keys():
        if 'iseg' in name:
            curidx = name.split('[')[1].split(']')[0]
            curidx = int(curidx)
            if curidx > iseg_lastidx:
                iseg_lastidx = curidx

    iseg_idx_A = int(iseg_lastidx*0.75)
    iseg_idx_B = int(iseg_lastidx*0.07)

    plt.plot(t,v[secs.get('nodeP[12]')])
    plt.plot(t,v[secs.get('nodeT[3]')])
    plt.plot(t,v[secs.get('nodeC[20]')])
    plt.plot(t,v[secs.get('iseg['+str(iseg_idx_A)+']')])
    plt.plot(t,v[secs.get('iseg['+str(iseg_idx_B)+']')])
    plt.plot(t,v[secs.get('soma[20]')])

    plt.ylabel("mV")
    plt.xlabel("ms")
    simdetails = filename.split('/resultfile_iSeg_')[-1].split('_v')[0]
    plt.title('Intracellular Waveform:\n'+simdetails)
    plt.legend(('distal peripheral node','stem node','central node','GIS compartment near stem','GIS compartment near soma','soma compartment'))

plt.xlim([9,19])
plt.savefig(filename+'_figure.pdf')
subprocess.Popen(filename+'_figure.pdf',shell=True)

