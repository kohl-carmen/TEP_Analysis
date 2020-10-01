# Preproc TMS/EEG Data

import mne
import numpy as np
from mayavi import mlab  
import matplotlib.pyplot as plt
from mne.forward import read_forward_solution
from mne.minimum_norm import (make_inverse_operator, apply_inverse,write_inverse_operator)
from surfer import Brain

# Set Dirs
eeg_partic='Beta02'
eeg_dir='C:\\Users\\ckohl\\Desktop\\Current\\Data\\TMS\\EEG\\'
mri_partic= eeg_partic
mri_dir='C:\\Users\\ckohl\\Desktop\\Current\\Data\\Virtual_Shared\\Pilot\\'+mri_partic+'\\MRI'
raw_name='Raw - EEG-20200912T162305Z-001\\actiCHamp_Plus_BC-TMS_BETA02_20200911_EOG000032.vhdr'
elec_dir='C:\\Users\\ckohl\\Desktop\\Current\\Data\\TMS\\'
elec_file='Exported Electrodes Edited.txt'

# Options
plot_steps=False
cut_pulse=[-2/1000, 18/1000]

# Load
raw=mne.io.read_raw_brainvision(eeg_dir+raw_name,eog=['EOG_aux1_vert','EOG_aux2_horz'],preload=False)

# Montage
import sys
sys.path.append(function_dir)
from my_mne_functions import import_dig_electrodes
montage=import_dig_electrodes(elec_dir+elec_file ,raw,plot=plot_steps)
raw.set_montage(montage)

# Events
events = mne.read_annotations(eeg_dir+raw_name[0:-4]+'vmrk')
raw.set_annotations(events)
events, event_ids = mne.events_from_annotations(raw)


# Interpolate Pulse
data, times = raw[:]  
count_trials=0
for event in events:
    if event[2]==13:
        count_trials+=1
        pulse_start=int(event[0]+cut_pulse[0]*raw.info['sfreq'] )
        pulse_end=int(event[0]+cut_pulse[1]*raw.info['sfreq'] )
        for chan in range(0, len(data)):
            #interpolate (slope=dy/dx, then y=i+xb)
            dy=(data[chan][pulse_end]-data[chan][pulse_start])
            dx=times[pulse_end]-times[pulse_start]
            slope=dy/dx
            x=range(pulse_start,pulse_end)
            data[chan][x]=data[chan][pulse_start]+times[range(0,len(x))]*slope
raw.load_data()            
raw[:]=data
del data

# Downsample
raw,events=raw.resample(1000,events=events)

#Filter
raw.filter(1, 80, fir_design='firwin')#,picks='eeg'
            
# Interpolate bad channels
if plot_steps:
    raw.plot_psd(fmax=60)
    raw.plot()
raw.info['bads'].extend(['TP9', 'Cz','Iz','TP8','T8'])
raw.interpolate_bads()

# Reref
raw.set_eeg_reference('average', projection=True)  

# Reject Blinks
eog_events = mne.preprocessing.find_eog_events(raw)
onsets = eog_events[:, 0] / raw.info['sfreq'] - 0.2
durations = [0.5] * len(eog_events)
descriptions = ['bad blink'] * len(eog_events)
blink_annot = mne.Annotations(onsets, durations, descriptions,orig_time=raw.info['meas_date'])
raw.set_annotations(blink_annot)
eeg_picks = mne.pick_types(raw.info, eeg=True)
if plot_steps:
    raw.plot(events=eog_events, order=eeg_picks)

# Reject Noise
reject_criteria=dict(eeg=100e-6)#100microVolts

# Epoch
epochs = mne.Epochs(raw, events, event_id=13, tmin=-0.2, tmax=0.5, baseline=(-0.2, -0.005),preload=True,
                    reject=reject_criteria, reject_by_annotation=True)

# Evoked
evoked= epochs.average()

if plot_steps:
    epochs.plot_drop_log()
    evoked.plot(spatial_colors=True)

# Save
mne.Epochs.save(epochs, fname=eeg_dir+eeg_partic+'_TEP_basic_MNE.fif')    