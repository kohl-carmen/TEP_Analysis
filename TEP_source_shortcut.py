# -*- coding: utf-8 -*-
"""

@author: ckohl
"""

## TEP Source Analysis - Single Pulse TMS - 
# takes data from TEP_preproc_single.m
# brief version


import mne
import numpy as np
from mayavi import mlab  
import matplotlib.pyplot as plt
from mne.forward import read_forward_solution
from mne.minimum_norm import (make_inverse_operator, apply_inverse,write_inverse_operator)
from surfer import Brain
from pptx import Presentation
import seaborn as sns
    print('indent something so that, if you doughnut press F5, it crashes out before doing any damage')


# =============================================================================  
# SET DIRS
# =============================================================================  

eeg_partic='Beta02'
eeg_dir='C:\\Users\\ckohl\\Desktop\\Current\\Data\\TMS\\EEG\\fieldtrip\\'
mri_partic= eeg_partic
mri_dir='C:\\Users\\ckohl\\Desktop\\Current\\Data\\Virtual_Shared\\Pilot\\'+mri_partic+'\\MRI'

eegfile=eeg_dir+'\\'+eeg_partic+'_fieldtrip1.set'
ppt=True
# =============================================================================  
## GET DATA
# =============================================================================
     
epochs=mne.read_epochs(fname=eeg_dir+eeg_partic+'_MNE_fieldtrip_epochs.fif',preload=True)
evoked= epochs.average()
dt=1

cov = mne.read_cov(mri_dir+'\\'+mri_partic+'\\single_cov.fif')
fwd=mne.read_forward_solution(mri_dir+'\\'+mri_partic+'\\single-fwd.fif')
src=mne.read_source_spaces(fname=mri_dir+'\\'+mri_partic+'\\single-src.fif')  
bem=mne.read_bem_solution(fname=mri_dir+'\\'+mri_partic+'\\single-bem.fif')  
trans= mne.read_trans(fname=mri_dir+'\\'+mri_partic+'-single-trans.fif')

# =============================================================================  
# INVERSESOLUTION
# =============================================================================  
inv = make_inverse_operator(epochs.info, fwd, cov, loose=0.001) 
#loose is super small because we essentially want the vertices to be normal to the cortical sueface, but we odn't wuite want to force it
snr = 3. 
lambda2 = 1. / snr ** 2 
# dSPM
stc_dspm = apply_inverse(evoked, inv, lambda2, method='dSPM', pick_ori="vector", return_residual=False, verbose=True) 
# MNE
stc_mne = apply_inverse(evoked, inv, lambda2, method='MNE', pick_ori="vector", return_residual=False, verbose=True) 

brain=stc_mne.plot(hemi='lh', subjects_dir=mri_dir,initial_time=.1)


#------------------
# If you already know where to look

# Find max in dSPM
# dSPM gives unitless outcomes, like z-scores. tha'ts great for establishing peaks, but not what we want in the end, so I define peaks based on dSPM and look at hose peak vertices in mne
tmin=40
tmax=60
t_oi_i = list(range(np.where(np.round(stc_dspm.times,3)==tmin/1000)[0][0], np.where(np.round(stc_dspm.times,3)==tmax/1000)[0][0]+1)) #time indeces of interes
data_dspm_normal = stc_dspm.normal(inv['src']).data 
# #this is Ryan's.this seems to be a very fancy way of asking fro the maximum. (fancy because these guys are vertices)
# first, we take only the time_oi sources, which gives a sources x time matrix. then, we take the vector norm (=vector magnitude) over each row (axis1, I think?) (https://www.educative.io/edpresso/what-is-the-nplinalgnorm-method-in-numpy)
# so that would give us one sort of magnitude per source (the way I understand this is that it's not really giving us the maximum per source, but the total activity - so I think the max source we get in the end of this might not be the absolute max max, but what contributed the most over a period of time)
# then we jsut take the max (spits out index of the max, not max value)
peak_vert = np.argmax(np.linalg.norm(data_dspm_normal[:, t_oi_i], axis=1)) #Vertex that contributes the most over window of time 
# so we found at which vertex was the most activity, now we want to find the max time spot in that vertex
peak_time = stc_dspm.times[t_oi_i[np.argmax(abs(data_dspm_normal[peak_vert, t_oi_i]))]] #Time at which abs(peak amplitude) of the selected vertex occurs 

#plot timecourse of this vertex in MNE
print('Peak vertex occurs at t =',peak_time) 
plt.ion()
fig, ax = plt.subplots()
plt.plot(stc_mne.times*1e3 ,  stc_mne.normal(inv['src']).data[peak_vert, :].T*(1e9), 'k-') 
plt.xlabel('Time (ms)') 
plt.ylabel('Amplitude (nAm)') 
plt.xlim((-100, 300)) 
ylims=ax.get_ylim()
plt.ylim(ylims)
plt.plot([0,0],ylims,':',color='grey')

# Find ID of selected vertex 
if peak_vert >= len(stc_dspm.vertices[0]): #which hemi
    pick_vertID = stc_dspm.vertices[1][peak_vert-len(stc_dspm.vertices[0])] 
else: 
    pick_vertID = stc_dspm.vertices[0][peak_vert]

# Set the components of all vertices besides the selected vertex to zero and ifnalte peak_oi
# (illsutration only obvs)
stc_thresh = stc_mne.copy() 
subthresh_inds = [i for i in range(stc_thresh.data.shape[0]) if i!=peak_vert] #indices of all sources besides peak
stc_thresh.data[subthresh_inds,:,:] = 1e-15 
stc_thresh.data[peak_vert,:,:] = stc_thresh.data[peak_vert,:,:] * 10 

surfer_kwargs = dict( 
    subjects_dir=mri_dir,  
    subject=mri_partic,  
    # clim=dict(kind='value', lims=[0,np.percentile(stc_mne.magnitude().data,95)*1e9,np.max(stc_mne.magnitude().data)*1e9]),  
    clim=dict(kind='value', lims=[np.percentile(stc_mne.magnitude().data,70),np.percentile(stc_mne.magnitude().data,95),np.max(stc_mne.magnitude().data)]),  
    initial_time= peak_time,   
    background='k', 
    foreground='w') 

brain = stc_thresh.plot(**surfer_kwargs) 
brain.add_foci(pick_vertID, coords_as_verts=True, hemi='lh', color='w', scale_factor=0.5, alpha=0.6) 


#------------------
# To explore time intervals
#plot basic brain
surfer_kwargs = dict( 
    clim=dict(kind='value', lims=[0, 1, 2]),  
    hemi='lh',
    subjects_dir=mri_dir,  
    subject=mri_partic,  
    initial_time= peak_time,   
    background='k', 
    foreground='k')
    # views=['lat','dors', 'med']) 
brain = stc_mne.plot(**surfer_kwargs) 

#put maxima on top
step=20#ms timestep
min=0
max=260

palette = sns.cubehelix_palette(light=.8, n_colors=int(max/step))
palette = sns.color_palette('YlOrRd', n_colors=int(max/step)+1)
plt.ion()
fig, ax = plt.subplots(3,int(np.round(max/step/3)+1))
for l in range(0,int(max/step)+1):
    t_oi_i = list(range(np.where(np.round(stc_dspm.times,3)==(min+l*step)/1000)[0][0], np.where(np.round(stc_dspm.times,3)==((min+l*step)+step)/1000)[0][0]+1)) #time indeces of interes
    data_dspm_normal = stc_dspm.normal(inv['src']).data 
    peak_vert = np.argmax(np.linalg.norm(data_dspm_normal[:, t_oi_i], axis=1)) #Vertex that contributes the most over window of time 
    peak_time = stc_dspm.times[t_oi_i[np.argmax(abs(data_dspm_normal[peak_vert, t_oi_i]))]] #Time at which abs(peak amplitude) of the selected vertex occurs 

    # Find ID of selected vertex 
    if peak_vert >= len(stc_dspm.vertices[0]): #which hemi
        pick_vertID = stc_dspm.vertices[1][peak_vert-len(stc_dspm.vertices[0])] 
    else: 
        pick_vertID = stc_dspm.vertices[0][peak_vert]

    #plot on brain
    brain.add_foci(pick_vertID, coords_as_verts=True, hemi='lh', color=palette[l], scale_factor=0.5) 

    #plot timecourse of this vertex in MNE
    fig.axes[l].plot(stc_mne.times*1e3 ,  stc_mne.normal(inv['src']).data[peak_vert, :].T*(1e9), '-',color=palette[l]) 
    fig.axes[l].set_title(str((min+l*step))+'-'+str(((min+l*step)+step))+' ms',fontsize='medium')
    fig.axes[l].set_xlim((-50, 300)) 
    ylim=[-0.03,0.03]
    fig.axes[l].set_ylim(ylim) 
    fig.axes[l].plot([0,0],ylim,':',color='grey')   
    fig.axes[l].axis('off')
    fig.axes[l].set_xlabel('time [s]')

    
    if l==int(max/step):
        for k in range(int(max/step),(3*int(np.round(max/step/3)+1))):
            fig.axes[k].axis('off')
# plt.tight_layout()

brain.animate(["m"] * 2)












# ---------------------------------
## ECD
from nilearn.plotting import plot_anat
from mne.evoked import combine_evoked
from mne.simulation import simulate_evoked
from mne.forward import make_forward_dipole

# evoked=evoked_full
evoked_full = evoked.copy()
evoked.crop(0.14, 0.16)
# Fit a dipole
dip = mne.fit_dipole(evoked, cov, bem, trans)[0]
# Plot the result in 3D brain with the MRI image.
dip.plot_locations(trans, mri_partic, subjects_dir=mri_dir, mode='orthoview')
#plot on scan
mni_pos = mne.head_to_mni(dip.pos, mri_head_t=trans,subject=mri_partic, subjects_dir=mri_dir)
mri_pos = mne.head_to_mri(dip.pos, mri_head_t=trans,subject=mri_partic, subjects_dir=mri_dir)

t1_fname = mri_dir+'\\'+ mri_partic+ '\\mri\\T1.mgz'
fig_T1 = plot_anat(t1_fname, cut_coords=mri_pos[0], title='Dipole loc.')
# #plot on standard
# from nilearn.datasets import load_mni152_template
# template = load_mni152_template()
# fig_template = plot_anat(template, cut_coords=mni_pos[0],title='Dipole loc. (MNI Space)')

#plot fied predicted by dipole with max goodness of fit, compare to data and take diff
fwd_dip, stc_dip = make_forward_dipole(dip, bem, evoked.info, trans)
pred_evoked = simulate_evoked(fwd_dip, stc_dip, evoked.info, cov=None, nave=np.inf)

# find time point with highest goodness of fit (gof)
best_idx = np.argmax(dip.gof)
best_time = dip.times[best_idx]
print('Highest GOF %0.1f%% at t=%0.1f ms with confidence volume %0.1f cm^3'
      % (dip.gof[best_idx], best_time * 1000,
         dip.conf['vol'][best_idx] * 100 ** 3))

# plot
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=[10., 3.4],
                         gridspec_kw=dict(width_ratios=[1, 1, 1],
                                          top=0.85))
vmin, vmax = -400, 400  # make sure each plot has same colour range

# first plot the topography at the time of the best fitting (single) dipole
plot_params = dict(times=best_time, ch_type='eeg', outlines='skirt',
                   colorbar=False, time_unit='s')
evoked.plot_topomap(time_format='Measured field', axes=axes[0], **plot_params)

# compare this to the predicted field
pred_evoked.plot_topomap(time_format='Predicted field', axes=axes[1],
                         **plot_params)

# Subtract predicted from measured data (apply equal weights)
diff = combine_evoked([evoked, pred_evoked], weights=[1, -1])
diff.plot_topomap(time_format='Difference', axes=axes[2:], **plot_params)
fig.suptitle('Comparison of measured and predicted fields '
             'at {:.0f} ms'.format(best_time * 1000.), fontsize=16)
fig.tight_layout()

# time course
# Estimate the time course of a single dipole with fixed position and orientation (the one that maximized GOF) over the entire interval
dip_fixed = mne.fit_dipole(evoked_full, cov, bem, trans,
                           pos=dip.pos[best_idx], ori=dip.ori[best_idx])[0]
dip_fixed.plot(time_unit='s')


#loop through diff time intervals
step=20#ms timestep
min=0
max=260
for l in range(0,int(max/step)):
    evoked= epochs.average()
    evoked.crop((min+l*step)/1000,((min+l*step)+step)/1000)
    # Fit a dipole
    dip = mne.fit_dipole(evoked, cov, bem, trans)[0]
    #plot on scan
    mni_pos = mne.head_to_mni(dip.pos, mri_head_t=trans,subject=mri_partic, subjects_dir=mri_dir)
    mri_pos = mne.head_to_mri(dip.pos, mri_head_t=trans,subject=mri_partic, subjects_dir=mri_dir)
    plt.ion()
    t1_fname = mri_dir+'\\'+ mri_partic+ '\\mri\\T1.mgz'
    fig_T1 = plot_anat(t1_fname, cut_coords=mri_pos[0], title=str((min+l*step))+'-'+str(((min+l*step)+step))+' ms')
    # #


