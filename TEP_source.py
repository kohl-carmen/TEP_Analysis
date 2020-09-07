# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:48:55 2020

@author: ckohl
"""

## TEP Spurce Analysis
# based on preproc2020mne.py


import mne
import numpy as np
from mayavi import mlab  
import matplotlib.pyplot as plt
from mne.forward import read_forward_solution
from mne.minimum_norm import (make_inverse_operator, apply_inverse,write_inverse_operator)
from surfer import Brain

    print('indent something so that, if you doughnut press F5, it crashes out before doing any damage')


# =============================================================================  
# =============================================================================
# SET DIRS - Manual
# =============================================================================
# =============================================================================  
# set partic
skip_preproc=True    
skip_source_prep=True
eeg_partic='Beta04'
plot_steps_preproc=False
plot_steps_source=False
use_prep_coreg=True
raw_dir='F:\\Brown\\TMS\\Pilot\\'
mri_partic= eeg_partic
#raw_dir=raw_dir+'\\'+eeg_partic+'_TEP_got_MNE'
shared_dir='C:\\Users\\ckohl\\Documents\\Virtual_Shared\\Pilot\\'+mri_partic+'\\MRI'



# =============================================================================  
# =============================================================================
##-LOAD-
# =============================================================================  
# =============================================================================
first_time=False
if first_time==True:
    #load eeglab epochs
    epochs=mne.io.read_epochs_eeglab(raw_dir+'\\'+eeg_partic+'_TEP_1k_runICA_filt100.set')
    #avg ref                            
    epochs.set_eeg_reference('average', projection=True)
    # MONTAGE
    ## --------------------------------
    #(this is not actually from the same session but better than standard)
    #digitised montage(this is a whoel script and includes several manual steps)
    import sys
    sys.path.append('C:\\Users\\ckohl\\Documents\\MNE Scripts\\')
    from my_mne_functions import *
    montage=import_dig_electrodes('C:\\Users\\ckohl\\Documents\\Virtual_Shared\\Pilot'+'\\'+eeg_partic+'\\EEG\\'+'\\Exported Electrodes faked.txt',epochs)
    epochs.set_montage(montage)
    #save
    mne.Epochs.save(epochs, fname=raw_dir+'\\'+eeg_partic+'_TEP_epochs.fif')    
else       
    epochs=mne.read_epochs(fname=raw_dir+'\\'+eeg_partic+'_TEP_epochs.fif',preload=True)
    epochs.set_eeg_reference('average', projection=True)



# =============================================================================  
# =============================================================================
#    # RE-PREPROCESS
#    og_epochs=epochs.copy()
#    epochs.filter(l_freq=1., h_freq=None)
#    # ICA (EOG)
#    ## --------------------------------     
#     #filter to remove slow drifts before ICA ( did a filter already but only at .1, now i should to 1. but i don't want that inmy data, so i'm gonna apply the filter on a copy, do ICA on the copy, and then apply it to my data afterwards. mne recommended) 
#    # A high-pass filter with 1 Hz cutoff frequency is recommended. However, because filtering is a linear operation, the ICA solution found from the filtered signal can be applied to the unfiltered signal (see 2 for more information), so we’ll keep a copy of the unfiltered Raw object around so we can apply the ICA solution to it later.
#    #https://mne.tools/dev/auto_tutorials/preprocessing/plot_40_artifact_correction_ica.html#tut-artifact-ica
#    from mne.preprocessing import (ICA, create_eog_epochs, create_ecg_epochs,corrmap)
#    #filt_raw = raw.copy()
#    #filt_raw.load_data().filter(l_freq=1., h_freq=None)
#    ica = ICA(n_components=15, random_state=97)
#    ica.fit(epochs,reject_by_annotation=True)
#    #plot ICs
#    ica.plot_sources(epochs)
#    ica.plot_components()
#    ica.plot_properties(epochs, picks=[0])
#    
#    ##to reject ICs manually
#    evoked=epochs
#    ica.plot_overlay(evoked, exclude=[0], picks='eeg')
#    ica.plot_overlay(evoked, exclude=[0], picks='C3')
#    #beta02 exclude=[0, 5] 
#    ica.exclude = [0] 
#    
#    
#    #or to reject automatically
#    #-----------------------
#    ica.exclude=[]
#    eog_indices, eog_scores = ica.find_bads_eog(epochs,reject_by_annotation=True)
#    ica.exclude = eog_indices
#    
#    
#    #thi is where it actually gets excluded (before just marked)
#    ica.apply(epochs)   
#    
#    
#    # REJECTION
#    ## --------------------------------   
#    #spectral plot
#    epochs.plot_psd(spatial_colors=True,fmax=50)   
#    
#    #manual
#    fig=raw.plot(n_channels=len(epochs.ch_names), scalings=dict(eeg=20e-6))
#    fig.canvas.key_press_event('a')  #press a to see options. this lets you define bad sections and you can call them diff things (mutliple choice the one you want to use)
#    
#    #save what I rejected manually
#    raw.annotations.save( fname=raw_dir+'\\'+eeg_partic+'_annotations_after_manual_art_rejection.fif')            
#    
#      
#    # INTERPOLATE
#    ## --------------------------------    
#    #interpolate dropped channelsback
#    epochs.interpolate_bads(reset_bads=False, verbose=False)
#
#    
#  fig=epochs.plot(n_epochs=1,n_channels=len(epochs.ch_names), scalings=dict(eeg=20e-6))
#    fig.canvas.key_press_event('a') 
# 
# =============================================================================  
# =============================================================================
# ERPs
# =============================================================================
# =============================================================================  
    
evoked= epochs.average()
#plot
if plot_steps_preproc==True:
    evoked.plot(spatial_colors=True)
    
    evoked.plot(picks=['C3'],xlim=[-.1, .3])
    evoked.plot_joint(times=[.01, .03, .05])


    times = np.linspace(0, 0.06, 10)
    evoked.plot_topomap(times=times, colorbar=True)
    
    mne.viz.plot_compare_evokeds(evoked, picks='C3')
                  
    mne.viz.plot_compare_evokeds(evoked, picks='eeg', axes='topo')  
    
    
    maps = mne.make_field_map(evoked, trans,subject=mri_partic, subjects_dir=shared_dir)
    evoked.plot_field(maps, time=0.07)


## ================================
# COVARIANCE
## ================================  
cov = mne.compute_covariance(epochs, tmax=0., method='auto',rank=None)
#save
#mne.write_cov(fname=raw_dir+'\\'+eeg_partic+'_cov_preproc2020.fif',cov=cov)
if plot_steps_preproc==True:
    evoked.plot_white(cov, time_unit='s', verbose=False)


### if you wana try different methods:
#method_params = dict(diagonal_fixed=dict(eeg=0.01))
#noise_covs = mne.compute_covariance(epochs, tmin=None, tmax=0, method='auto',
#                                return_estimators=True, verbose=True, n_jobs=1,
#                                projs=None, rank=None,
#                                method_params=method_params)
#evoked.plot_white(noise_covs, time_unit='s')
#    
#mne.viz.plot_cov(noise_cov, raw.info)
    
    
    
if skip_source_prep ==False:
    # =============================================================================  
    # =============================================================================
    # FORWARD SOLUTION
    # =============================================================================
    # =============================================================================  

    ## FORWARD OPERATOR
    ###to get fwd, we need bem,trans, and source space
    
    ## ================================
    # BEM
    ## ================================ 
    model=mne.make_bem_model(mri_partic,conductivity=(0.3, 0.006, 0.3), subjects_dir=shared_dir,ico=4)   
    bem = mne.make_bem_solution(model)      #inner skull, outer kull,skin
    if plot_steps_source==True:
        mne.viz.plot_bem(subject=mri_partic, subjects_dir=shared_dir, brain_surfaces='white',
                     orientation='coronal');
    # mine looks weird but apparently it's fine?https://github.com/mne-tools/mne-python/issues/5874
    #save
    #mne.write_bem_solution(fname=shared_dir+'\\'+mri_partic+'\\2020-bem.fif',bem=bem)
                     
                 
                                 
                 
    ## ================================
    # COREGISTER
    ## ================================
    use_prep_coreg=True
    if use_prep_coreg==False:    
        mne.gui.coregistration(inst=raw_dir+'\\'+eeg_partic+'_raw_preprocessed_in_preproc2020.fif', subject=mri_partic, subjects_dir=shared_dir )
       
    trans_file=raw_dir+'\\'+mri_partic+'-trans_2020.fif'
    trans= mne.read_trans(raw_dir+'\\'+mri_partic+'-trans_2020.fif')
    if plot_steps_source==True:    
        fig = mne.viz.plot_alignment(epochs.info, trans, subject=mri_partic, dig=True,
                         eeg=['original', 'projected'], meg=[],
                         coord_frame='head', subjects_dir=shared_dir)
    

    
    ## ================================
    # SOURCE SPCAE
    ## ================================
    
    src=mne.setup_source_space(subject=mri_partic,subjects_dir=shared_dir,spacing='oct6',add_dist=False)
    #mne.SourceSpaces.save(src,shared_dir+'\\'+mri_partic+'\\2020-src.fif')
    #src=mne.read_source_spaces(fname=shared_dir'\\'+mri_partic++'\\2020-src.fif')    

    if plot_steps_source==True:   
        brain = Brain(hemi='lh', surf='orig', subjects_dir=shared_dir,subject_id=mri_partic)
        surf = brain.geo['lh']
        vertidx = np.where(src[0]['inuse'])[0]
        mlab.points3d(surf.x[vertidx], surf.y[vertidx],surf.z[vertidx], color=(1, 1, 0), scale_factor=1.5)
        
    
    


    ## ================================
    # FORWARD
    ## ================================
    #fwd = mne.make_forward_solution(raw_dir+'\\'+eeg_partic+'_raw_preprocessed_in_preproc2020.fif', trans=trans, src=src, bem=bem,
    #                                meg=False, # include MEG channels
    #                                eeg=True, # include EEG channels
    #                                mindist=5.0, # ignore sources <= 5mm from inner skull
    #                                n_jobs=1) # number of jobs to run in parallel   
    fwd = mne.make_forward_solution(raw.info, trans=trans_file, src=src, bem=bem)    
    
    #Convert to surface orientation for cortically constrained inverse modeling
    
    fwd = mne.convert_forward_solution(fwd, surf_ori=True)
    
    
    mne.write_forward_solution(shared_dir+'\\'+mri_partic+'\\2020-fwd.fif',fwd)
    #leadfield = fwd['sol']['data']
    #print("Leadfield size : %d sensors x %d dipoles" % leadfield.shape)
    
    
    
    # SENSITIVITY/FIELD MAP
    #--------------------------------
    eeg_map = mne.sensitivity_map(fwd, ch_type='eeg', mode='fixed')
    
    if plot_steps_source==True:  
        clim = dict(kind='percent', lims=(0.0, 50, 95), smoothing_steps=3)  # let's see single dipoles
        brain = eeg_map.plot(subject=mri_partic, time_label='EEG sensitivity', surface='inflated',
                              subjects_dir=shared_dir, clim=clim, smoothing_steps=8, alpha=0.9);
        view = 'lat'
        brain.show_view(view)
    
    
        # FIELD MAP
        #---------------
        # Visualizing field lines based on coregistration
        from mne import make_field_map
        
        #make_field_map?
        maps = mne.make_field_map(evoked, trans=trans, subject=mri_partic,
                              subjects_dir=shared_dir, n_jobs=1)
        # explore several points in time
        field_map = evoked.plot_field(maps, time=.07);
    
    
else:
    
    fwd=mne.read_forward_solution(shared_dir+'\\'+mri_partic+'\\2020-fwd.fif')
    src=mne.read_source_spaces(fname=shared_dir+'\\'+mri_partic+'\\2020-src.fif')  
    bem=mne.read_bem_solution(fname=shared_dir+'\\'+mri_partic+'\\2020-bem.fif')  
    trans= mne.read_trans('C:\\Users\\ckohl\\Documents\\Virtual_Shared\\Pilot'+'\\'+eeg_partic+'\\EEG\\'+eeg_partic+'-trans_2020.fif')




# =============================================================================  
# =============================================================================
# INVERSESOLUTION
# =============================================================================
# =============================================================================  






## MNE
# =============================================================================
# =============================================================================  
        
evoked.set_eeg_reference('average',projection=True)
inverse_operator = make_inverse_operator(evoked.info, fwd, cov,loose=0.2, depth=0.8, verbose=False)


method = "MNE"
snr = 3.
lambda2 = 1. / snr ** 2
#stc = apply_inverse(evoked, inverse_operator, lambda2, method=method)
stc=mne.minimum_norm.apply_inverse(evoked,inverse_operator, lambda2, method,pick_ori='normal')


##basic plot at 70ms
#brain = stc.plot(surface='inflated', hemi='lh', subjects_dir=shared_dir,initial_time=.1)
#
#
##brain.scale_data_colormap(fmin=8, fmid=12, fmax=15, transparent=True)
#brain.show_view('lateral');
#
#
#
#
## plottin from EEG_pilot_source.py
## =============================================================================
#
#
## animation
## =============================================================================
## time dilation stretches animation in time (the bigger the longer)
#
#fig= stc.plot(initial_time=0, hemi='split', views=['lat','dors', 'med'], subjects_dir=shared_dir,subject=mri_partic)
##fig = stc.plot(initial_time=0,subjects_dir=shared_dir,subject=mri_partic,hemi='lh',views=['lat'],colormap='mne')
#fig.save_movie('04.gif',time_dilation=60,tmin=0,tmax=.1)
#fig.close()
#
#
#
#
#
#
#
#
#
## using label
## =============================================================================
##predefined. all sources in postcentral gyrus
#
#parc='aparc.a2009s'
#labels = mne.read_labels_from_annot(subject=mri_partic,parc=parc,subjects_dir=shared_dir) ## DKTatlas40 or a2009s
#
#label_oi='G_postcentral-lh'
#SI_label = [label for label in labels if label.name == label_oi][0]
#
##plot source on brain
## =====================
##brain = Brain(subject_id=mri_partic,subjects_dir=shared_dir,surf='orig',hemi='lh', background='white', size=(800, 600))
##brain.add_annotation(parc)
##brain.add_label(SI_label, borders=False,color=[1, 0, 0])
#
###plot label_oi & sources
##normals = src[1]['nn'][src[1]['vertno']]  
##brain = Brain(hemi='lh', surf='inflated', subjects_dir=shared_dir,subject_id=mri_partic)
##surf = brain.geo['lh']    
##mlab.quiver3d(surf.x[src[0]['vertno']], surf.y[src[0]['vertno']],surf.z[src[0]['vertno']], normals[:, 0], normals[:, 1], normals[:, 2], color=(0.0, 0.0, 0.0)  ,scale_factor=1.5)
##brain.add_label(SI_label, borders=False,color=[1, 0, 0])
#
##show whole brain parcellation
##brain = Brain(subject_id=mri_partic,subjects_dir=shared_dir,surf='orig',hemi='both', background='white', size=(800, 600))
##brain.add_annotation(parc)
#
#mean_lh = stc.extract_label_time_course(SI_label, inverse_operator['src'], mode='mean')
#
### plot waveform
## =====================
#
### data
#data_x=stc.times
#data_y=mean_lh.T
##Am ton Am
#data_y=[num*1e+9 for num in data_y]
#
### get jones waveform
#jones=open('C:/Users/ckohl/Desktop/Current/jones_supra_extracted.txt','r')#ead()
#jones_x=[]
#jones_y=[]
#for line in jones.readlines():
#    sep=line.find('/t')
#    jones_x.append(float(line[0:sep]))
#    jones_y.append(float(line[sep+2:]))
##plt.plot(jones_x,jones_y)    
#
## fit x axis
#if max(data_x) <5: #data is in seconds
#    jones_x=[num /1000 for num in jones_x]
#
##fit jones scale
#import scipy.optimize
#
#def min_this(scale,args):
#    import numpy as np
#    abs_mag_data=max(data_y)-min(data_y)
#    abs_mag_jones=max(jones_y)-min(jones_y)
#    difference_to_minimise=np.abs(np.abs(abs_mag_jones)*scale-np.abs(abs_mag_data))
#    return(difference_to_minimise)                         
#scale=1;
#res = scipy.optimize.minimize(min_this, scale, args=[data_y,jones_y], method='nelder-mead') #,options={'maxiter': 10000, 'disp': True})
#jones_y_scaled=[num*res.x for num in jones_y]
#
##plot
#fig, ax = plt.subplots()
#plt.plot(data_x,data_y,'k')
#plt.plot(jones_x,jones_y_scaled,':',color='grey')
#ax.set_ylabel('Amplitude(nAm)')
#ax.set_xlabel('Time(s)')
#ax.set_title('Source Average')
#ax.legend(['Source','Jones'])
#textstr='Jones scaled by: '+ str(res.x)
#props = dict(boxstyle='round', facecolor='gray', alpha=0.9)
#ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
#
#
#ylims=ax.get_ylim()
#plt.xlim(-0.05,0.3)
#plt.ylim(ylims)
#plt.plot([0,0],ylims,':',color='grey')
#plt.plot([0.02,0.02],ylims,':',color='grey')
#plt.plot([0.04,0.04],ylims,':',color='grey')
#
#plt.show()
#
#
#
#
#
#
#
#
#
#
#
#
#
## find maximum activation
## =============================================================================
#
#
#
## using get_peak function
## =======================
##don't have to restrict time where to look for peak but can
#tmin=.05
#tmax=.1
#vertno_max, time_max = stc.get_peak(hemi='lh',tmin=tmin,tmax=tmax)
#surfer_kwargs = dict(
#    subjects_dir=shared_dir, views=['lat','dors'],
#    initial_time=time_max, time_unit='s', size=(800, 800), smoothing_steps=5,surface='orig')
#brain = stc.plot(**surfer_kwargs)
#brain.add_foci(vertno_max, coords_as_verts=True, hemi='lh', color='red',
#               scale_factor=1, alpha=0.9)
#brain.add_text(0.1, 0.9, 'Peak between '+str(tmin*1000)+' and '+str(tmax*1000)+' ms: '+str(round(time_max*1000))+' ms', 'title',font_size=14)
#
#
#
#
## using my own: max over time
## =======================
#
#
##just get left hemisphere:
#lh_data=stc.data[:len(stc.vertices[0])]; #right: stc.data[:len(stc.vertices[1])]
#
##find maximum activation in LH o0ver time and 
#max_per_time=[]
#max_per_time_i=[]
#for t in range(0,lh_data.shape[1]): #loop through time points by index
#    max_per_time.append(abs(lh_data[:,t]).max())
#    max_per_time_i.append(stc.vertices[0][np.argmax(abs(lh_data[:,t]))])
#
##plot max over time (to find at what point was most activated)
#fig, ax = plt.subplots()
#plt.plot(stc._times,max_per_time)  
#plt.plot(stc._times[np.argmax(max_per_time)],max_per_time[np.argmax(max_per_time)],'or')  
#ax.set_ylabel('Amplitude(Am)')
#ax.set_xlabel('Time(s)')
#ax.set_title('LH max Am over time')
#textstr='Overall Max at '+ str(round(stc._times[np.argmax(max_per_time)]*1000)) +' ms'
#props = dict(boxstyle='round', facecolor='red', alpha=0.9)
#ax.text(0.65, 0.1, textstr, transform=ax.transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
#plt.show()
#
##plot location of all max (for each time)
#
#surfer_kwargs = dict(subjects_dir=shared_dir, subject_id=mri_partic, surf='inflated',hemi='lh',views='lateral')
#brain=Brain(**surfer_kwargs)
## go through time and plot dots. don't plot before time 0 
#max_t=[]
#count=0
#for t in range(0,lh_data.shape[1]): 
#    if stc._times[t]>0 and stc.times[t]<.2:
#        if stc._times[t]<.17:
#            brain.add_foci(max_per_time_i[t], coords_as_verts=True, hemi='lh', color='orange',scale_factor=0.5, alpha=0.8)
#            count+=1
#        else:
#            brain.add_foci(max_per_time_i[t], coords_as_verts=True, hemi='lh', color='yellow',scale_factor=0.5, alpha=0.5)
#            count+=1
#        if stc._times[t]==stc._times[np.argmax(max_per_time)]:
#            max_t=t
#brain.add_foci(max_per_time_i[max_t], coords_as_verts=True, hemi='lh', color='red',scale_factor=1, alpha=0.8)
#brain.add_text(0.1, 0.9, 'Maximum at '+str(round(stc._times[max_t]*1000))+' ms', 'title',font_size=14)
#



# lets look at peak activation (only one max vertex)
# =======================
#set time window to find peak within certain times
tmin=.00
tmax=.1
vertno_max, time_max = stc.get_peak(hemi='lh',tmin=tmin,tmax=tmax)
#also get indices
vertno_max_i, time_max_i = stc.get_peak(hemi='lh',tmin=tmin,tmax=tmax,vert_as_index=True, time_as_index=True)

#show peak on brain
surfer_kwargs = dict(
    subjects_dir=shared_dir, views='lateral',
    initial_time=time_max, time_unit='s', size=(800, 800), smoothing_steps=5)
brain = stc.plot(**surfer_kwargs)
brain.add_foci(vertno_max, coords_as_verts=True, hemi='lh', color='black',
               scale_factor=1, alpha=0.9)
brain.add_text(0.1, 0.9, 'Peak between '+str(tmin*1000)+' and '+str(tmax*1000)+' ms: '+str(round(time_max*1000))+' ms', 'title',font_size=14)

#take activation at only max time and see which sources should be included in peak, based on a cutoff
lh_data=stc.data[:len(stc.vertices[0])]; #right: stc.data[len(stc.vertices[1]:)]
max_data=lh_data[:,time_max_i]

#max_data=abs(max_data)
#take absolute? i guess not because I would just canel stuff out?
# so we have to distinguish if the max is a peak or a trough
#if lh_data[vertno_max_i,time_max_i]<0:
#    #max is trough
#    midway_in_value=max_data.min()+abs(max_data.min()-max_data.max())/2
#    midway_in_value_from0=max_data.min()-max_data.min()/2
#    cutoff=max_data.max()#midway_in_value_from0
#    #now find all sources where its above cutoff
#    source_mask=[]
#    for i in max_data:
#        if i <cutoff:
#            source_mask.append(0)
#        else:
#            source_mask.append(1)
#else:
#    #max is peak
#    #take top 50% - in value not in amount of point
#    midway_in_value=max_data.max()-abs(max_data.min()-max_data.max())/2
#    midway_in_value_from0=max_data.max()-max_data.max()/2
#    topthird_from0=max_data.max()-max_data.max()/3*1
#    cutoff=max_data.max()#midway_in_value_from0
#    #now find all sources where its above cutoff
#    source_mask=[]
#    for i in max_data:
#        if i <cutoff:
#            source_mask.append(0)
#        else:
#            source_mask.append(1)


#plot all sources selected as part of peak activation based on cutff
vertices_oi= stc.vertices[0][np.where(np.array(source_mask)==1)]  
vertices_oi=[vertno_max_i]
#plot all sources
#brain = Brain(hemi='lh', surf='inflated', subjects_dir=shared_dir,subject_id=mri_partic,initial_time=time_max)
brain=stc.plot(hemi='lh', subjects_dir=shared_dir,initial_time=time_max,surface='orig')
surf = brain.geo['lh']
vertidx = np.where(src[0]['inuse'])[0]
mlab.points3d(surf.x[vertidx], surf.y[vertidx],surf.z[vertidx], color=(.5,.5,.5), scale_factor=1.5)
#now plot the sources in our ROI
mlab.points3d(surf.x[vertices_oi], surf.y[vertices_oi],surf.z[vertices_oi], color=( 0,0,0), scale_factor=3)

#make waveform based on selected sources
data_indeces=[]

for v in vertices_oi:
    data_indeces.append(np.where(stc.vertices[0]==v))
    
data_indeces=np.array(data_indeces)    
data_in_vertices_oi=stc.data[data_indeces]   
data_in_vertices_oi=np.reshape(data_in_vertices_oi,[len(vertices_oi),len(stc._times)])
data_oi=np.mean(data_in_vertices_oi,0)


#Am to nAM
data_oi=[num*1e+9 for num in data_oi]


## now let's plot it with jones
data_x=stc._times
data_y=data_oi

jones=open('C:/Users/ckohl/Desktop/Current/jones_supra_extracted.txt','r')#ead()
jones_x=[]
jones_y=[]
for line in jones.readlines():
    sep=line.find('/t')
    jones_x.append(float(line[0:sep]))
    jones_y.append(float(line[sep+2:]))
#plt.plot(jones_x,jones_y)    

# fit x axis
if max(data_x) <5: #data is in seconds
    jones_x=[num /1000 for num in jones_x]


#fit scale
import scipy.optimize

def min_this(scale,args):
    import numpy as np
    abs_mag_data=max(data_y)-min(data_y)
    abs_mag_jones=max(jones_y)-min(jones_y)
    difference_to_minimise=np.abs(np.abs(abs_mag_jones)*scale-np.abs(abs_mag_data))
    return(difference_to_minimise)                         
scale=1;
res = scipy.optimize.minimize(min_this, scale, args=[data_y,jones_y], method='nelder-mead') #,options={'maxiter': 10000, 'disp': True})
jones_y_scaled=[num*res.x for num in jones_y]

#plot
fig, ax = plt.subplots()
plt.plot(data_x,data_y,'k')
jones_y_scaled=[num*-1 for num in jones_y_scaled]
plt.plot(jones_x,jones_y_scaled,':',color='grey')
ax.set_ylabel('Amplitude(nAm)')
ax.set_xlabel('Time(s)')
ax.set_title('Source Average')
ax.legend(['Source','Jones'])
textstr='Jones scaled by: '+ str(res.x)
props = dict(boxstyle='round', facecolor='gray', alpha=0.9)
ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
ylims=ax.get_ylim()
plt.xlim(-0.05,0.3)
plt.ylim(ylims)
plt.plot([0,0],ylims,':',color='grey')
plt.plot([0.02,0.02],ylims,':',color='grey')
plt.plot([0.04,0.04],ylims,':',color='grey')

plt.show()



#
#
#
## now let's do somethign similar (defining which sources around peak should be included), but only include them if they're eiwthin 3cm of peak source
## =======================
##this only works with an inflated brain. i think it takes the coords based on whether its been told if its inflater in the first definition of 'brain'. and if i dont inflate, that te spatial distance is different
#
##find max
#tmin=.00
#tmax=.1
#vertno_max, time_max = stc.get_peak(hemi='lh',tmin=tmin,tmax=tmax)
#vertno_max_i, time_max_i = stc.get_peak(hemi='lh',tmin=tmin,tmax=tmax,vert_as_index=True, time_as_index=True)
#brain=stc.plot(hemi='lh', subjects_dir=shared_dir,initial_time=time_max,surface='orig')
#surf = brain.geo['lh']
#vertidx = np.where(src[0]['inuse'])[0]
##plot all sources
#mlab.points3d(surf.x[vertidx], surf.y[vertidx],surf.z[vertidx], color=(.5,.5,.5), scale_factor=1.5)
##plot max
#vertix_oi=stc.vertices[0][vertno_max_i]
#mlab.points3d(surf.x[vertix_oi], surf.y[vertix_oi],surf.z[vertix_oi], color=( 1,0,0), scale_factor=3)
#
##define cutoff(like above)
#peak_is_negative=False
#lh_data=stc.data[:len(stc.vertices[0])]; #right: stc.data[len(stc.vertices[1]:)]
#max_data=lh_data[:,time_max_i]
#if lh_data[vertno_max_i,time_max_i]<0:
#    #max is trough
#    midway_in_value=max_data.min()+abs(max_data.min()-max_data.max())/2
#    midway_in_value_from0=max_data.min()-max_data.min()/2
#    cutoff=midway_in_value_from0
#    peak_is_negative=True
#
#else:
#    #max is peak
#    #take top 50% - in value not in amount of point
#    midway_in_value=max_data.max()-abs(max_data.min()-max_data.max())/2
#    midway_in_value_from0=max_data.max()-max_data.max()/2
#    topthird_from0=max_data.max()-max_data.max()/3*1
#    cutoff=midway_in_value_from0
#            
##find other sources based on cutoff and spatial distance
#vertices_oi=[]
#for source in vertidx:
#    if np.sum(abs(np.array([surf.x[source],surf.y[source],surf.z[source]]) - np.array([surf.x[vertix_oi], surf.y[vertix_oi],surf.z[vertix_oi]]))) <30:
#        #all sources within 3cm of max source
#        if peak_is_negative==False:
#            if stc.data[np.where(stc.vertices[0]==source),time_max_i]>cutoff:
#                vertices_oi.append(source)
#        else:
#            if stc.data[np.where(stc.vertices[0]==source),time_max_i]<cutoff:
#                vertices_oi.append(source)
##plot selected sources on brain            
#mlab.points3d(surf.x[vertices_oi], surf.y[vertices_oi],surf.z[vertices_oi], color=( 0,0,0), scale_factor=3)
#mlab.points3d(surf.x[vertix_oi], surf.y[vertix_oi],surf.z[vertix_oi], color=( 1,0,0), scale_factor=3)
#
##get waveform from those sources
#vertices_oi.append(vertix_oi)  
#data_indeces=[]
#for v in vertices_oi:
#    data_indeces.append(np.where(stc.vertices[0]==v))
#data_indeces=np.array(data_indeces)   
#data_in_vertices_oi=stc.data[data_indeces]   
#data_in_vertices_oi=np.reshape(data_in_vertices_oi,[len(vertices_oi),len(stc._times)])
#data_oi=np.mean(data_in_vertices_oi,0)
##Am to nAM
#data_oi=[num*1e+9 for num in data_oi]
#
#
#    
### now let's plot it with jones
#data_x=stc._times
#data_y=data_oi
#
#jones=open('C:/Users/ckohl/Desktop/Current/jones_supra_extracted.txt','r')#ead()
#jones_x=[]
#jones_y=[]
#for line in jones.readlines():
#    sep=line.find('/t')
#    jones_x.append(float(line[0:sep]))
#    jones_y.append(float(line[sep+2:]))
##plt.plot(jones_x,jones_y)    
#
## fit x axis
#if max(data_x) <5: #data is in seconds
#    jones_x=[num /1000 for num in jones_x]
#
#
##fit scale
#import scipy.optimize
#
#def min_this(scale,args):
#    import numpy as np
#    abs_mag_data=max(data_y)-min(data_y)
#    abs_mag_jones=max(jones_y)-min(jones_y)
#    difference_to_minimise=np.abs(np.abs(abs_mag_jones)*scale-np.abs(abs_mag_data))
#    return(difference_to_minimise)                         
#scale=1;
#res = scipy.optimize.minimize(min_this, scale, args=[data_y,jones_y], method='nelder-mead') #,options={'maxiter': 10000, 'disp': True})
#jones_y_scaled=[num*res.x for num in jones_y]
#
##plot
#fig, ax = plt.subplots()
#plt.plot(data_x,data_y,'k')
##jones_y_scaled=[num*-1 for num in jones_y_scaled]
#plt.plot(jones_x,jones_y_scaled,':',color='grey')
#ax.set_ylabel('Amplitude(nAm)')
#ax.set_xlabel('Time(s)')
#ax.set_title('Source Average')
#ax.legend(['Source','Jones'])
#textstr='Jones scaled by: '+ str(res.x)
#props = dict(boxstyle='round', facecolor='gray', alpha=0.9)
#ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
#
#
#ylims=ax.get_ylim()
#plt.xlim(-0.05,0.3)
#plt.ylim(ylims)
#plt.plot([0,0],ylims,':',color='grey')
#plt.plot([0.02,0.02],ylims,':',color='grey')
#plt.plot([0.04,0.04],ylims,':',color='grey')
#
#plt.show()
#
#   
#   
## what I'm doing here is nonsene. first of ll, if theres a few sources, avg or sum?
##second, if the area is across a sulcus, they should really cancel each other out in direction?so maybe just pickone
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
### ECD
## =============================================================================
## =============================================================================  
#
#
#########################################
#from mne.forward import make_forward_dipole
#from nilearn.plotting import plot_anat
#from nilearn.datasets import load_mni152_template
#
#evoked_full = evoked.copy()
##evoked=evoked_full
#evoked.crop(0.07, 0.1)
#
## Fit a dipole
#dip = mne.fit_dipole(evoked, cov, bem, trans)[0]
#
## Plot the result in 3D brain with the MRI image.
#dip.plot_locations(trans,mri_partic, shared_dir, mode='orthoview')
#
## Plot the result in 3D brain with the MRI image using Nilearn
## In MRI coordinates and in MNI coordinates (template brain)
#
#mni_pos = mne.head_to_mni(dip.pos, mri_head_t=trans,
#                          subject=mri_partic, subjects_dir=shared_dir)
#
#mri_pos = mne.head_to_mri(dip.pos, mri_head_t=trans,
#                          subject=mri_partic, subjects_dir=shared_dir)
#
#
#t1_fname = shared_dir+'\\'+ mri_partic+'\\'+'mri'+'\\'+'T1.mgz'
#fig_T1 = plot_anat(t1_fname, cut_coords=mri_pos[0], title='Dipole loc.')
#
#template = load_mni152_template()
#fig_template = plot_anat(template, cut_coords=mni_pos[0],
#                         title='Dipole loc. (MNI Space)')
#
#
#
#######
#from mne.simulation import simulate_evoked
#fwd_d, stc_d = make_forward_dipole(dip, bem, evoked.info, trans)
#pred_evoked = simulate_evoked(fwd_d, stc_d, evoked.info, cov=None, nave=np.inf)
#
## find time point with highest GOF to plot
#best_idx = np.argmax(dip.gof)
#best_time = dip.times[best_idx]
#print('Highest GOF %0.1f%% at t=%0.1f ms with confidence volume %0.1f cm^3'
#      % (dip.gof[best_idx], best_time * 1000,
#         dip.conf['vol'][best_idx] * 100 ** 3))
## remember to create a subplot for the colorbar
#fig, axes = plt.subplots(nrows=1, ncols=4, figsize=[10., 3.4])
#vmin, vmax = -400, 400  # make sure each plot has same colour range
#
## first plot the topography at the time of the best fitting (single) dipole
#plot_params = dict(times=best_time, ch_type='eeg', outlines='skirt',
#                   colorbar=False, time_unit='s')
#evoked.plot_topomap(time_format='Measured field', axes=axes[0], **plot_params)
#
## compare this to the predicted field
#pred_evoked.plot_topomap(time_format='Predicted field', axes=axes[1],
#                         **plot_params)
#
## Subtract predicted from measured data (apply equal weights)
##diff = combine_evoked([evoked, -pred_evoked], weights='equal')
##plot_params['colorbar'] = True
##diff.plot_topomap(time_format='Difference', axes=axes[2], **plot_params)
##plt.suptitle('Comparison of measured and predicted fields '
#             #'at {:.0f} ms'.format(best_time * 1000.), fontsize=16)
#####
#
#dip_fixed = mne.fit_dipole(evoked_full, cov, bem, trans,
#                           pos=dip.pos[best_idx], ori=dip.ori[best_idx])[0]
#dip_fixed.plot(time_unit='s')
#
#dip_fixed.plot()
#
#
#
#
#
#
#
## make own waveform plot
##get waveform from those sources
#data_oi=dip_fixed.data[0]
##Am to nAM
#data_oi=[num*1e+9 for num in data_oi]
#
#
#    
### now let's plot it with jones
#data_x=dip_fixed.times
#data_y=data_oi
#
#jones=open('C:/Users/ckohl/Desktop/Current/jones_supra_extracted.txt','r')#ead()
#jones_x=[]
#jones_y=[]
#for line in jones.readlines():
#    sep=line.find('/t')
#    jones_x.append(float(line[0:sep]))
#    jones_y.append(float(line[sep+2:]))
##plt.plot(jones_x,jones_y)    
#
## fit x axis
#if max(data_x) <5: #data is in seconds
#    jones_x=[num /1000 for num in jones_x]
#
#
##fit scale
#import scipy.optimize
#
#def min_this(scale,args):
#    import numpy as np
#    abs_mag_data=max(data_y)-min(data_y)
#    abs_mag_jones=max(jones_y)-min(jones_y)
#    difference_to_minimise=np.abs(np.abs(abs_mag_jones)*scale-np.abs(abs_mag_data))
#    return(difference_to_minimise)                         
#scale=1;
#res = scipy.optimize.minimize(min_this, scale, args=[data_y,jones_y], method='nelder-mead') #,options={'maxiter': 10000, 'disp': True})
#jones_y_scaled=[num*res.x for num in jones_y]
#
##plot
#fig, ax = plt.subplots()
#plt.plot(data_x,data_y,'k')
#jones_y_scaled=[num*-1 for num in jones_y_scaled]
#plt.plot(jones_x,jones_y_scaled,':',color='grey')
#ax.set_ylabel('Amplitude(nAm)')
#ax.set_xlabel('Time(s)')
#ax.set_title('Dipole')
#ax.legend(['Source','Jones'])
#textstr='Jones scaled by: -'+ str(res.x)
#props = dict(boxstyle='round', facecolor='gray', alpha=0.9)
#ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
#plt.xlim(-0.1,0.3)
#plt.show()
#
#   













# lets look at peak activation (only one max vertex)
# =======================
#set time window to find peak within certain times
tmin=.00
tmax=.1


peak_vertex, peak_time = stc.get_peak(hemi='lh', vert_as_index=True, time_as_index=True)
peak_vertex_surf = stc.lh_vertno[peak_vertex]
peak_value = stc.lh_data[peak_vertex, peak_time]

surfer_kwargs = dict(
    subjects_dir=shared_dir, views='lateral',
    initial_time=time_max, time_unit='s', size=(800, 800), smoothing_steps=5)


brain = stc.plot(**surfer_kwargs)


# We add the new peak coordinate (as vertex index) as an annotation dot
brain.add_foci(peak_vertex_surf, coords_as_verts=True, hemi='lh', color='blue')
# We add a title as well, stating the amplitude at this time and location
brain.add_text(0.1, 0.9, 'Peak coordinate', 'title', font_size=14)



vertno_max, time_max = stc.get_peak(hemi='lh',tmin=tmin,tmax=tmax)
#also get indices
vertno_max_i, time_max_i = stc.get_peak(hemi='lh',tmin=tmin,tmax=tmax,vert_as_index=True, time_as_index=True)

#show peak on brain
surfer_kwargs = dict(
    subjects_dir=shared_dir, views='lateral',
    initial_time=time_max, time_unit='s', size=(800, 800), smoothing_steps=5)
brain = stc.plot(**surfer_kwargs)
brain.add_foci(vertno_max, coords_as_verts=True, hemi='lh', color='black',
               scale_factor=1, alpha=0.9)
brain.add_text(0.1, 0.9, 'Peak between '+str(tmin*1000)+' and '+str(tmax*1000)+' ms: '+str(round(time_max*1000))+' ms', 'title',font_size=14)

#take activation at only max time and see which sources should be included in peak, based on a cutoff
lh_data=stc.data[:len(stc.vertices[0])]; #right: stc.data[len(stc.vertices[1]:)]
max_data=lh_data[:,time_max_i]


#plot all sources selected as part of peak activation based on cutff
vertices_oi= stc.vertices[0][np.where(np.array(source_mask)==1)]  
vertices_oi=[vertno_max_i]
#plot all sources
#brain = Brain(hemi='lh', surf='inflated', subjects_dir=shared_dir,subject_id=mri_partic,initial_time=time_max)
brain=stc.plot(hemi='lh', subjects_dir=shared_dir,initial_time=time_max,surface='orig')
surf = brain.geo['lh']
vertidx = np.where(src[0]['inuse'])[0]
mlab.points3d(surf.x[vertidx], surf.y[vertidx],surf.z[vertidx], color=(.5,.5,.5), scale_factor=1.5)
#now plot the sources in our ROI
mlab.points3d(surf.x[vertices_oi], surf.y[vertices_oi],surf.z[vertices_oi], color=( 0,0,0), scale_factor=3)

#make waveform based on selected sources
data_indeces=[]

for v in vertices_oi:
    data_indeces.append(np.where(stc.vertices[0]==v))
    
data_indeces=np.array(data_indeces)    
data_in_vertices_oi=stc.data[data_indeces]   
data_in_vertices_oi=np.reshape(data_in_vertices_oi,[len(vertices_oi),len(stc._times)])
data_oi=np.mean(data_in_vertices_oi,0)


#Am to nAM
data_oi=[num*1e+9 for num in data_oi]


## now let's plot it with jones
data_x=stc._times
data_y=data_oi

jones=open('C:/Users/ckohl/Desktop/Current/jones_supra_extracted.txt','r')#ead()
jones_x=[]
jones_y=[]
for line in jones.readlines():
    sep=line.find('/t')
    jones_x.append(float(line[0:sep]))
    jones_y.append(float(line[sep+2:]))
#plt.plot(jones_x,jones_y)    

# fit x axis
if max(data_x) <5: #data is in seconds
    jones_x=[num /1000 for num in jones_x]


#fit scale
import scipy.optimize

def min_this(scale,args):
    import numpy as np
    abs_mag_data=max(data_y)-min(data_y)
    abs_mag_jones=max(jones_y)-min(jones_y)
    difference_to_minimise=np.abs(np.abs(abs_mag_jones)*scale-np.abs(abs_mag_data))
    return(difference_to_minimise)                         
scale=1;
res = scipy.optimize.minimize(min_this, scale, args=[data_y,jones_y], method='nelder-mead') #,options={'maxiter': 10000, 'disp': True})
jones_y_scaled=[num*res.x for num in jones_y]

#plot
fig, ax = plt.subplots()
plt.plot(data_x,data_y,'k')
jones_y_scaled=[num*-1 for num in jones_y_scaled]
plt.plot(jones_x,jones_y_scaled,':',color='grey')
ax.set_ylabel('Amplitude(nAm)')
ax.set_xlabel('Time(s)')
ax.set_title('Source Average')
ax.legend(['Source','Jones'])
textstr='Jones scaled by: '+ str(res.x)
props = dict(boxstyle='round', facecolor='gray', alpha=0.9)
ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
ylims=ax.get_ylim()
plt.xlim(-0.05,0.3)
plt.ylim(ylims)
plt.plot([0,0],ylims,':',color='grey')
plt.plot([0.02,0.02],ylims,':',color='grey')
plt.plot([0.04,0.04],ylims,':',color='grey')

plt.show()





    t=100
    t_win_inds = np.where(abs(stc.times-t/1000/2)*2<(t/1000))[0][:] #Restrict vertex selection to the time window [0,t] 
    data_dspm_normal = stc.data 
    peak_vert = np.argmax(np.linalg.norm(data_dspm_normal[:, t_win_inds], axis=1)) #Vertex that contributes the most over window of time 
    peak_time = stc.times[t_win_inds[np.argmax(abs(data_dspm_normal[peak_vert, t_win_inds]))]] #Time at which abs(peak amplitude) of the selected vertex occurs 
     
    print('Peak vertex occurs at t =',peak_time) 
     
    # Plot vertex MNE timeseries 
    plt.figure() 
    plt.plot(1e3 * stc.times, 1e9 * stc.data[peak_vert, :].T, 'k-') 
    plt.xlabel('time (ms)') 
    plt.ylabel('current dipole (nAm)') 
    plt.xlim((0, 200)) 
    plt.axhline(0,c='k',ls=':') 
    plt.show() 
    #plt.savefig('n20_source_waveform.pdf')"

    # Set the components of all vertices besides the selected vertex to zero 
    stc_thresh = stc.copy() 
    subthresh_inds = [i for i in range(stc_thresh.data.shape[0]) if i!=peak_vert] 
    stc_thresh.data[subthresh_inds,:,:] = 1e-15 
    stc_thresh.data[peak_vert,:,:] = stc_thresh.data[peak_vert,:,:] * 1e9 
     
    # Find ID of selected vertex 
    if peak_vert >= len(stc.vertices[0]): 
        pick_vertID = stc.vertices[1][peak_vert-len(stc.vertices[0])] 
    else: 
        pick_vertID = stc.vertices[0][peak_vert]

    # Plot freesurfer brain with single vertex activation (40ms) 
    surfer_kwargs = dict( 
        #hemi='lh',  
        subjects_dir=shared_dir,  
        subject=mri_partic,  
        #clim=dict(kind='value', lims=[0,np.percentile(stc_mne.magnitude().data,95)*1e9,np.max(stc_mne.magnitude().data)*1e9]),  
        #brain_alpha=0.85,  
        #overlay_alpha=1,  
        initial_time=1e3 * peak_time, 
        views='dor', 
        time_unit='ms',  
        size=(1000, 1000), 
        #glyph='arrow3d', 
        colormap='Blues', 
        background='w', 
        foreground='k') 
    brain_40 = stc_thresh.plot(**surfer_kwargs) 
    brain_40.add_foci(pick_vertID, coords_as_verts=True, hemi='lh', color='yellow', scale_factor=0.2, alpha=0.6) 
    #brain_40.save_image('n20_source_40.jpg')



