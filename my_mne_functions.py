# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:16:51 2019

@author: ckohl
"""

# MNE Functions to make my life easier. Put them all in here so I can import them all at once. Do not like.




def marker_edit(raw):
# Function to edit EEG markers
# -  Spits out overview of all markers
# -  Makes new burst markers (optional)
#    Burst defined as 3 'S 13' markers within 50ms (assumes 'S 15' markers in between)
#    Creates new labels 'S 131', 'S 132', 'S 133) for each pulse in burst.
#    Does not edit original markers.
# -  Deletes reset markers 'S 15' (optional)
# 'C:/Users/ckohl/Documents/MNE Scripts/marker_edit.py'

    TMS_str='Stimulus/S 13'
    Reset_str='Stimulus/S 15'
    
    # Show markers
    all_events=[];
    all_events_count=[];
    for annot in raw.annotations:
        if any(annot['description']  in s for s in all_events) == False:
            all_events.append(annot['description'])
            all_events_count.append(1)
        else:
            all_events_count[all_events.index(annot['description'])]=all_events_count[all_events.index(annot['description'])]+1
    print('--- Events ---')        
    print('ID\t\tCount')
    for e in range(0,len(all_events)):
        print(all_events[e]+'\t'+ str(all_events_count[e]))
        
    # TMS
    if any(TMS_str in s for s in all_events)    :
        print('\nFound TMS markers ('+TMS_str+').')
        new_tms=input('Would you like to create new burst markers (y/n)?')
        if new_tms=='y':
            print('Creating new markers for bursts (S 131/132/133).')
            burst_count=0;
            new_onsets=[];
            new_durations=[];
            new_descriptions=[];
            for annot_i in range(0,len(raw.annotations)-4):
                if raw.annotations[annot_i]['description']==TMS_str: #if this trigger is TMS
                    if annot_i==0 or raw.annotations[annot_i-1]['description']!=TMS_str: #and last one wasnt
                        if raw.annotations[annot_i+2]['description']==TMS_str and raw.annotations[annot_i+4]['description']==TMS_str and raw.annotations[annot_i+4]['onset']-raw.annotations[annot_i]['onset']<.05 : # is burst
                            burst_count+=1
                            new_onsets.extend([raw.annotations[annot_i]['onset'],raw.annotations[annot_i+2]['onset'],raw.annotations[annot_i+4]['onset']])
                            new_durations.extend([raw.annotations[annot_i]['duration'],raw.annotations[annot_i+2]['duration'],raw.annotations[annot_i+4]['duration']])
                            new_descriptions.extend(['S 131','S 132', 'S133'])
            raw.annotations.append(onset=new_onsets, duration=new_durations,description=new_descriptions)    
            print('Found and created ' + str(burst_count) +' bursts.')               
    else:
        print('No TMS markers ('+TMS_str+') found.')
    
    if any(Reset_str in s for s in all_events)    :    
        dels15=input('Would you like to delete ' + Reset_str +' markers (y/n}?')   
        if dels15=='y':
            s15_i=[]
            for annot_i in range(0,len(raw.annotations)):
                if raw.annotations[annot_i]['description']==Reset_str:
                    s15_i.append(annot_i)
            raw.annotations.delete(s15_i)
            print('Reset markers removed.')
        else:
            print('Reset markers left in.')    
                
    return raw            
    



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
    






def import_dig_electrodes(electr_file_str,raw):
    # Function to import digitised electrode files
# - wants the .txt file!
# - checks if everything we need is in the file (63 channels,3 fiducials)
# - makes inputs for mne.channels.DigMontage
#   - dictionary of channels ('label': [X, Y, Z])
#   - array of headshape coords (n_channels x 3)
#   - one vector for each fiducial (1 x 3)
# - plots montage
# - let's you manually redefine unclear channels if applicable
#   - asks if you want to delete or re-label each unclear electrode
#   - plots updated montage
# - returns montage
       
    import mne
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from mpl_toolkits.mplot3d import Axes3D
    electr_file=open(electr_file_str,'r') 
    print('\n------------------------------------')
    print('-- Scanning Electrode File -- ')   
    print('------------------------------------')
    electr_lines=electr_file.readlines()
    comment_lines=[]
    fiducial_lines=[]
    channel_lines=[]
    head_lines=[]
    for line in range(0,len(electr_lines)):
        if electr_lines[line][0]=='#':
            comment_lines.append(line)
        elif electr_lines[line][0:4]=='Nasi' or electr_lines[line][0:4]=='Left' or electr_lines[line][0:4]=='Righ' :
            fiducial_lines.append(line)
        elif electr_lines[line][0]=='h'or electr_lines[line][0]=='H':
            head_lines.append(line)
        else:
            channel_lines.append(line)
    print('Found: \n'+str(len(comment_lines))+'\t Header lines \n'+str(len(fiducial_lines))+
          '\t Fiducials \n'+str(len(channel_lines))+'\t Channels \n'+str(len(head_lines))+
          '\t Head Shape Points \n'+str(len(electr_lines)-len(comment_lines)-len(fiducial_lines)
          -len(channel_lines)-len(head_lines)) +'\t Unidentified lines')     
    cont=input('Ok (y/n)?')
    if cont=='n':
        sys.exit()
        
              
    # make headshape array
    Headshape=np.empty((len(head_lines),3))
    for point in range(0,len(head_lines)):  
        first_tab_i=electr_lines[head_lines[point]].find('\t')
        second_tab_i=electr_lines[head_lines[point]].find('\t',first_tab_i+1, len(electr_lines[head_lines[point]]))
        third_tab_i=electr_lines[head_lines[point]].find('\t',second_tab_i+1, len(electr_lines[head_lines[point]]))
        fourth_tab_i=electr_lines[head_lines[point]].find('\t',third_tab_i+1, len(electr_lines[head_lines[point]]))
        fifth_tab_i=electr_lines[head_lines[point]].find('\t',fourth_tab_i+1, len(electr_lines[head_lines[point]]))
        last_tab_i=electr_lines[head_lines[point]].find('\n')  
        X=float(electr_lines[head_lines[point]][third_tab_i+1:fourth_tab_i])/1000
        Y=float(electr_lines[head_lines[point]][fourth_tab_i+1:fifth_tab_i])/1000
        Z=float(electr_lines[head_lines[point]][fifth_tab_i+1:last_tab_i])/1000
        Headshape[point,:]=[X,Y,Z]
    
    # make fiducial variables
        Nasion=[]
        LPA=[]
        RPA=[]
    for fids in range(0,len(fiducial_lines)):
        first_tab_i=electr_lines[fiducial_lines[fids]].find('\t')
        second_tab_i=electr_lines[fiducial_lines[fids]].find('\t',first_tab_i+1, len(electr_lines[fiducial_lines[fids]]))
        third_tab_i=electr_lines[fiducial_lines[fids]].find('\t',second_tab_i+1, len(electr_lines[fiducial_lines[fids]]))
        fourth_tab_i=electr_lines[fiducial_lines[fids]].find('\t',third_tab_i+1, len(electr_lines[fiducial_lines[fids]]))
        fifth_tab_i=electr_lines[fiducial_lines[fids]].find('\t',fourth_tab_i+1, len(electr_lines[fiducial_lines[fids]]))
        last_tab_i=electr_lines[fiducial_lines[fids]].find('\n')  
        X=float(electr_lines[fiducial_lines[fids]][third_tab_i+1:fourth_tab_i])/1000
        Y=float(electr_lines[fiducial_lines[fids]][fourth_tab_i+1:fifth_tab_i])/1000
        Z=float(electr_lines[fiducial_lines[fids]][fifth_tab_i+1:last_tab_i])/1000
        if electr_lines[fiducial_lines[fids]][0:4]=='Nasi':
            Nasion=[X,Y,Z]
        elif electr_lines[fiducial_lines[fids]][0:4]=='Left':
            LPA=[X,Y,Z]
        elif electr_lines[fiducial_lines[fids]][0:4]=='Righ':
            RPA=[X,Y,Z]
        else:
            print('Fiducial Labels not recognised')
            sys.exit()
    
    
    # make electrode dict
    print('\n------------------------------------')
    print('-- Matching Electrodes -- ')   
    print('------------------------------------')
    Electr_dict={}
    Define_dict={}
    Label_list=[]
    no_match=0
    double=0
    double_list=[]
    missing_list=[]
    for channel in range(0,len(channel_lines)):
        first_tab_i=electr_lines[channel_lines[channel]].find('\t')
        second_tab_i=electr_lines[channel_lines[channel]].find('\t',first_tab_i+1, len(electr_lines[channel_lines[channel]]))    
        third_tab_i=electr_lines[channel_lines[channel]].find('\t',second_tab_i+1, len(electr_lines[channel_lines[channel]]))    
        fourth_tab_i=electr_lines[channel_lines[channel]].find('\t',third_tab_i+1, len(electr_lines[channel_lines[channel]]))    
        fifth_tab_i=electr_lines[channel_lines[channel]].find('\t',fourth_tab_i+1, len(electr_lines[channel_lines[channel]]))    
        last_tab_i=electr_lines[channel_lines[channel]].find('\n')  
        label=electr_lines[channel_lines[channel]][0:first_tab_i]
        #make sure spelling matches (this was needed before. not sure if it still is
        if label=='32':
            label='FCZ'
        elif label=='31':
            label='AFZ'  
        #match    
        label=label.upper()
        found=False
        for ch in raw.info['ch_names']:
            if ch.upper()==label:
                label=ch
                found=True
        if found==False:       
           print('No label match found: '+label)
           no_match+=1 
           
        X=float(electr_lines[channel_lines[channel]][third_tab_i+1:fourth_tab_i])/1000
        Y=float(electr_lines[channel_lines[channel]][fourth_tab_i+1:fifth_tab_i])/1000
        Z=float(electr_lines[channel_lines[channel]][fifth_tab_i+1:last_tab_i])/1000
        
        #delete doubles
        double=False
        for existing in Label_list:
            if label==existing:
                double=True   
                double_list.append(label)
                

        if double==False and found==True:
            Electr_dict[label]=[X,Y,Z]
            Label_list.append(label)
        else:
            Define_dict[label]=[X,Y,Z]
            
   
    print('Found (and deleted) '+str(len(double_list)) +' doubles:')
    for d in double_list:
        print('\t' +d)
    
    if len(Label_list)<63:
        for ch in raw.info['ch_names']:
            found=False
            for existing in Label_list:
                if existing==ch:
                    found=True
            if found==False:
                missing_list.append(ch)
    print(str(len(missing_list))+ ' channels missing from electrode file:')
    for m in missing_list:
        print('\t' +m)
    print(str(len(Label_list))+  ' channels matched.')
        
    
    # plot electrodes
    print('\n------------------------------------')
    print('-- Plotting Montage -- ')   
    print('------------------------------------')
    fig = plt.figure()
    ax1=fig.add_subplot(121)
    plt.imshow(mpimg.imread('C:/Users/ckohl/Desktop/Current/montage.JPG'))
    ax1.title.set_text('Montage')
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    ax2 = fig.add_subplot(122, projection='3d')
    for i in range(len(Electr_dict)+len(Define_dict)): # plot each point + it's index as text above
      try:
          label = Label_list[i]
          x= Electr_dict[label][0]
          y= Electr_dict[label][1]
          z= Electr_dict[label][2]
          leg1=ax2.scatter(x,y,z , color='b',s=10)
          ax2.text(x,y,z, '%s' % (label), size=10, zorder=1, color='k')
      except:
          label=list(Define_dict.keys())[i-len(Electr_dict)]
          x= Define_dict[label][0]
          y= Define_dict[label][1]
          z= Define_dict[label][2]  
          leg2=ax2.scatter(x,y,z , color='r',s=10)
          ax2.text(x,y,z, '%s' % (label), size=10, zorder=1, color='r')
    leg3=ax2.scatter(Nasion[0],Nasion[1],Nasion[2], color='g',s=50)    
    ax2.text(Nasion[0],Nasion[1],Nasion[2], '%s' % ('Nasion'), size=10, zorder=1, color='g')
    ax2.scatter(LPA[0],LPA[1],LPA[2], color='g',s=50)    
    ax2.text(LPA[0],LPA[1],LPA[2], '%s' % ('LPA'), size=10, zorder=1, color='g')
    ax2.scatter(RPA[0],RPA[1],RPA[2], color='g',s=50)    
    ax2.text(RPA[0],RPA[1],RPA[2], '%s' % ('RPA'), size=10, zorder=1, color='g')
    ax2.title.set_text('Imported Electrodes')
    ax2._axis3don = False
    try:
        plt.legend((leg1, leg2, leg3),('Regognised Electrodes','Unclear Electrodes','Fiducials'))
    except:
        try:
            plt.legend((leg1, leg3),('Regognised Electrodes','Fiducials'))
        except:
            plt.legend((leg1),('Regognised Electrodes'))
    fig.tight_layout()
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()
    print('Press Enter to continue (in figure).')
    while True:
        if plt.waitforbuttonpress():
            break
    
    
    # Define unclear channels
    if len(Define_dict)>0:
        print('\n------------------------------------')
        print('-- Define Electrodes -- ')   
        print('------------------------------------')
        print ('\nLets go through the ones that havent been matched: \n')
        for i in range(len(Define_dict)):
            decision=False
            while decision==False:
                display=input('--' +list(Define_dict.keys())[i]+' -- \n (Press y to interact with figure, or Enter to continue)')
                if display=='y':
                    print('\tPress enter to exit figure.')
                    #replot
                    while True:
                        if plt.waitforbuttonpress():
                            break
                dr=input(list(Define_dict.keys())[i]+' - Delete or Re-label? (d/r):')           
                if dr=='d':
                    print('\t'+list(Define_dict.keys())[i]+' deleted')
                    decision=True
                elif dr=='r':
                    decision=True
                    matched=False
                    while matched==False:
                        new_label=input('\t New label: ')
                        new_label=new_label.upper()
                        found=False
                        for ch in raw.info['ch_names']:
                            if ch.upper()==new_label:
                                new_label=ch
                                found=True
                        if found==True:
                            matched=True
                            Label_list.append(new_label)
                            Electr_dict[new_label]=Define_dict[list(Define_dict.keys())[i]]
                            print('\t' +new_label+ ' has been re-labeled and added to montage.')
                        else:
                            print('\t No match for the suggested label '''+new_label +' '' has been found. Please re-enter.')           
                else:
                    print('\t Error: Input must either be ''d'' (delete)  or ''r'' (re-label)')     
        #plot redefined montage
        print(str(len(Label_list))+  ' channels matched.')
        print('\n------------------------------------')
        print('-- Plotting Redefined Montage -- ')   
        print('------------------------------------')
        ax2.cla()
        for i in range(len(Electr_dict)): # plot each point + it's index as text above
            label = Label_list[i]
            x= Electr_dict[label][0]
            y= Electr_dict[label][1]
            z= Electr_dict[label][2]
            leg1=ax2.scatter(x,y,z , color='b',s=10)
            ax2.text(x,y,z, '%s' % (label), size=10, zorder=1, color='k')
        leg3=ax2.scatter(Nasion[0],Nasion[1],Nasion[2], color='g',s=50)    
        ax2.text(Nasion[0],Nasion[1],Nasion[2], '%s' % ('Nasion'), size=10, zorder=1, color='g')
        ax2.scatter(LPA[0],LPA[1],LPA[2], color='g',s=50)    
        ax2.text(LPA[0],LPA[1],LPA[2], '%s' % ('LPA'), size=10, zorder=1, color='g')
        ax2.scatter(RPA[0],RPA[1],RPA[2], color='g',s=50)    
        ax2.text(RPA[0],RPA[1],RPA[2], '%s' % ('RPA'), size=10, zorder=1, color='g')
        ax2.title.set_text('Imported Electrodes')
        ax2._axis3don = False
        plt.legend((leg1, leg3),('Regognised Electrodes','Fiducials'))
        fig.tight_layout()
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()
        print('Press Enter to continue (in figure).')
        while True:
            if plt.waitforbuttonpress():
                break
    
            
    plt.close()
    montage=mne.channels.make_dig_montage(hsp=Headshape, nasion=Nasion, lpa=LPA,rpa=RPA,ch_pos=Electr_dict)
  
    return(montage)