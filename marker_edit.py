# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 12:40:59 2019

@author: ckohl
"""
# Function to edit EEG markers
# -  Spits out overview of all markers
# -  Makes new burst markers (optional)
#    Burst defined as 3 'S 13' markers within 50ms (assumes 'S 15' markers in between)
#    Creates new labels 'S 131', 'S 132', 'S 133) for each pulse in burst.
#    Does not edit original markers.
# -  Deletes reset markers 'S 15' (optional)
def marker_edit(raw):

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
    