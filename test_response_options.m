% Tests response options
% thumb movement vs eye movement

% for classification see Categorise_eye_movements.m

% Recorded on Beta02 - 11.11.2020
% EOG/EMG only


clear
close all

% load eeglab
eeglab_dir='C:\Users\ckohl\Documents\MATLAB\eeglab2020_0';
cd(eeglab_dir)
eeglab



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TMS session from Sep11, Beta02, single pulses
raw_path='C:\Users\ckohl\Downloads\';
filename= 'actiCHamp_Plus_BC-TMS_BETA02_20201111_EOG_response_test000036.vhdr';
dt=25;

EEG = pop_loadbv(raw_path, filename, [], []);
channels=EEG.chanlocs;
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw

%% First, let's see what triggerd there
count_triggers=[];
trigger_types={};
boundary_count=0;
for event = 1:length(EEG.event)-1
    if ~any(ismember(EEG.event(event).type ,trigger_types))
        trigger_types{end+1}=EEG.event(event).type;
        count_triggers(end+1)=1;
    else
        for trig=1:length(trigger_types)
            if EEG.event(event).type(1:2) == trigger_types{trig}(1:2)
                count_triggers(trig)=count_triggers(trig)+1;
            end
        end
    end
end

fprintf('Found a total of %i triggers and %i trigger types: \n',length(EEG.event),length(trigger_types))
for trig=1:length(trigger_types)
    fprintf('\t %s (%i) \n', trigger_types{trig},count_triggers(trig))
end
disp('----')
    

%% LEFT v RIGHT
EEG_keep=EEG;
EEG = pop_epoch( EEG_keep, { 'right', 'left'}, [-2 2], 'epochinfo', 'yes');
% need long epochs to get a proper mode baseline later
EEG = pop_rmbase( EEG, [-999   -500]);
% drop EEG channels
EEG=pop_select(EEG, 'nochannel',[1:63]);
EEG= pop_resample(EEG,50);
dt=.05;

left = pop_epoch( EEG, { 'left'}, [-2 2], 'epochinfo', 'yes');
right = pop_epoch( EEG, {'right'}, [-2 2], 'epochinfo', 'yes');
channel_oi=2;
rows=max([size(left.data,3),size(right.data,3)]);

figure
count=0;    
for trial=1:rows
    count=count+1;
    try  
        subplot(rows,2,count)
        plot(left.times,left.data(channel_oi,:,trial),'k','Linewidth',2) 
    end
    count=count+1;
    try
        subplot(rows,2,count)
        plot(right.times,right.data(channel_oi,:,trial),'k','Linewidth',2)
    end
end
subplot(rows,2,1)
title('Left')
subplot(rows,2,2)
title('Right')




%% UP v DOWN
EEG = pop_epoch( EEG_keep, { 'up', 'down'}, [-2 2], 'epochinfo', 'yes');
% need long epochs to get a proper mode baseline later
EEG = pop_rmbase( EEG, [-999   -500]);
% drop EEG channels
EEG=pop_select(EEG, 'nochannel',[1:63]);
EEG= pop_resample(EEG,50);
dt=.05;

up = pop_epoch( EEG, { 'up'}, [-2 2], 'epochinfo', 'yes');
down = pop_epoch( EEG, {'down'}, [-2 2], 'epochinfo', 'yes');
channel_oi=1;
rows=max([size(up.data,3),size(down.data,3)]);

figure
count=0;    
for trial=1:rows
    count=count+1;
    try  
        subplot(rows,2,count)
        plot(up.times,up.data(channel_oi,:,trial),'k','Linewidth',2) 
    end
    count=count+1;
    try
        subplot(rows,2,count)
        plot(down.times,down.data(channel_oi,:,trial),'k','Linewidth',2)
    end
end
subplot(rows,2,1)
title('Up')
subplot(rows,2,2)
title('Down')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EMG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TMS session from Sep11, Beta02, single pulses
raw_path='C:\Users\ckohl\Downloads\';
filename= 'actiCHamp_Plus_BC-TMS_BETA02_20201111_EOG_response_test000037.vhdr';
dt=25;

EEG = pop_loadbv(raw_path, filename, [], []);
channels=EEG.chanlocs;
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw

%% First, let's see what triggerd there
count_triggers=[];
trigger_types={};
boundary_count=0;
for event = 1:length(EEG.event)-1
    if ~any(ismember(EEG.event(event).type ,trigger_types))
        trigger_types{end+1}=EEG.event(event).type;
        count_triggers(end+1)=1;
    else
        for trig=1:length(trigger_types)
            if EEG.event(event).type(1:2) == trigger_types{trig}(1:2)
                count_triggers(trig)=count_triggers(trig)+1;
            end
        end
    end
end

fprintf('Found a total of %i triggers and %i trigger types: \n',length(EEG.event),length(trigger_types))
for trig=1:length(trigger_types)
    fprintf('\t %s (%i) \n', trigger_types{trig},count_triggers(trig))
end
disp('----')
    



%% LEFT v RIGHT
EEG = pop_epoch( EEG, { 'right', 'left'}, [-2 2], 'epochinfo', 'yes');
% need long epochs to get a proper mode baseline later
EEG = pop_rmbase( EEG, [-999   -500]);
% drop EEG channels
EEG=pop_select(EEG, 'nochannel',[1:63]);
EEG= pop_resample(EEG,50);
dt=.05;

left = pop_epoch( EEG, { 'left'}, [-2 2], 'epochinfo', 'yes');
right = pop_epoch( EEG, {'right'}, [-2 2], 'epochinfo', 'yes');

rows=max([size(left.data,3),size(right.data,3)]);

figure
count=0;    
for trial=1:rows
    count=count+1;
    try  
        subplot(rows,2,count)
        hold on
        plot(left.times,left.data(1,:,trial),'r','Linewidth',2) 
        plot(left.times,left.data(2,:,trial),'g','Linewidth',2)
    end
    count=count+1;
    try
        subplot(rows,2,count)
        hold on
        plot(right.times,right.data(1,:,trial),'r','Linewidth',2)
        plot(right.times,right.data(2,:,trial),'g','Linewidth',2)
    end
end
subplot(rows,2,1)
title('Left')
subplot(rows,2,2)
title('Right')

