% Tests response options
% thumb movement vs eye movement

% Recorded on Beta02 - 11.11.2020
% EOG/EMG only


clear
close all

% load eeglab
eeglab_dir='C:\Users\ckohl\Documents\MATLAB\eeglab2020_0';
cd(eeglab_dir)
eeglab



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD RAW
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPROC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% epoch
EEG = pop_epoch( EEG, { 'right', 'left'}, [-2 2], 'epochinfo', 'yes');
% need long epochs to get a proper mode baseline later
EEG = pop_rmbase( EEG, [-999   -500]);

% drop EEG channels
EEG=pop_select(EEG, 'nochannel',[1:63]);
EEG= pop_resample(EEG,50);
dt=.05;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CATEGORISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Review=true;
channel_oi=2;
threshold=100;
class=nan(size(EEG.data,3),1);
base=nan(size(squeeze(EEG.data(channel_oi,:,:))));
edited_base=zeros(size(squeeze(EEG.data(channel_oi,:,:))));
Class_string={'right','left'};
Class_colours={[.625 .875 .625],[.906 .625 .344]};
for trial = 1:size(EEG.data,3)
    channel_oi=2;
    %find extremes
    base(:,trial)=abs(EEG.data(channel_oi,:,trial)- mean(EEG.data(channel_oi,:,trial)))>threshold;
    %base is just identifyign extreme values, so let's keep only the
    %longest extreme interval    
    base_i=1;
    longest_int=[];
    while base_i < length(base(:,trial))
        if base(base_i,trial)
           stop=min(find(base([base_i+1:end],trial)==false));
           if length(base_i:base_i+stop-1)>length(longest_int)
               longest_int=base_i:base_i+stop-1;
           end
           base_i=base_i+stop;
        else
            base_i=base_i+1;
        end
    end
    % keep inly longest int
    edited_base(longest_int,trial)=true;             

%     subplot(5,6,trial)
%     yyaxis left
%     plot(EEG.data(channel_oi,:,trial),'Color',[.5 .5 .5])
%     hold on
%     yyaxis right
%     plot(base(:,trial),'r-')
%     plot(edited_base(:,trial),'k-','Linewidth',2 )
    
    % classify
    movement=mean(EEG.data(channel_oi,longest_int,trial));
    background=mean(EEG.data(channel_oi,:,trial));
    
    if movement > background
        class(trial)=1;
    else
        class(trial)=2;
    end
    
    if Review & (trial/9==round(trial/9) | trial==size(EEG.data,3))
        figure('units','normalized','outerposition',[0 0 1 1])
        count=0;
        for review_trial=trial-8:trial 
            count=count+1;
            sb=subplot(3,3,count);
            hold on
            title(strcat('Trial ',num2str(review_trial)))
            
            % Data
            plot(EEG.times, EEG.data(channel_oi,:,review_trial),'-','Color',[.2 .2 .2],'Linewidth',2)
            yyaxis right
            % all extreme values
            plot(EEG.times,base(:,review_trial),'-','Color',Class_colours{class(review_trial)});
            % identified saccade
            ar=area(EEG.times,edited_base(:,review_trial));
            ar.EdgeColor=Class_colours{class(review_trial)};
            ar.FaceColor=Class_colours{class(review_trial)};
            ar.FaceAlpha=.5;
            % classification
            display_string=Class_string{class(review_trial)};
            tx=text(EEG.times(1)+ EEG.times(end)/10,.4,display_string);
            tx.Color=Class_colours{class(review_trial)};
            tx.FontSize=15;      
        end

        decision=input('All ok? (Yes: Enter, No: 1)');
        if decision==1
            change_trials=input('Which trials are wrong? (Enter vector)\n');
            for change=change_trials
                class(change)=class(change)-1;
                if class(change)==0
                    class(change)=2;
                end
            end
        end
    end
end

