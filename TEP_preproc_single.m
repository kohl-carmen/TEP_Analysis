%TEP Preprocessing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

partic=2; %[2,4]
TESAICA=0; %if 0, runICA for blink, if 1 TESAICA automatic
refilt=0;
plot_steps=0;
plot_all_elecs=0;
S1_80=1;

% load eeglab
eeglab_dir='C:\Users\ckohl\Documents\MATLAB\eeglab2020_0';
cd(eeglab_dir)
eeglab

Partic=[2,4];
dt=25;
electr_oi='C3';

pulse_zero_time=[-2,5];%[-.1,2];
recharge_zero_time=[-.5,.5];
plot_times=[-5 20];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD RAW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TMS session from Sep11, Beta02, single pulses
raw_path='C:\Users\ckohl\Desktop\Current\TMS\CurrentData\EEG';

if S1_80
    filename= 'actiCHamp_Plus_BC-TMS_BETA02_20200911_EOG000032.vhdr';
else %S1_100
    filename= 'actiCHamp_Plus_BC-TMS_BETA02_20200911_EOG000033.vhdr';
end
% 
% % M1 stimulation - thresholding 
% % (corresponding to pulses 1-41 in Brainsight file)
% filename= 'actiCHamp_Plus_BC-TMS_BETA02_20200911_EOG000031.vhdr';
% % S1 stimulation - 80% threshold 
% % (corresponding to pulses 42-151 in Brainsight file)
% filename= 'actiCHamp_Plus_BC-TMS_BETA02_20200911_EOG000032.vhdr';
% % S1 stimulation - 100% threshold
% % (corresponding to pulses 152-201 in Brainsight file)
% filename= 'actiCHamp_Plus_BC-TMS_BETA02_20200911_EOG000033.vhdr';
% % S1 stimulation - 100% threshold - opposite direction of current - no neuronavigation - don't use this
% % (corresponding to pulses 202-206 in Brainsight file)
% filename= 'actiCHamp_Plus_BC-TMS_BETA02_20200911_EOG000034.vhdr';
% % S1 stimulation - 100% threshold - opposite direction of current - with neuronavigation, but no good?
% % (corresponding to pulses 207-212 in Brainsight file)
% filename= 'actiCHamp_Plus_BC-TMS_BETA02_20200911_EOG000035.vhdr';

EEG = pop_loadbv('C:\Users\ckohl\Desktop\Current\TMS\CurrentData\EEG\', filename, [], []);


channels=EEG.chanlocs;%for later
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw

%% delete empty intervals
% if partic==2
%     %for S1_80, but no longer necessry
% 	EEG = eeg_eegrej( EEG, [341 3225725;4621908 6973467;9305976 11206768;13518123 15378218;17828242 17962334;18376857 20565749;23383042 25310059;27935367 28506000]);
% end
%EEG.data=detrend(EEG.data);

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
            if EEG.event(event).type(1:4) == trigger_types{trig}(1:4)
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
    

%find electrode
for chan= 1:length(EEG.chanlocs)
    if length(EEG.chanlocs(chan).labels)==2
        if EEG.chanlocs(chan).labels==electr_oi
            electr_oi_i=chan;
        end
    end
end

%% epoch
EEG = pop_epoch( EEG, { 'S 13'}, [-2 2], 'epochinfo', 'yes');
EEG = pop_rmbase( EEG, [-500   -100]);
keep_EOG=EEG.data(64:65,:,:);

%plot
if plot_steps==true
    preprocfig=figure;
    nr_trials_to_plot=5;
    %select a few trials at random to look at
    trials_to_plot=sort(randi(size(EEG.data,3),1,nr_trials_to_plot));
    for trial =1:nr_trials_to_plot
        subplot(nr_trials_to_plot,1,trial)
        hold on
        plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trials_to_plot(trial)))
    end
end
%zoom into artifact
% if plot_steps==true
%     figure;
%     zoom_plot_times=[-1 20];
%     nr_trials_to_plot=5;
%     %select a few trials at random to look at
%     trials_to_plot=sort(randi(size(EEG.data,3),1,nr_trials_to_plot));
%     for trial =1:nr_trials_to_plot
%         hold on
%         plot(EEG.times(find(EEG.times==zoom_plot_times(1)):find(EEG.times==zoom_plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==zoom_plot_times(1)):find(EEG.times==zoom_plot_times(2))],trials_to_plot(trial)))
%     end
%     xlabel('Time (ms)')
%     ylabel('Ampltidue (µV)')
%     ylim([-200 200])
% end

if plot_all_elecs==true
    allelecfig=figure;
    for elec=1:size(EEG.data,1)
        trial=5;
        subplot(round(size(EEG.data,1)/4),4,elec)
        hold on
        set(gca,'XTickLabel',[], 'YTickLabel', [])
        plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(elec,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))
        eleclabel=text(min(xlim)+1,double(min(ylim)+(max(ylim)-min(ylim))/2),EEG.chanlocs(elec).labels)
    end
end

keep_raw_EEG=EEG;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%detection: if the value of a sample is further than 100 sds away from the
%mean of the some base
tms_time_wndw=[-20 80];% this is where I'll look for pulses
scale=100;% this is how much bigger i expect the signal to be at pulseonset
refract=10;% this is how long after a detected pulse I won't look for a new one
base=[-500: -100];
base_i=[find(EEG.times==base(1)): find(EEG.times==base(end))];
pulse_onset=[];
found_pulse=zeros(1,size(EEG.data,3));

for trial=1:size(EEG.data,3)
    basemean=mean(EEG.data(electr_oi_i,base_i,trial));
    basesd=std(EEG.data(electr_oi_i,base_i,trial));
    just_found_one=0;
    for t=1:length(EEG.times)
        if just_found_one>refract*dt
            just_found_one=0;
        elseif just_found_one >0
            just_found_one=just_found_one+1;
        end
        if EEG.times(t)>tms_time_wndw(1) & EEG.times(t)<tms_time_wndw(2) & just_found_one==0
            if  abs(EEG.data(electr_oi_i,t+1,trial))> basemean+basesd*scale%EEG.data(electr_oi_i,t+1,trial)> EEG.data(electr_oi_i,t+1,trial).*10
                pulse_onset(trial)=EEG.times(t);
                just_found_one=just_found_one+1;
                found_pulse(trial)=1;
            end
        end
    end
    if found_pulse(trial)==0
        fprintf(' Couldnt find pulse in trial %d\n',trial)
    end
end
%sanity check - find times that are different from the rest
pulsemean=mean(pulse_onset);
pulsesd=2;%std(pulse_onset).*2; % I just put no further away than 2ms rather than an sd caus sds are tiny
redo_trials=[];
redo_trials=[redo_trials; find(pulse_onset>pulsemean+pulsesd | pulse_onset<pulsemean-pulsesd )];
if ~isempty(redo_trials)
    fprintf('%d trials were found to have odd pulse times \n',size(redo_trials,1))
    pulse_onset(redo_trials,:)
else
    fprintf('Pulses where found for all trials.\nMean pulse time is:\t%d\t%d\t%d\n',round(pulsemean))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate Pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zero=false;
interp=true;
count=0;
for trial=1:size(EEG.data,3)
    pulse=[find(EEG.times==pulse_onset(trial))+round(pulse_zero_time(1)*dt): find(EEG.times==pulse_onset(trial))+round(pulse_zero_time(2)*dt)];
    
    if interp==true
        for electr=1:size(EEG.data,1)
        %     clf
        %     plot(EEG.times(999*dt:1020*dt), EEG.data(electr,[999*dt:1020*dt],1))
                hold on
                y=EEG.data(electr,:,trial);
                x=EEG.times;
                x([pulse])=[];
                y([pulse])=[];
                xx=EEG.times([pulse]);   
                yy=interp1(x,y,xx,'pchip');%pchip is cubib
                EEG.data(electr,[pulse],trial)=yy;
        %     plot(xx, yy,'.')
         end
    elseif zero==true
        EEG.data(:,[pulse],trial)=0;
    %plot(EEG.times(pulse2),zeros(size(pulse1)))
    end
    if plot_steps==true
        figure(preprocfig)
        if any(trial==trials_to_plot)
            count=count+1;
            subplot(nr_trials_to_plot,1,count)
            hold on
            plot(pulse_onset(trial),EEG.data(electr_oi_i,[find(EEG.times==pulse_onset(trial))],trial),'ro')
            plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))
            ylim([min(EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial)),max(EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))])
        end
    end
    if plot_all_elecs==true
        figure(allelecfig)
        if trial==size(EEG.data,3)
            for elec=1:size(EEG.data,1)
                trial=5;               
                subplot(round(size(EEG.data,1)/4),4,elec)
                hold on
                set(gca,'XTick',[], 'YTick', [])
                plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(elec,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))
                ylim([min(EEG.data(elec,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial)),max(EEG.data(elec,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))])       
                set(eleclabel,'Position',[min(xlim)+1,double(min(ylim)+(max(ylim)-min(ylim))/2)])
            end
        end
    end
end
EEG_minus_pulse=EEG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find recharge artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tms_time_wndw=[4 100];% this is where I'll look for pulses
if partic==2
    scale=6;% this is how much bigger i expect the signal to be at pulseonset
elseif partic==4
    scale=4;
end
refract=10;% this is how long after a detected pulse I won't look for a new one
base=[-500: -100];
base_i=[find(EEG.times==base(1)): find(EEG.times==base(end))];
recharge_onset=[];
missing=[];
a=[1 1/5];
b=[1/5 0];
for trial=1:size(EEG.data,3)
%     The
    temp=filter(b,a,EEG.data(electr_oi_i,:,trial));
    temp=temp-temp(find(EEG.times==-5));
    basesd=std(temp(base_i));
%     if basesd>10
%         basesd=10;
%     end
    found_recharge=[0];
    just_found_one=0;
    for t=1:length(EEG.times)
        if just_found_one>refract*dt
            just_found_one=0;
        elseif just_found_one >0
            just_found_one=just_found_one+1;
        end
%         t=t+1
%         EEG.times(t)
%         EEG.data(electr_oi_i,t+1,trial)
        if EEG.times(t)>tms_time_wndw(1) & EEG.times(t)<tms_time_wndw(2) & just_found_one==0
            if  abs(temp(t+1))> basesd*scale %& abs(EEG.data(electr_oi_i,t+1,trial))-abs(EEG.data(electr_oi_i,t,trial))>20
                recharge_onset(trial)=EEG.times(t);
                found_recharge=1;
                just_found_one=just_found_one+1;
            end
        end
    end
    if sum(found_recharge)<1
        missing=[missing,trial];
        fprintf(' Couldnt find recharge in trial %d\n',trial)
    end
end

% check if those times make sense
rechargemean=median(recharge_onset);
rechargesd=[2];%std(pulse_onset).*2; % I just put no further away than 2ms rather than an sd caus sds are tiny
redo_trials=[];
redo_trials=[redo_trials; find(recharge_onset>rechargemean+rechargesd | recharge_onset<rechargemean-rechargesd )];
% if not, change scale and go through the whole thing again fo those trial
if ~isempty(redo_trials)
    fprintf('%d trials were found to have odd recharge times \n',size(redo_trials))
    %recharge_onset(redo_trials,:)
    figure
    for i=1:length(redo_trials)
        subplot(length(redo_trials),1,i)
        hold on
        plot(EEG.times([2000*dt:2100*dt]),EEG.data(electr_oi_i,[2000*dt:2100*dt],redo_trials(i)))
        plot(recharge_onset(redo_trials(i)),[0 ], 'ro')
    end
    old_redo_trials=redo_trials;
    fprintf('Decreasing sensitivity\n')
    count=0;
    while ~isempty(redo_trials) & scale<20
        count=count+1;
        scale=scale+1;
        fprintf('Iteration %d - Scale %d\n',count,scale)
        for trial=1:size(redo_trials,1)
            basemean=mean(EEG.data(electr_oi_i,base_i,redo_trials(trial)));
            basesd=std(EEG.data(electr_oi_i,base_i,redo_trials(trial)));
            found_recharge=[0];
            just_found_one=0;
            for t=1:length(EEG.times)
                if just_found_one>refract*dt
                    just_found_one=0;
                elseif just_found_one >0
                    just_found_one=just_found_one+1;
                end
        %         t=t+1
        %         EEG.times(t)
        %         EEG.data(electr_oi_i,t+1,trial)
                if EEG.times(t)>tms_time_wndw(1) & EEG.times(t)<tms_time_wndw(2) & just_found_one==0
                    if  abs(EEG.data(electr_oi_i,t+1,redo_trials(trial)))> basemean+basesd*scale
                        if any(EEG.data(electr_oi_i,t+1:t+4,redo_trials(trial))~=0) & any(EEG.data(electr_oi_i,t-2:t,redo_trials(trial))~=0)
                            recharge_onset(redo_trials(trial))=EEG.times(t);
                            found_recharge=1;
                            just_found_one=just_found_one+1;
                        end
                    end
                end
            end
        end
        redo_trials=[];
        redo_trials=[redo_trials; find(recharge_onset>rechargemean+rechargesd | recharge_onset<rechargemean-rechargesd )];
       
        clf
        for i=1:length(old_redo_trials)
            subplot(length(old_redo_trials),1,i)
            hold on
            plot(EEG.times([2000*dt:2100*dt]),EEG.data(electr_oi_i,[2000*dt:2100*dt],old_redo_trials(i)))
            plot(recharge_onset(old_redo_trials(i)),[0], 'ro')
        end
        pause(1)

    end
end
if ~isempty(redo_trials)
    fprintf('%d odd trials remaining\nSetting recharge time to mean\d',length(redo_trials))
    recharge_onset(redo_trials,:)=nan;
    rechargemean=nanmean(recharge_onset);
    for i =1:length(redo_trials)
        recharge_onset(redo_trials(i),:)=rechargemean;
    end
    rechargemean=mean(recharge_onset);
    fprintf('\nRecharges where assigned for all trials.\nMean times are:\t%d\t%d\t%d\n',round(rechargemean))
else
    rechargemean=mean(recharge_onset);
    fprintf('\nRecharges where found for all trials.\nMean times are:\t%d\t%d\t%d\n',round(rechargemean))
end
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Median filter Recharges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
medianfilt=true;
zero=false;
interp=false;
count=0;
for trial=1:size(EEG.data,3)
    recharge=[find(EEG.times==recharge_onset(trial))+round(recharge_zero_time(1)*dt): find(EEG.times==recharge_onset(trial))+round(recharge_zero_time(2)*dt)];
     if interp==true
        for electr=1:size(EEG.data,1)
        %     clf
        %     plot(EEG.times(999*dt:1020*dt), EEG.data(electr,[999*dt:1020*dt],1))
                hold on
                y=EEG.data(electr,:,trial);
                x=EEG.times;
                x([recharge])=[];
                y([recharge])=[];
                xx=EEG.times([recharge]);   
                yy=interp1(x,y,xx,'pchip');%pchip is cubib
                EEG.data(electr,[recharge],trial)=yy;
        %     plot(xx, yy,'.')
         end
    elseif zero==true
        EEG.data(:,[recharge],trial)=0;
        %plot(EEG.times(pulse2),zeros(size(pulse1)))
    elseif medianfilt==true
        for electr=1:size(EEG.data,1)
            filt_ord =length(recharge);%number of datapoints considered left and right
            EEG.data(electr,[recharge],trial)=medfilt1(EEG.data(electr,[recharge],trial),filt_ord,'truncate');
        end
    end
    
    if plot_steps==true
        figure(preprocfig)
        if any(trial==trials_to_plot)
            count=count+1;
            subplot(nr_trials_to_plot,1,count)
            hold on
            plot(recharge_onset(trial),EEG.data(electr_oi_i,[find(EEG.times==recharge_onset(trial)),find(EEG.times==recharge_onset(trial)),find(EEG.times==recharge_onset(trial))],trial),'bo')        
            plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))
            ylim([min(EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial)),max(EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))])
        end
    end
    if plot_all_elecs==true
        figure(allelecfig)
        if trial==size(EEG.data,3)
            for elec=1:size(EEG.data,1)
                trial=5;               
                subplot(round(size(EEG.data,1)/4),4,elec)
                set(gca,'XTick',[], 'YTick', [])
                plot(recharge_onset(trial),EEG.data(electr_oi_i,[find(EEG.times==recharge_onset(trial)),find(EEG.times==recharge_onset(trial)),find(EEG.times==recharge_onset(trial))],trial),'bo')        
                plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(elec,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))
                ylim([min(EEG.data(elec,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial)),max(EEG.data(elec,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))])       
                %text(min(xlim)+1,double(min(ylim)+(max(ylim)-min(ylim))/2),EEG.chanlocs(elec).labels)
                set(eleclabel,'Position',[min(xlim)+1,double(min(ylim)+(max(ylim)-min(ylim))/2)])
            end         
        end
    end
       
end
if plot_steps==true
    figure(preprocfig)
    legend('raw','pulseonset','raw-pulse','rechargeonset','raw-pulse-recharge')
end

before_basic_cleaning=EEG;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Cleaning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG=eeg_checkset(EEG);
if plot_steps==true
    preprocfig2=figure;
     for trial=1:nr_trials_to_plot
         subplot(nr_trials_to_plot,1,trial)
         plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trials_to_plot(trial)))
     end
end

%% downsample
EEG.data(size(EEG.data,1)-1:size(EEG.data,1),:,:)=keep_EOG;%EOG needs to be downsampled too
EEG= pop_resample(EEG,1000);
keep_EOG=EEG.data(size(EEG.data,1)-1:size(EEG.data,1),:,:);%take back out
dt=1;
if plot_steps==true
    figure(preprocfig2)
     for trial=1:nr_trials_to_plot
         subplot(nr_trials_to_plot,1,trial)
         hold on
         plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trials_to_plot(trial)))
     end
end

%% filter
EEG=tesa_filtbutter(EEG,1,100,4,'bandpass')

if plot_steps==true
    figure(preprocfig2)
     for trial=1:nr_trials_to_plot
         subplot(nr_trials_to_plot,1,trial)
         hold on
         plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trials_to_plot(trial)))
     end
end


%% remove bad channels
figure; pop_spectopo(EEG, 1, [-2000  1999], 'EEG' , 'percent', 20, 'freq', [6 10 22], 'freqrange',[2 80],'electrodes','on');

if partic==2
    bad=[14 20 29 47 54 57];
end
EEG=pop_select(EEG, 'nochannel',bad);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero out pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%zero out for ICA
%it's downsampled now, so tha'ts different!
zero=true;
interp=false;

for trial=1:size(EEG.data,3)
    pulse=[find(EEG.times==round(pulse_onset(trial)))+round(pulse_zero_time(1)*dt): find(EEG.times==round(pulse_onset(trial)))+round(pulse_zero_time(2)*dt)];   
    if interp==true
        for electr=1:size(EEG.data,1)
        %     clf
        %     plot(EEG.times(999*dt:1020*dt), EEG.data(electr,[999*dt:1020*dt],1))
                hold on
                y=EEG.data(electr,:,trial);
                x=EEG.times;
                x([pulse1])=[];
                y([pulse1])=[];
                xx=EEG.times([pulse]);   
                yy=interp1(x,y,xx,'pchip');%pchip is cubib
                EEG.data(electr,[pulse],trial)=yy;
        %     plot(xx, yy,'.')
         end
    elseif zero==true
        EEG.data(:,[pulse],trial)=0;
    %plot(EEG.times(pulse2),zeros(size(pulse1)))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_steps==true
    figure(preprocfig2)
     for trial=1:nr_trials_to_plot
         subplot(nr_trials_to_plot,1,trial)
         hold on
         plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trials_to_plot(trial)))
         xlim([-5 100])
     end
end
keep_here=EEG;
if TESAICA
    % addpath to be able to run fast ica
    addpath('C:\Users\ckohl\Documents\MATLAB\FastICA_2.5')
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
%     EEG = pop_tesa_compselect( EEG,'comps',15,'figSize','small','plotTimeX',[-200 500],'plotFreqX',[1 100],'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','off','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','off','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','off','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off','elecNoise','off','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    EEG = pop_tesa_compselect( EEG,'compCheck','off','comps',[],'figSize','small','plotTimeX',[-200 500],'plotFreqX',[1 100],'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','on','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','on','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off','elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
else
    
    EEG=pop_runica(EEG, 'icatype', 'runica')
    pop_selectcomps(EEG, [1:35] );
    EEG = pop_subcomp( EEG, [1], 0);
end

if plot_steps==true
    figure(preprocfig2)
     for trial=1:nr_trials_to_plot
         subplot(nr_trials_to_plot,1,trial)
         hold on
         plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trials_to_plot(trial)))
     end
end

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% manual rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%put raw EOG back in
EEG.data(size(EEG.data,1)-1:size(EEG.data,1),:,:)=keep_EOG;
%pop_eegplot( EEG, 1, 1, 1);
if partic==2
    if S1_80
        EEG = pop_rejepoch( EEG, [2 7 10 25 26 28 31 42 46 50 51 55 56 59 60 61 65 69 70 71 75 76 97 104] ,0);
    else%SI100%
        EEG = eeg_eegrej( EEG, [1 12000;28000 32000;52000 64000;72000 76000;100000 112000;124000 128000;136000 140000;160000 168000;176000 188000]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate pulses and Filter again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zero=false;
interp=true;
for trial=1:size(EEG.data,3)
    pulse1=[find(EEG.times==round(pulse_onset(trial)))+round(pulse_zero_time(1)*dt): find(EEG.times==round(pulse_onset(trial)))+round(pulse_zero_time(2)*dt)];    
    if interp==true
        for electr=1:size(EEG.data,1)
        %     clf
        %     plot(EEG.times(999*dt:1020*dt), EEG.data(electr,[999*dt:1020*dt],1))
                hold on
                y=EEG.data(electr,:,trial);
                x=EEG.times;
                x([pulse])=[];
                y([pulse])=[];
                xx=EEG.times([pulse]);   
                yy=interp1(x,y,xx,'pchip');%pchip is cubib
                EEG.data(electr,[pulse],trial)=yy;
        %     plot(xx, yy,'.')
         end
    elseif zero==true
        EEG.data(:,[pulse],trial)=0;
    %plot(EEG.times(pulse2),zeros(size(pulse1)))
    end
end
if refilt
    EEG=tesa_filtbutter(EEG,1,100,4,'bandpass')
        %need new trials to plot
        trials_to_plot=sort(randi(size(EEG.data,3),1,nr_trials_to_plot));
        if plot_steps==true
             for trial=1:nr_trials_to_plot
                 subplot(nr_trials_to_plot,1,trial)
                 hold on
                 plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trials_to_plot(trial)))
             end
        end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate Removed Channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(EEG.data,1) < 65
    EEG=pop_interp(EEG,channels,'spherical');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Re-reference to Average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_reref( EEG, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TESAICA
    name='TESAICA';
else
    name='runICA';
end
if refilt
    name=strcat(name,'_refilt');
else
    name=strcat(name,'_nofilt');
end

if ~S1_80
    nameplus='100_';
else
    nameplus='';
end

%% save for Matlab
cd(raw_path)
save(strcat('Beta0',num2str(partic),'_TEP_',nameplus,name),'EEG')
%% transform for MNE-Python
EEG = pop_saveset( EEG, 'filename',strcat('Beta0',num2str(partic),'_TEP_',nameplus,name,'.set'),'filepath',raw_path);

