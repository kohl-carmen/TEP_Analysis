%% TEP Preprocessing
% % most basic preprocessing
% - interpolate pulse
% - median filter recharge (if this wasn't cut out with the pulse)
% - downsample
% - filter
% - reject bad channels
% - reject bad trials
% - interpolate bad channels
% - avg ref


%% INIT

clear
close all

partic=2; %[2,4]
S1_80=1;

% load eeglab
eeglab_dir='C:\Users\ckohl\Documents\MATLAB\eeglab2020_0';
cd(eeglab_dir)
eeglab

Partic=[2,4];
dt=25;
electr_oi='C3';

pulse_zero_time=[-2,17];%[-.1,2];
recharge_zero_time=[-1 3];
plot_times=[-50 200];



%% LOAD RAW
%% TMS session from Sep11, Beta02, single pulses
raw_path='C:\Users\ckohl\Desktop\Current\Data\TMS\EEG\';

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

EEG = pop_loadbv(strcat(raw_path,'Raw - EEG-20200912T162305Z-001'), filename, [], []);

channels=EEG.chanlocs;%for later
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw


%find electrode oi
for chan= 1:length(EEG.chanlocs)
    if length(EEG.chanlocs(chan).labels)==2
        if EEG.chanlocs(chan).labels==electr_oi
            electr_oi_i=chan;
        end
    end
end

%% epoch
EEG = pop_epoch( EEG, { 'S 13'}, [-.5 .5], 'epochinfo', 'yes');
EEG = pop_rmbase( EEG, [-200   -10]);
keep_EOG=EEG.data(64:65,:,:);


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
interp_type='linear';%'pchip');
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
                yy=interp1(x,y,xx,interp_type);%pchip is cubib
                EEG.data(electr,[pulse],trial)=yy;
        %     plot(xx, yy,'.')
         end
    elseif zero==true
        EEG.data(:,[pulse],trial)=0;
    %plot(EEG.times(pulse2),zeros(size(pulse1)))
    end
end


if pulse_zero_time(2)<11
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

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Cleaning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG=eeg_checkset(EEG);

%% downsample
EEG.data(size(EEG.data,1)-1:size(EEG.data,1),:,:)=keep_EOG;%EOG needs to be downsampled too
EEG= pop_resample(EEG,1000);

keep_EOG=EEG.data(size(EEG.data,1)-1:size(EEG.data,1),:,:);%take back out
dt=1;


%% filter
EEG=tesa_filtbutter(EEG,1,80,4,'bandpass')


%% remove bad channels
% figure; pop_spectopo(EEG, 1, [-2000  1999], 'EEG' , 'percent', 20, 'freq', [6 10 22], 'freqrange',[2 80],'electrodes','on');

if partic==2
    bad=[14 20 29 47 54 57];
end
EEG=pop_select(EEG, 'nochannel',bad);


%% manual rejection
%put raw EOG back in
EEG.data(size(EEG.data,1)-1:size(EEG.data,1),:,:)=keep_EOG;
%pop_eegplot( EEG, 1, 1, 1);
if partic==2
    if S1_80
       EEG = pop_rejepoch( EEG, [1 2 4:2:10 11 12 16 18 19 20 24 26 27 29 32 38 43 47 51 52 54 56 57 60 61 62 65 66 70 71 72 76 77 78 95 98 99 102:3:105] ,0);
    else%SI100%
        EEG = pop_rejepoch( EEG, [1:3 6 8 9 14 15:17 19 26 27 28 33 35 41 42 45 46:47] ,0);
    end
end

%% Interpolate Removed Channels
if size(EEG.data,1) < 65
    EEG=pop_interp(EEG,channels,'spherical');
end

%% Re-reference to Average
% don't really need EOG now
EEG=pop_select(EEG, 'nochannel',[64, 65]);
EEG = pop_reref( EEG, [] );%EEG = pop_reref( EEG, [],'exclude',[64 65] );

%% check
% eeglab redraw
figure; pop_timtopo(EEG, [-100  300], [20   50  150], ' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save for Matlab
% cd(raw_path)
%  save(strcat('Beta0',num2str(partic),'_TEP_basic_rmv17'),'EEG')
% %% transform for MNE-Python
% EEG = pop_saveset( EEG, 'filename',strcat('Beta0',num2str(partic),'_TEP_basic_rmv17.set'),'filepath',raw_path);

close all


