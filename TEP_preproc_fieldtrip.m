http://www.fieldtriptoolbox.org/tutorial/continuous/
raw_path='C:\Users\ckohl\Desktop\Current\TMS\CurrentData\EEG';

addpath C:\Users\ckohl\Documents\MATLAB\fieldtrip-20200914
ft_defaults

%load
cfg = [];
% cfg.dataset = strcat(raw_path,'\',filename);
cfg.dataset = 'C:\Users\ckohl\Desktop\Current\TMS\CurrentData\EEG\Beta02_TEP_TESAICA_refilt.set'
data_eeg=ft_preprocessing(cfg)

%% DEFINE TRIAL MATRIX
trigger ={'S 13'};
cfg = [];
cfg.dataset                 = strcat(raw_path,'\Beta02_raw_cut.set')
cfg.continuous              = 'yes';
cfg.trialdef.prestim        = .5;         % prior to event onset
cfg.trialdef.poststim       = 1.5;        % after event onset
cfg.trialdef.eventtype      = 'Stimulus'; % see above
cfg.trialdef.eventvalue     = trigger ;
cfg = ft_definetrial(cfg);                % make the trial definition matrix

trl = cfg.trl;

%% LOAD AND REREF
cfg.channel = {'all' '-TP10' '-EOG_aux1_vert' '-EOG_aux2_horz'}; % indicate the channels we would like to read and/or exclude.
cfg.reref = 'yes';        % We want to rereference our data
cfg.refchannel = {'all'}; % Here we specify our reference channels
cfg.implicitref = 'TP10';    % Here we can specify the name of the implicit reference electrode used during the acquisition
data_tms_raw = ft_preprocessing(cfg);

cd('C:\Users\ckohl\Desktop\')
save('data_tms_raw','data_tms_raw','-v7.3');

%% BROWSE
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 -0.001];
ft_databrowser(cfg, data_tms_raw);

%% ERP
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 -0.001];

data_tms_avg = ft_timelockanalysis(cfg, data_tms_raw);

clear data_tms_raw;


% plot all channels
close all

for i=1:numel(data_tms_avg.label)                   % Loop through all channels
    figure;
    plot(data_tms_avg.time, data_tms_avg.avg(i,:)); % Plot this channel versus time
    xlim([-0.1 0.6]);     % Here we can specify the limits of what to plot on the x-axis
    ylim([-23 15]);       % Here we can specify the limits of what to plot on the y-axis
    title(['Channel ' data_tms_avg.label{i}]);
    ylabel('Amplitude (uV)')
    xlabel('Time (s)');
end

% closer look at C3

channel = 'C3';

figure;
i = find(strcmp(channel, data_tms_avg.label));
plot(data_tms_avg.time, data_tms_avg.avg(i,:));   % Plot data
xlim([-0.1 0.6]);    % Here we can specify the limits of what to plot on the x-axis
ylim([-23 15]);      % Here we can specify the limits of what to plot on the y-axis
title(['Channel ' data_tms_avg.label{i}]);
ylabel('Amplitude (uV)')
xlabel('Time (s)');

xlim([-0 0.020]);
ylim([-60 100]);

figure;
channel_idx = find(strcmp(channel, data_tms_avg.label));
plot(data_tms_avg.time, data_tms_avg.avg(channel_idx,:));  % Plot all data
xlim([-0.1 0.6]);    % Here we can specify the limits of what to plot on the x-axis
ylim([-60 100]);     % Here we can specify the limits of what to plot on the y-axis
title(['Channel ' data_tms_avg.label{channel_idx}]);
ylabel('Amplitude (uV)')
xlabel('Time (s)');

hold on; % Plotting new data does not remove old plot

% Specify time-ranges to higlight
ringing  = [-0.0002 0.0044];
muscle   = [ 0.0044 0.015 ];
decay    = [ 0.015  0.200 ];
recharge = [ 0.4994 0.5112];

colors = 'rgcm';
labels = {'ringing','muscle','decay','recharge'};
artifacts = [ringing; muscle; decay; recharge];

for i=1:numel(labels)
  highlight_idx = [nearest(data_tms_avg.time,artifacts(i,1)) nearest(data_tms_avg.time,artifacts(i,2)) ];
  plot(data_tms_avg.time(highlight_idx(1):highlight_idx(2)), data_tms_avg.avg(channel_idx,highlight_idx(1):highlight_idx(2)),colors(i));
end
legend(['raw data', labels]);

%% REJECT ARTIFACTS
% Ringing/Step Response artifact
cfg                         = [];
cfg.method                  = 'marker'; % The alternative is 'detect' to detect the onset of pulses
cfg.dataset                 = strcat(raw_path,'\Beta02_raw_cut.set')
cfg.prestim                 = .001;     % First time-point of range to exclude
cfg.poststim                = .006;     % Last time-point of range to exclude
cfg.trialdef.eventtype      = 'Stimulus';
cfg.trialdef.eventvalue     = trigger ;
cfg_ringing = ft_artifact_tms(cfg);     % Detect TMS artifacts

% Here we use a negative value because the recharging artifact starts AFTER TMS-pulse onset
cfg.prestim   = -.011;
cfg.poststim  = .013;
cfg_recharge  = ft_artifact_tms(cfg); % Detect TMS artifacts

% Combine into one structure
cfg_artifact = [];
cfg.dataset  = strcat(raw_path,'\Beta02_raw_cut.set');
cfg_artifact.artfctdef.ringing.artifact = cfg_ringing.artfctdef.tms.artifact; % Add ringing/step response artifact definition
cfg_artifact.artfctdef.recharge.artifact   = cfg_recharge.artfctdef.tms.artifact; % Add recharge artifact definition

cfg_artifact.artfctdef.reject = 'partial'; % Can also be 'complete', or 'nan';
cfg_artifact.trl = trl; % We supply ft_rejectartifact with the original trial structure so it knows where to look for artifacts.
cfg_artifact.artfctdef.minaccepttim = 0.0001; % This specifies the minimumm size of resulting trials. You have to set this, the default is too large for thre present data, resulting in small artifact-free segments being rejected as well.
cfg = ft_rejectartifact(cfg_artifact,data_tms_raw); % Reject trials partially

data_tms_segmented=cfg;

%% BROWSE
% Browse through the segmented data
cfg = [];
cfg.artfctdef = cfg_artifact.artfctdef; % Store previously obtained artifact definition
cfg.continuous = 'yes'; % Setting this to yes forces ft_databrowser to represent our segmented data as one continuous signal
ft_databrowser(cfg, data_tms_segmented);

%% ICA
%% Perform ICA on segmented data
cfg = [];
cfg.demean = 'yes';
cfg.method = 'fastica';        % FieldTrip supports multiple ways to perform ICA, 'fastica' is one of them.
cfg.fastica.approach = 'symm'; % All components will be estimated simultaneously.
cfg.fastica.g = 'gauss';

comp_tms = ft_componentanalysis(cfg, data_tms_segmented);

save('comp_tms','comp_tms','-v7.3');


cfg = [];
comp_tms_avg = ft_timelockanalysis(cfg, comp_tms);


figure;
cfg = [];
cfg.viewmode = 'butterfly';
cfg.channel   ='fastica011';
ft_databrowser(cfg, comp_tms_avg);

cfg=[];
cfg.elec=data_tms_raw.elec;
layout=ft_prepare_layout(cfg);

figure;
cfg           = [];
cfg.component = [1:63];
cfg.comment   = 'no';
cfg.layout    =  layout;% 'easycapM10'; % If you use a function that requires plotting of topographical information you need to supply the function with the location of your channels
ft_topoplotIC(cfg, comp_tms);


cfg          = [];
cfg.demean   = 'no'; % This has to be explicitly stated as the default is to demean.
cfg.unmixing = comp_tms.unmixing; % Supply the matrix necessay to 'unmix' the channel-series data into components
cfg.topolabel = comp_tms.topolabel; % Supply the original channel label information

comp_tms          = ft_componentanalysis(cfg, data_tms_segmented);

%%  REMOVE COMPONENTS
% skipping this for not
data_tms_clean_segmented=data_tms_segmented;

cfg                = [];
cfg.vartrllength   = 2;
cfg.preproc.demean = 'no'; % Demeaning is still applied on the segments of the trials, rather than the entire trial. To avoid offsets within trials, set this to 'no'

data_tms_clean_avg = ft_timelockanalysis(cfg, data_tms_clean_segmented);

% Plot all channels
for i=1:numel(data_tms_clean_avg.label) % Loop through all channels
    figure;
    plot(data_tms_clean_avg.time, data_tms_clean_avg.avg(i,:),'b'); % Plot all data
    xlim([-0.1 0.6]); % Here we can specify the limits of what to plot on the x-axis
    title(['Channel ' data_tms_clean_avg.label{i}]);
    ylabel('Amplitude (uV)')
    xlabel('Time (s)');
end

%% INTERPOLATE
% Apply original structure to segmented data, gaps will be filled with nans
cfg     = [];
cfg.trl = trl;
data_tms_clean = ft_redefinetrial(cfg, data_tms_clean_segmented); % Restructure cleaned data


% Replacing muscle artifact with nans
muscle_window = [0.006 0.013]; % The window we would like to replace with nans
muscle_window_idx = [nearest(data_tms_clean.time{1},muscle_window(1)) nearest(data_tms_clean.time{1},muscle_window(2))]; % Find the indices in the time vector corresponding to our window of interest
for i=1:numel(data_tms_clean.trial) % Loop through all trials
  data_tms_clean.trial{i}(:,muscle_window_idx(1):muscle_window_idx(2))=nan; % Replace the segment of data corresponding to our window of interest with nans
end

% Interpolate nans using cubic interpolation
cfg = [];
cfg.method = 'pchip'; % Here you can specify any method that is supported by interp1: 'nearest','linear','spline','pchip','cubic','v5cubic'
cfg.prewindow = 0.01; % Window prior to segment to use data points for interpolation
cfg.postwindow = 0.01; % Window after segment to use data points for interpolation
data_tms_clean = ft_interpolatenan(cfg, data_tms_clean); % Clean data

% compute the TEP on the cleaned data
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.05 -0.001];

data_tms_clean_avg = ft_timelockanalysis(cfg, data_tms_clean);

%compare to raw
for i=1:numel(data_tms_avg.label) % Loop through all channels
    figure;
    plot(data_tms_avg.time, data_tms_avg.avg(i,:),'r'); % Plot all data
    hold on;
    plot(data_tms_clean_avg.time, data_tms_clean_avg.avg(i,:),'b'); % Plot all data
    xlim([-0.1 0.6]); % Here we can specify the limits of what to plot on the x-axis
    ylim([-23 15]); % Here we can specify the limits of what to plot on the y-axis
    title(['Channel ' data_tms_avg.label{i}]);
    ylabel('Amplitude (uV)')
    xlabel('Time (s)');
    legend({'Raw' 'Cleaned'});
end

%% DOWNSAMPLE
cfg = [];
cfg.resamplefs = 1000;
cfg.detrend = 'no';
cfg.demean = 'yes'; % Prior to downsampling a lowpass filter is applied, demeaning avoids artifacts at the edges of your trial
data_tms_clean = ft_resampledata(cfg, data_tms_clean);

save('data_tms_clean','data_tms_clean','-v7.3');

%% FILTER
data=data_tms_clean;
cfg                = [];
cfg.hpfilter       = 'yes';        % enable high-pass filtering
cfg.lpfilter       = 'yes';        % enable low-pass filtering
cfg.hpfreq         = 1;           % set up the frequency for high-pass filter
cfg.lpfreq         = 80;          % set up the frequency for low-pass filter
cfg.dftfilter      = 'no';        % enable notch filtering to eliminate power line noise
% cfg.dftfreq        = [50 100 150]; % set up the frequencies for notch filtering
cfg.baselinewindow = [-0.2 -0.01];    % define the baseline window
data               = ft_preprocessing(cfg,data);

%% REJECT
cfg          = [];
cfg.method   = 'trial';
cfg.latency = [-100 400]/1000;
data        = ft_rejectvisual(cfg,data);
% channel  Iz FT7
% trial 1 2 4 7 8 9 11 16 17 18 19 21 26 29 32 40 43 45 47 51 54 56 57 62 63 66 70 71 72 75 78 79 84 88 90 98 102 105 109
cfg        = [];
cfg.metric = 'zvalue';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
data       = ft_rejectvisual(cfg,data);
% the following channels were removed: Fp1, Fp2, AF8, Fpz
% the following trials were removed: 12, 39, 42
data_after_rej=data;
save('data_after_rej','data_after_rej','-v7.3');

%% INTERPOLATE CHANNELS
cfg=[];
ft_channelselection('all')
cfg=[];
cfg.method='spline';
cfg.elec=data.elec;
cfg.badchannel={'Iz','FT7','Fp1','Fp2','AF8','Fpz'};
data_new=ft_channelrepair(cfg,data);
data_after_rej=data_new;
save('data_after_rej','data_after_rej','-v7.3');


%% SAVE
%saving this as an appropriate format is near impossible so I'm putting it 
% into eeglab by hand. Which is just super. I'm thrilled.

%get eeglab
eeglab_dir='C:\Users\ckohl\Documents\MATLAB\eeglab2020_0';
cd(eeglab_dir)
eeglab

%load something with the right channels and stuff
load('C:\Users\ckohl\Desktop\Current\TMS\CurrentData\EEG\Beta02_TEP_100_TESAICA_nofilt.mat')
eeglab redraw

%edit EEG structure by hand
EEG.comments='made in TEP_preproc_fieldtrip but based on Beta02_TEP_100_TESAICA_nofilt';
EEG.nbchan=length(data_after_rej.label);
EEG.trials=length(data_after_rej.trial);
EEG.pnts=length(data_after_rej.trial{1});
EEG.srate=data_after_rej.fsample;
EEG.xmin=data_after_rej.time{1}(1);
EEG.max=data_after_rej.time{1}(end);
EEG.times=data_after_rej.time{1}*1000;
% REORDER CHANNELS
%when channels are interpolated, they get stuck in the end, so to match
%this to what eeglab expects, we need to shuffle a bit
DATA=[];
if EEG.nbchan==length(data_after_rej.label)
    for og_chan=1:EEG.nbchan
        og_label=EEG.chanlocs(og_chan).labels;
        match=[];
        ft_chan=0;
        while isempty(match) & ft_chan<=length(data_after_rej.label)
            ft_chan=ft_chan+1;
            if og_label(1:2)==data_after_rej.label{ft_chan}(1:2)
                if og_label==data_after_rej.label{ft_chan}
                    match=ft_chan;
                end
            end
        end
        if isempty(match)
            fprintf('No match found for %s\n',og_label)
            brk
        end
        Sort_i(og_chan)=match;
        for trial=1:length(data_after_rej.trial)
            DATA(og_chan,:,trial)=data_after_rej.trial{trial}(match,:);
        end
    end
else
    fprinft('Numbers of trials do not match \n')
end
EEG.data=DATA;
EEG.icachansind=[];
EEG.chanlocs=EEG.chanlocs(1:63);
% EEG.event=[];
EEG.epoch=[];
EEG.reject=[];
EEG.etc=[];
EEG.icaBadComp1=[];

%event latencies
% lat=501;
% for trial=2:length(data_after_rej.trial)
%     lat(trial,:)=lat(trial-1,:)+2000;
% end
EEG=eeg_checkset(EEG);
eeglab redraw

%% transform for MNE-Python
EEG = pop_saveset( EEG, 'filename',strcat('Beta02_fieldtrip.set'),'filepath',raw_path);






