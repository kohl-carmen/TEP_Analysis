%% Let's see if TEPs are related to beta somehow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

partic=2;
data_path='F:\Brown\TMS\Pilot\';
Partic=[2,4];
dt=1;
TESAICA=0; %if 0, runICA for blink, if 1 TESAICA automatic (not done here but pick which preprocessed file you want)
electr_oi='C3';


fieldtrip_dir='C:\Users\ckohl\Documents\fieldtrip-20190802\fieldtrip-20190802';
eeglab_dir='C:\Users\ckohl\Documents\MATLAB\eeglab2019_0';
rmpath(fieldtrip_dir)
addpath(eeglab_dir)
eeglab
close all
%% load data
if TESAICA
    name='TESA';
else
    name='run';
end
data=load(strcat(data_path,strcat('Beta0',num2str(partic),'_TEP_1k_',name,'ICA_filt100')));
EEG=data.EEG;

bandfreq=[15 29];


%find electrode
for chan= 1:length(EEG.chanlocs)
    if length(EEG.chanlocs(chan).labels)==2
        if EEG.chanlocs(chan).labels==electr_oi
            electr_oi_i=chan;
        end
    end
end

%% PPT
ppt=1;
h = actxserver('PowerPoint.Application');
Presentation = h.Presentation.Add;


%% Ryan has a beta toolbox, but for now, I'll just do it my own way
rmpath(eeglab_dir)
addpath C:\Users\ckohl\Documents\fieldtrip-20190802\fieldtrip-20190802
ft_defaults
keep=EEG;
%EEG2=eeglab2fieldtrip(EEG,'timelockanalysis','none');
EEG=eeglab2fieldtrip(EEG,'preprocessing','none');

figure
hold on
for trial=1:length(EEG.trial)
    x=EEG.trial{trial}(5,:);
    fs=1000;
    y = fft(x);
    n = length(x);          % number of samples
    f = (0:n-1)*(fs/n);     % frequency range
    power(trial,:) = abs(y).^2/n;    % power of the DFT

    plot(f,power(trial,:))
end
plot(f,mean(power),'k','Linewidth',2)
xlabel('Frequency')
ylabel('Power')
xlim([0 50])

% see if certain electrodes have higer beta power than others
figure
clf
hold on
c=parula(length(EEG.elec.label));
for elec=1:length(EEG.elec.label)
    power=[];
    for trial=1:length(EEG.trial)
        x=EEG.trial{trial}(elec,:);
        fs=1000;
        y = fft(x);
        n = length(x);          % number of samples
        f = (0:n-1)*(fs/n);     % frequency range
        power(trial,:) = abs(y).^2/n;    % power of the DFT
    end
    power=mean(power);
    plot(f,power,'Color',c(elec,:),'Linewidth',2)
    xtext=bandfreq(end)%;randi(bandfreq);
    tx=text(xtext,power(f==xtext),EEG.elec.label{elec});
    tx.Color=c(elec,:);
    tx.FontWeight='bold';
    tx.FontSize=10;
    if EEG.elec.label{elec}(1)=='C' & EEG.elec.label{elec}(2)=='3'
        C3power=power;
    end
        
end
plot(f,C3power,'Color','k','Linewidth',3)
xtext=bandfreq(end)%;randi(bandfreq);
tx=text(xtext,C3power(f==xtext),'C3');
tx.Color='k';
tx.FontWeight='bold';
tx.FontSize=10;
xlabel('Frequency')
ylabel('Power')
xlim([0 50])
xlim(bandfreq)
print('-dpng','-r150',strcat('temp','.png'));
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500


 

%% TFR
beta_only=1; %if beta_only, TFR is looking at bea, if not, it's to have a general look

cfg = [];	                
cfg.method     = 'wavelet';                
cfg.width      = 7; %we want at least 4 apparently
cfg.output     = 'pow';	
if beta_only
    cfg.foi        = bandfreq(1):1:bandfreq(end);%10:1:30;%for short segments, use 3 as starting req for toi and foi we can be pretty generous. nothing to do with the time window and stuff	
    cfg.toi        = [-.75:.01:0];%EEG.time{1};
else
    cfg.foi        = 10:1:30;%for short segments, use 3 as starting req for toi and foi we can be pretty generous. nothing to do with the time window and stuff	
    cfg.toi        = EEG.time{1};
end
cfg.keeptrials = 'yes';
cfg.channel    = 'C3';
TFR = ft_freqanalysis(cfg, EEG);
cfg.keeptrials = 'no';
TFR_avg=ft_freqanalysis(cfg, EEG);

%baseline
% this is a bit tricky. Lots of stuff happens at zero so I don't want that
% to be in my baseline (and it smears into before zero!). But I care about
% the entire period up to zero equally, so I can't really define a clear
% base time. So I'll probably jsut take the very first time point and do a
% relative baseline?
% if beta_only
%     cfg=[];
%     cfg.baseline=
%     cfg.parameter='powspctrm';
%     TFR=fr_freqbaseline
% end
% 
%        
%% plot avg (to see max freq
figure('units','normalized','outerposition', [0 0 1 1]);
cfg=[];      
    cfg.parameter='powspctrm';
    cfg.colormap=jet;
    cfg.colorbar='yes';
    cfg.channel='C3';
%         cfg.baseline=[-1 0];
%         cfg.baselinetype ='relative';
    cfg.title = strcat('Avg');
ft_singleplotTFR(cfg,TFR_avg)
print('-dpng','-r150',strcat('temp','.png'));
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500
clear TFR_avg


%% plot trials
figure('units','normalized','outerposition', [0 0 1 1]);
count=0;
for trial=1:length(EEG.trial)
    count=count+1;
    if count>9
        print('-dpng','-r150',strcat('temp','.png'));
        blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
        Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
        Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500

        count=1;
        figure('units','normalized','outerposition', [0 0 1 1]);
    end
    
    subplot(3,3,count)
    cfg=[];      
        cfg.parameter='powspctrm';
        cfg.colormap=jet;
        cfg.colorbar='yes';
        cfg.channel='C3';
%         cfg.baseline=[-1 0];
%         cfg.baselinetype ='relative';
        cfg.trials= [trial];
        cfg.title = strcat('Trial',num2str(trial));
     ft_singleplotTFR(cfg,TFR)
     hold on
     plot([TFR.time(1), TFR.time(end)],[bandfreq(1), bandfreq(1)],'k','Linewidth',2)
     plot([TFR.time(1), TFR.time(end)],[bandfreq(2), bandfreq(2)],'k','Linewidth',2)
end
print('-dpng','-r150',strcat('temp','.png'));
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500



%% Let's only  look at peak frequency
if partic==2
    freq_oi=[18:22];
end

freq_oi_i=[find(TFR.freq==freq_oi(1)) : find(TFR.freq==freq_oi(end))];
data=squeeze(mean(TFR.powspctrm(:,:,freq_oi_i,:),3));
%delete nans (there actually shouldnt be any)
time=TFR.time;
del_i=isnan(data(1,:));
time(del_i)=[];
data(:,del_i)=[];

%do my own baseline
for trial=1:length(EEG.trial)
    data(trial,:)=data(trial,:)./mean(data(trial,:));
end

%find overall median
datacont=reshape(data,1,size(data,1)*size(data,2));
med=median(datacont);

%find events
cutoff=4*med;
count=0;
trialsoi=[];
for trial=1:length(EEG.trial)
    if any(data(trial,:)>cutoff) % (if median per trial):median(data(trial,:))*6)
        count=count+1;
        trialsoi=[trialsoi,trial];
    end

end

% see if there's one or more and find peaks
event_count=[];
peak=[];
for trial = trialsoi
    overcut=find(data(trial,:)>cutoff);
    event_count(trial)=1;
    overstart=1;%overcut(1);
    for i=1:length(overcut)-1
        if overcut(i+1)~=overcut(i)+1       
            %if theres a new period of high beta startin here, find the peak for the
            %previous one 
            [temp, temp_i]=max(data(trial,overcut(overstart:i)));
            temp_i_i=overcut(overstart:i);
            peak(trial,event_count(trial))=temp_i_i(temp_i);
            overstart=i+1;%overcut(i+1);
            event_count(trial)=event_count(trial)+1;
        end
    end
    [temp, temp_i]=max(data(trial,overcut(overstart:end)));
    temp_i_i=overcut(overstart:end);
    peak(trial,event_count(trial))=temp_i_i(temp_i);
end
peak(peak==0)=nan;
         

for trial=trialsoi(1:4)
    figure('units','normalized','outerposition', [0 0 1 1]);

    subplot(3,1,1)
    title(strcat(num2str(event_count(trial)),' event(s) found'))
    hold on
    plot(time, data(trial,:),'b','Linewidth',2)
    plot([time(1),time(end)],[cutoff, cutoff],'k--','Linewidth',2)
    xlim([TFR.time(1), TFR.time(end)+0.05])
    legtext={'Power','Event Cutoff'};
    %plot peaks
    yl=ylim();
    for p=1:size(peak,2)
        if ~isnan(peak(trial,p))
            plot([TFR.time(peak(trial,p)),TFR.time(peak(trial,p))],yl,'r--','Linewidth',2)
            legtext{end+1}='Beta Peak';
        end
    end
    legend(legtext)
    
    
    subplot(3,1,[2:3])
 
     cfg=[];
        
        %cfg.layout=layout;
        cfg.parameter='powspctrm';
        cfg.colormap=jet;
        cfg.colorbar='yes';
        cfg.channel='C3';
         cfg.baseline=[-1 0];
         cfg.baselinetype ='relative';
%        cfg.zlim=[0 200]
       cfg.trials= [trial];
       cfg.title = strcat('Trial',num2str(trial));
       %Grand_base=Grand_no_trials{1};
       %Grand_base.powspctrm=Grand_no_trials{1}.powspctrm-repmat(mean(Base_S{1},3),1,1,size(Grand_no_trials{1}.powspctrm,3));
%  figure
 ft_singleplotTFR(cfg,TFR)
 hold on
 plot([TFR.time(1), TFR.time(end)],[bandfreq(1), bandfreq(1)],'k--','Linewidth',2)
 plot([TFR.time(1), TFR.time(end)],[bandfreq(2), bandfreq(2)],'k--','Linewidth',2)
 plot([TFR.time(1), TFR.time(end)],[freq_oi(1), freq_oi(1)],'k-','Linewidth',2)
 plot([TFR.time(1), TFR.time(end)],[freq_oi(end), freq_oi(end)],'k-','Linewidth',2)

end


%find beta peaks in data
for trial=trialsoi
    figure
    eeg=EEG.trial{trial}(electr_oi_i,:);
    plot(EEG.time{1},eeg,'k')
    hold on
    yl=ylim();
    for p=1:size(peak,2)
        if ~isnan(peak(trial,p))
            plot([TFR.time(peak(trial,p)),TFR.time(peak(trial,p))],yl,'r--','Linewidth',2)
            legtext{end+1}='Beta Peak';
        end
    end
    xlim([TFR.time(1), TFR.time(end)])

end




%% ok time for ryan
% x ia the data -> Time-by-trial matrix of timeseries trials for detection/non-detection prestimulus MEG
x=[];
for trial=1:length(EEG.trial)
     x(:,trial)=EEG.trial{trial}(electr_oi_i,1:find(EEG.time{1}==0))';
%     x(:,trial)=EEG.trial{trial}(electr_oi_i,find(EEG.time{1}==-0.75):find(EEG.time{1}==0))';
end
classLabels=ones(length(EEG.trial),1);
eventBand = [15 29]; %Frequency range of spectral events
fVec = 15:29; %Vector of fequency values over which to calculate TFR
Fs = 1000; %Sampling rate of time-series
findMethod = 2; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = true; %Generate standard visualization plots for event features across all subjects/sessions
tVec = (-1/Fs:1/Fs:0);
[specEvents,TFRs,timeseries] = spectralevents(eventBand,fVec,Fs,findMethod,vis,x,classLabels); %Run spectral event analysis

% 
% [TFR,tVec,fVec] = spectralevents_ts2tfr(x(:,1),fVec,Fs,7)
% specEvents = spectralevents_find(1, eventBand, 6, tVec, fVec, TFR, 1)
% spectralevents_vis(specEvents,x(:,1),TFR, tVec,fVec)
% specEvents.Events.Events

%%  plot Ryans TFRs to sanity check
newpow=[];
for i=1:72
    for f=1:length(fVec)
        newpow(i,1,f,:)=TFRs{1}(f,:,i);
    end
end
       
TFR_temp=TFR;
TFR_temp.powspctrm=newpow;
TFR_temp.time=[-1000:0];
TFR_temp.freq=fVec;

figure('units','normalized','outerposition', [0 0 1 1]);
count=0;
for trial=1:length(EEG.trial)
    count=count+1;
    if count>9
        print('-dpng','-r150',strcat('temp','.png'));
        blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
        Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
        Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500

        count=1;
        figure('units','normalized','outerposition', [0 0 1 1]);
    end
    
    subplot(3,3,count)
    cfg=[];      
        cfg.parameter='powspctrm';
        cfg.colormap=jet;
        cfg.colorbar='yes';
        cfg.channel='C3';
%         cfg.baseline=[-1 0];
%         cfg.baselinetype ='relative';
        cfg.trials= [trial];
        cfg.title = strcat('Trial',num2str(trial));
     ft_singleplotTFR(cfg,TFR_temp)
     hold on
%      plot([TFR.time(1), TFR.time(end)],[bandfreq(1), bandfreq(1)],'k','Linewidth',2)
%      plot([TFR.time(1), TFR.time(end)],[bandfreq(2), bandfreq(2)],'k','Linewidth',2)
end
print('-dpng','-r150',strcat('temp','.png'));
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500


%let's plot the events because weird
%so it seems like it just finds events at all frequencies?? not just freq
%oi given?

event_trials=[];
event_max_times=[];
event_onsets=[];
event_offsets=[];
ryan_assumed_time=0:1/Fs:1;
ryan_assumed_time=round(ryan_assumed_time,3);
actual_time=EEG.time{1}( find(EEG.time{1}==tVec(1)/dt):find(EEG.time{1}==tVec(end)/dt));
for event=1:length(specEvents.Events.Events.trialind)
    %only those events that are in beta freq
    if specEvents.Events.Events.lowerboundFspan(event)>=eventBand(1) & specEvents.Events.Events.upperboundFspan(event)<= eventBand(end)
        event_trials=[event_trials,specEvents.Events.Events.trialind(event)];
        %have to convert each time to the actual time (ryan doesn't seem to
        %care about timeing?
        maxtime=specEvents.Events.Events.maximatiming(event);
        maxtime=actual_time(find(ryan_assumed_time==round(maxtime,3)));
        event_max_times=[event_max_times, maxtime];
        
        onsettime=specEvents.Events.Events.onsettiming(event);
        onsettime=actual_time(find(ryan_assumed_time==round(onsettime,3)));
        event_onsets=[event_onsets,onsettime];
        
        offsettime=specEvents.Events.Events.offsettiming(event);
        offsettime=actual_time(find(ryan_assumed_time==round(offsettime,3)));
        event_offsets=[event_offsets,offsettime];
        
    end
end


%% plot ryans events my way
newpow=[];
for i=1:72
    for f=1:length(fVec)
        newpow(i,1,f,:)=TFRs{1}(f,:,i);
    end
end
     
TFR_temp=TFR;
TFR_temp.powspctrm=newpow;
TFR_temp.time=[-1000:0];
TFR_temp.freq=fVec;

event_counter=0;
for trial=unique(event_trials)
    figure('units','normalized','outerposition', [0 0 1 1]);

    subplot(3,1,1)    
    title(strcat(num2str(sum(event_trials==trial)),' event(s) found'))
    hold on
    plot(time, data(trial,:),'b','Linewidth',2)
    plot(actual_time,timeseries{1}(:,trial)','k--','Linewidth',2)
    %xlim([TFR.time(1), TFR.time(end)+0.05])
    
    yl=ylim();
    %mark events
    for event=1:sum(event_trials==trial)
        event_counter=event_counter+1;
        plot([event_max_times(event_counter),event_max_times(event_counter)],yl,'r--','Linewidth',2)
%         plot([event_onsets(event_counter),event_onsets(event_counter)],yl,'r--','Linewidth',1)
%         plot([event_offsets(event_counter),event_offsets(event_counter)],yl,'r--','Linewidth',1)
    end
    legtext={'Power','Data'};
    
    
    subplot(3,1,[2:3])
 
     cfg=[];
        
        %cfg.layout=layout;
        cfg.parameter='powspctrm';
        cfg.colormap=jet;
        cfg.colorbar='yes';
        cfg.channel='C3';
         cfg.baseline=[-1 0];
         cfg.baselinetype ='relative';
%        cfg.zlim=[0 200]
       cfg.trials= [trial];
       cfg.title = strcat('Trial',num2str(trial));
       %Grand_base=Grand_no_trials{1};
       %Grand_base.powspctrm=Grand_no_trials{1}.powspctrm-repmat(mean(Base_S{1},3),1,1,size(Grand_no_trials{1}.powspctrm,3));
%  figure
 ft_singleplotTFR(cfg,TFR_temp)
 hold on
 plot([TFR_temp.time(1), TFR_temp.time(end)],[bandfreq(1), bandfreq(1)],'k--','Linewidth',2)
 plot([TFR_temp.time(1), TFR_temp.time(end)],[bandfreq(2), bandfreq(2)],'k--','Linewidth',2)
end



