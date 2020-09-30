%% Let's see if TEPs are related to beta somehow
% I think this is an updated version of TEP_beta but not sure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

partic=2;
data_path='C:\Users\ckohl\Desktop\Current\Data\TMS\EEG\';
Partic=[2,4];
dt=1;
TESAICA=0; %if 0, runICA for blink, if 1 TESAICA automatic (not done here but pick which preprocessed file you want)
electr_oi='C3';


%% PPT
ppt=1;
h = actxserver('PowerPoint.Application');
Presentation = h.Presentation.Add;


eeglab_dir='C:\Users\ckohl\Documents\MATLAB\eeglab2020_0';
cd(eeglab_dir)
eeglab
close all
%% load data
% if TESAICA
%     name='TESA';
% else
%     name='run';
% end
% data=load(strcat(data_path,strcat('Beta0',num2str(partic),'_TEP_1k_',name,'ICA_filt100')));
data=load(strcat(data_path,strcat('Beta0',num2str(partic),'_TEP_basic_forbeta')))
EEG=data.EEG;
% EEG = pop_epoch( EEG, { 'S  1' }, [-1 0], 'epochinfo', 'yes');
EEG = pop_epoch( EEG, { 'S 13' }, [-1 0], 'epochinfo', 'yes');

bandfreq=[15 29];

%find electrode
for chan= 1:length(EEG.chanlocs)
    if length(EEG.chanlocs(chan).labels)==2
        if EEG.chanlocs(chan).labels==electr_oi
            electr_oi_i=chan;
        end
    end
end


eventBand=bandfreq;
fVec=1:40;
Fs=1000;
findMethod=2;
vis=0;
X{1}=squeeze(EEG.data(electr_oi_i,:,:));
classLabels{1}=1;
tVec_assumed=linspace(1/Fs,1,Fs);
[specEv_struct,TFRs,X] = spectralevents(eventBand,fVec,Fs,findMethod,vis,X,classLabels);

event_trial=[];
event_max=[];
event_onset=[];
event_offset=[];

sub_count=0;       
figure('units','normalized','outerposition', [0 0 1 1]);
for trial=unique(specEv_struct.Events.Events.trialind)'  
    sub_count=sub_count+1;
    if sub_count==1
       clf
       hold on
    end
    subplot(2,2,sub_count)

        
% %         subplot(2,1,1)
%             %% ryan time
%             imagesc([tVec(1) tVec(end)],eventBand,TFRs{1}(eventBand(1):eventBand(end),:,trial))
% %             imagesc([tVec(1) tVec(end)],[fVec(1) fVec(end)],TFRs{1}(:,:,trial))
%             colormap jet
%             cb = colorbar;         
%             % Overlay locations of event peaks and the waveform corresponding with each trial
%             hold on
%             max_t=specEv_struct.Events.Events.maximatiming(specEv_struct.Events.Events.trialind==trial);
%             max_f=specEv_struct.Events.Events.maximafreq(specEv_struct.Events.Events.trialind==trial);
%             plot(max_t,max_f,'w.','Markersize',30) %Add points at event maxima
%             
%               yyaxis right
%             plot(tVec,X{1}(:,trial),'w')
%             plot(EEG.times,X{1}(:,trial),'w')
            
            %% actual time
        imagesc([EEG.times(1) EEG.times(end)],eventBand,TFRs{1}(eventBand(1):eventBand(end),:,trial))
%             imagesc([EEG.times(1) EEG.times(end)]],[fVec(1) fVec(end)],TFRs{1}(:,:,trial))
        colormap jet
        cb = colorbar;         
        % Overlay locations of event peaks and the waveform corresponding with each trial
        hold on
        max_t=specEv_struct.Events.Events.maximatiming(specEv_struct.Events.Events.trialind==trial);
        max_f=specEv_struct.Events.Events.maximafreq(specEv_struct.Events.Events.trialind==trial);
        max_t_realtime=[];
        for i=1:length(max_t)
            max_t_realtime(i)=EEG.times(find(round(tVec_assumed,3)==round(max_t(i),3)));
        end
         plot(max_t_realtime,max_f,'w.','Markersize',30) %Add points at event maxima
        
             %% for later - only keep the event closest to time 0
             event_trial=[event_trial,trial];
             event_max=[event_max,max(max_t_realtime)];
             
             
             onset=specEv_struct.Events.Events.onsettiming(specEv_struct.Events.Events.trialind==trial);
             offset=specEv_struct.Events.Events.offsettiming(specEv_struct.Events.Events.trialind==trial);
             onset_realtime=[];
             offset_realtime=[];
             for i=1:length(max_t)
                 onset_realtime(i)=EEG.times(find(round(tVec_assumed,3)==round(onset(i),3)));
                 offset_realtime(i)=EEG.times(find(round(tVec_assumed,3)==round(offset(i),3)));
             end
             event_onset=[event_onset, max(onset_realtime)];
             event_offset=[event_offset,max(offset_realtime)];
          
        yyaxis right
        plot(EEG.times,X{1}(:,trial),'w','Linewidth',2)
           
        title(strcat('Trial ',num2str(trial), '- ',num2str(length(max_t)), 'Events'))
            
        if sub_count==4 | trial==specEv_struct.Events.Events.trialind(end)
            print('-dpng','-r150',strcat('temp','.png'));
            blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
            Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
            Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500

            sub_count=0;
        end
           
end
close all    


%% now let's find the peak in the timeseries based on events

count=0;
for trial=event_trial
    if count==0
        figure('units','normalized','outerposition', [0 0 1 1]);
    end
    count=count+1;
    subplot(5,5,count)
    % plot timeseries
    plot(EEG.times,X{1}(:,trial),'Color',[.5 .5 .5])    
    hold on
    % plot where event is
    plot(event_onset(event_trial==trial):event_offset(event_trial==trial),X{1}(find(EEG.times==event_onset(event_trial==trial)):find(EEG.times==event_offset(event_trial==trial)),trial))
    % plot where max was detected
    plot(event_max(event_trial==trial), X{1}(find(EEG.times==event_max(event_trial==trial)),trial),'ro')
    % plot the actual trough here
    [trough,trough_i]=min(X{1}(find(EEG.times==event_onset(event_trial==trial)):find(EEG.times==event_offset(event_trial==trial)),trial));
    temp=1:length(EEG.times);
    temp=temp(find(EEG.times==event_onset(event_trial==trial)):find(EEG.times==event_offset(event_trial==trial)));
    trough_i=temp(trough_i);
    plot(EEG.times(trough_i),X{1}(trough_i,trial),'r*')
    
    title(strcat('Trial',num2str(trial)))
    
    
    if count==25 | trial==event_trial(end)
        count=0;
        print('-dpng','-r150',strcat('temp','.png'));
        blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
        Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
        Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500
    end
end

%% do it max_peak_locked?
%don't like these peaks so I'll exclude them from the trough-locked signal
%(for the actual TEP analysis, this won't matter cause I'll just go by
%event times, not trough times)
if partic==2
    bad=[12 36 37 60 63 67];
elseif partic==4
    bad=[14 38 43 62 73 92 93 95 101 105 133 137 138];
end
trough_lock=[];
t_interval=500;
trough_t=-t_interval:t_interval;
count=0;
figure
hold on
sanity_check=0;%if one, i'll just pick a random timeperiod in the trial, take the minimum and lock to that, to see if there's anything special about locking to beta
for trial=event_trial
    if all(trial~=bad)
        title('Trough Lock')
        count=count+1;
        [trough,trough_i]=min(X{1}(find(EEG.times==event_onset(event_trial==trial)):find(EEG.times==event_offset(event_trial==trial)),trial));
                if sanity_check==1
                    title('Sanity Check')
                    ran=randi(length(EEG.times)-(event_offset(event_trial==trial)-event_onset(event_trial==trial))-1);
                    [trough,trough_i]=min(X{1}(ran:ran+(event_offset(event_trial==trial)-event_onset(event_trial==trial)),trial));
                end
        
        temp=1:length(EEG.times);
        temp=temp(find(EEG.times==event_onset(event_trial==trial)):find(EEG.times==event_offset(event_trial==trial)));
                if sanity_check==1
                    temp=1:length(EEG.times);
                    temp=temp((ran:ran+(event_offset(event_trial==trial)-event_onset(event_trial==trial))));
                end
        trough_i=temp(trough_i);
        
        try
            trough_lock(count,:)=X{1}(trough_i-t_interval: trough_i+t_interval,trial);
        catch
            %might need to pad
            temp_data=[nan(1,t_interval),X{1}(:,trial)',nan(1,t_interval)];
            trough_lock(count,:)=temp_data(trough_i-t_interval+t_interval:trough_i+t_interval+t_interval);
        end
        plot(trough_t, trough_lock(count,:),'Color',[.5 .5 .5])
    end
end
plot(trough_t, nanmean(trough_lock),'Color','k','Linewidth',2)
print('-dpng','-r150',strcat('temp','.png'));
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500
         
    

%% now the actual analysis
EEG=data.EEG;
%put distance between TEP and beta in bins

Bins=[0 50 100 200 400 600 800 1000];
Bins=[0 100 500 1000];
Bins=[0 1000];
TEP_in_bins=[];

for trial=1:size(EEG.data,3)
    if any(trial==event_trial)
        for bin=1:length(Bins)-1
            if abs(event_offset(event_trial==trial))> Bins(bin) & abs(event_offset(event_trial==trial))<= Bins(bin+1)
                TEP_in_bins.(strcat('Bin',num2str(Bins(bin)),'to',num2str(Bins(bin+1))))(trial,:)=EEG.data(electr_oi_i,:,trial);
            end
        end
    else
        TEP_in_bins.None(trial,:)=EEG.data(electr_oi_i,:,trial);
    end       
end     
        
%delete empties
for bin=1:length(Bins)-1
    TEP_in_bins.(strcat('Bin',num2str(Bins(bin)),'to',num2str(Bins(bin+1))))(sum(TEP_in_bins.(strcat('Bin',num2str(Bins(bin)),'to',num2str(Bins(bin+1)))),2)==0,:)=[];
end
TEP_in_bins.None(sum(TEP_in_bins.None,2)==0,:)=[];
 



figure
clf
hold on
plot_times=[-100 400];
c=summer(length(Bins)-1);
if size(c,1)==2
    c(1,:)=[.875 .875 .219];
    c(2,:)=[.25 .5 0];
end
legtext={};
count=0;
for bin=1:length(Bins)-1
    count=count+1;
    p_t=EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2)));
    p_d=mean(TEP_in_bins.(strcat('Bin',num2str(Bins(bin)),'to',num2str(Bins(bin+1)))));
    p_d=p_d([find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))]);
    plot(p_t,p_d,'Color',c(bin,:),'Linewidth',2)
    legtext{count}=strcat('Bin',num2str(Bins(bin)),'to',num2str(Bins(bin+1)),' -- (',num2str(size(TEP_in_bins.(strcat('Bin',num2str(Bins(bin)),'to',num2str(Bins(bin+1)))),1)),')');
end
p_d=mean(TEP_in_bins.None);
p_d=p_d([find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))]);
plot(p_t,p_d,'Color','k','Linewidth',2)
legtext{count+1}=strcat('None -- (',num2str(size(TEP_in_bins.None,1)),')');
legend(legtext,'Location','southeast')


print('-dpng','-r150',strcat('temp','.png'));
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500
 
        
        
        
        
