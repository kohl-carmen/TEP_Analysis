
%% Plot single pulse TEPs
% takes data from TEP_preproc_single

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%pick what you want
S1_80=1; %if not:S1_100
ICA=1;
    TESAICA=0; %if not: runica
refilt=1; %if not: no second filter applied
partic=2;
ppt=1;

data_path='C:\Users\ckohl\Desktop\Current\TMS\CurrentData\EEG\';
eeglab_dir='C:\Users\ckohl\Documents\MATLAB\eeglab2020_0';
cd(eeglab_dir)
eeglab

if ppt
    h = actxserver('PowerPoint.Application');
    Presentation = h.Presentation.Add;
end

Partic=[2,4];
if ~ICA
    name='noica';
else
    name='runICA';
    if TESAICA
        name='TESAICA';
    end
end
if refilt
    name=strcat(name,'_refilt');
else
    name=strcat(name,'_nofilt');
end
nameplus='';
if ~S1_80
    nameplus='100_';
end
%load
load(strcat(data_path,'Beta0',num2str(partic),'_TEP_',nameplus,name));

event='S 13';
electr_oi='C3';

%find electrode
for chan= 1:length(EEG.chanlocs)
    if length(EEG.chanlocs(chan).labels)==2
        if EEG.chanlocs(chan).labels==electr_oi
            electr_oi_i=chan;
        end
    end
end

%reepoch
EEG = pop_epoch( EEG, { event }, [-.2 1], 'epochinfo', 'yes');
EEG = pop_rmbase( EEG, [-200   0]);
%get rid of EOG
EEG=pop_select(EEG, 'nochannel',[64 65]);

figure; pop_timtopo(EEG, [-50  500], [10 20 100 200]);

if ppt
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(1);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,850,550);%10,20,700,500
end
%% ERP IMAGE
figure; 
pop_erpimage(EEG,1, [electr_oi_i],[[]],strcat(EEG.chanlocs(electr_oi_i).labels, '- Trials: ',num2str(size(EEG.data,3))),10,1,{},[],'' ,'yerplabel','\muV','erp','on','limits',[-50 299 NaN NaN NaN NaN NaN NaN] ,'cbar','on','topo', { [electr_oi_i] EEG.chanlocs EEG.chaninfo } );
if ppt
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(1);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,850,550);%10,20,700,500
end

%% All Electrodes
subplarray=[nan,nan,nan,nan,1,60,2,nan,nan,nan,nan;nan,nan,nan,52,38,30,39,53,nan,nan,nan;nan,11,46,3,32,17,33,4,47,12,nan;nan,54,25,40,21,31,22,41,26,55,nan;nan,13,48,5,34,18,35,6,49,14,nan;29,56,27,42,23,61,24,43,28,57,nan;nan,15,50,7,36,19,37,8,51,16,nan;nan,nan,nan,58,44,62,45,59,nan,nan,nan;nan,nan,nan,nan,9,63,10,nan,nan,nan,nan;nan,nan,nan,nan,nan,20,nan,nan,nan,nan,nan];                
Data=EEG;
figure('units','normalized','outerposition',[0 0 1 1])
clf
hold on
% time=[-1*dt+1:1.*dt];
time=EEG.times;
y_lims=[0 0];
colour=[.4 .4 .4];
for c=1:length(EEG.chanlocs)
    sub_i=find(subplarray'==c);
    sub=subplot(10,11,sub_i);     
    hold on
%   title(EEG.chanlocs(c).labels)
    plot(zeros(size([-10*10^6:1*10^6:10*10^6])),[-10*10^6:1*10^6:10*10^6],'Color',[.5 .5 .5])
    plot(time,zeros(size(time)),'Color',[.5 .5 .5])
    plot(time,mean(Data.data(c,:,:),3));
        %make standard error
        SE_upper=[];
        SE_lower=[];
        for i=1:size(Data.data,2)
            se=std(Data.data(c,i,:))./sqrt(length(Data.data(c,i,:)));
            SE_upper(i)=mean(Data.data(c,i,:))+se;
            SE_lower(i)=mean(Data.data(c,i,:))-se;
        end 
        tempx=[time,fliplr(time)];
        tempy=[SE_upper,fliplr(SE_lower)];
        A=fill(tempx,tempy,'k');
        A.EdgeColor='none';
        A.FaceColor=colour;
        A.FaceAlpha=.2;
    xlim([-50 300])
    y_lims(1)=min([y_lims(1),min(mean(Data.data(c,:,:),3))]);
    y_lims(2)=max([y_lims(2),max(mean(Data.data(c,:,:),3))]);
    set(gca,'visible','off')
    b=annotation('textbox','String',Data.chanlocs(c).labels);
    b.Position=sub.Position;
    b.Position(1)=b.Position(1)+.03;
%   b.FontSize=14;
    b.EdgeColor='none';
end

for c=1:length(EEG.chanlocs)
    sub_i=find(subplarray'==c);
    subplot(10,11,sub_i) 
    ylim(y_lims)
    % colour coding
    thirds=linspace(0,y_lims(2)-y_lims(1),4);
    this_diff=max(mean(Data.data(c,:,:),3))-min(mean(Data.data(c,:,:),3));
%     if this_diff<thirds(2)
%         plot(time,mean(Data.data(c,:,:),3),'Linewidth',2,'Color',[.344 .906 .344]);
%     elseif this_diff<thirds(3)
%         plot(time,mean(Data.data(c,:,:),3),'Linewidth',2,'Color',[1 .625 .25]);
%     else
%          plot(time,mean(Data.data(c,:,:),3),'Linewidth',2,'Color',[.75 0 0]);
%     end
    plot(time,mean(Data.data(c,:,:),3),'Linewidth',1,'Color',[.3 .3 .3]);
end        
if ppt
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(1);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,850,550);%10,20,700,500
end


%% Topos
% fig1.Renderer='Painters';
dataoi=EEG;
topo_time=[0:10 :200]
topo_time_interval=topo_time(2)-topo_time(1);
topo_samples=topo_time+200;%topo_time*dt + 1*dt;
max_clim=[0 0];

for i=[1:length(topo_samples)]
    temp=mean(mean(dataoi.data(:,topo_samples(i)-topo_time_interval:topo_samples(i),:),3),2);
    max_clim(1)=min(max_clim(1),min(temp));
    max_clim(2)=max(max_clim(2),max(temp));
end
figure
hold on
for i=1:length(topo_time)-1
    subplot(4,5,i)
    topoplot(mean(mean(dataoi.data(:,topo_samples(i)-topo_time_interval:topo_samples(i),:),3),2),EEG.chanlocs,'maplimits',max_clim./3);
    title(strcat(num2str(topo_time(i)-topo_time_interval),'-',num2str(topo_time(i)),'ms'));
end

if ppt
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(1);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,850,550);%10,20,700,500
end


%% ERPs with topos
dataoi=EEG;
figure
clf
time=EEG.times;
subplot(2,5,[6:10]);
hold on
% plot(eeg_time,mean(dataoi.data(:,eeg_samples,:),3));
plot(time,mean(dataoi.data(electr_oi_i,:,:),3),'k','Linewidth',2);
xlim([-50 250])
%  make standard error
SE_upper=[];
SE_lower=[];
for i=1:size(dataoi.data,2)
    se=std(dataoi.data(electr_oi_i,i,:))./sqrt(length(dataoi.data(electr_oi_i,i,:)));
    SE_upper(i)=mean(dataoi.data(electr_oi_i,i,:))+se;
    SE_lower(i)=mean(dataoi.data(electr_oi_i,i,:))-se;
end 
tempx=[time,fliplr(time)];
tempy=[SE_upper,fliplr(SE_lower)];
A=fill(tempx,tempy,'k');
A.EdgeColor='none';
A.FaceColor=colour;
A.FaceAlpha=.2;
% for i=1:size(dataoi.data,3)/10
%     plot(eeg_time,(dataoi.data(electr_oi_i,eeg_samples,i)),'Color',[.5 .5 .5]);
% end
y=ylim;
x=xlim;
ylim(y)
xlim(x)
xlabel('Time');
ylabel('Amplitude');
set(gca,'Clipping','Off')
topo_time=[10 20 50 100 150];%[50 70 100 150 200];%[0:2];%ms
topo_time_interval=10;%topo_time(2)-topo_time(1);
topo_samples=topo_time+200;%topo_time*dt + 1*dt;
topo_loc=[-25 40 100 160 220];
max_clim=[0 0];

for i=[1:length(topo_samples)]
    temp=mean(mean(dataoi.data(:,topo_samples(i)-topo_time_interval:topo_samples(i)+topo_time_interval,:),3),2);
    max_clim(1)=min(max_clim(1),min(temp));
    max_clim(2)=max(max_clim(2),max(temp));
end
for i=1:length(topo_samples)
    subplot(2,5,i)
    topoplot(mean(mean(dataoi.data(:,topo_samples(i)-topo_time_interval:topo_samples(i)+topo_time_interval,:),3),2),EEG.chanlocs,'maplimits',max_clim./3);
%     cbar('horiz',0,round(max_clim/1000),3)
%     title(strcat(num2str(topo_time(i)-topo_time_interval),'-',num2str(topo_time(i)+topo_time_interval),'ms'));
    title(strcat(num2str(topo_time(i)),'ms'))
    subplot(2,5,[6:10]);
    h=line([topo_time(i)-topo_time_interval,topo_loc(i)],      [mean(EEG.data(electr_oi_i,topo_samples(i)-topo_time_interval,:),3),y(2)*2]);
    h.Color='k';
    h=line([topo_time(i)+topo_time_interval,topo_loc(i)],      [mean(EEG.data(electr_oi_i,topo_samples(i)+topo_time_interval,:),3),y(2)*2]);
    h.Color='k';
end

if ppt
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(1);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,850,550);%10,20,700,500
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare to SI_100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tep1=EEG;
%get SI_100
load(strcat(data_path,'Beta0',num2str(partic),'_TEP_100_',nameplus,name));
EEG = pop_epoch( EEG, { event }, [-.2 1], 'epochinfo', 'yes');
EEG = pop_rmbase( EEG, [-200   0]);
EEG=pop_select(EEG, 'nochannel',[64 65]);
tep2=EEG;

%% ERPs with topos
dataoi=tep1;
dataoi2=tep2;
figure
clf
time=EEG.times;
subplot(3,5,[6:10]);
hold on
% plot(eeg_time,mean(dataoi.data(:,eeg_samples,:),3));
colour2=[1 .625 .25];
plot(time,mean(dataoi2.data(electr_oi_i,:,:),3),'Color',colour2,'Linewidth',2);
plot(time,mean(dataoi.data(electr_oi_i,:,:),3),'k','Linewidth',2);
legend('SI-100','SI-80','Location','southeast')
xlim([-50 250])
%  make standard error
SE_upper=[];
SE_lower=[];
for i=1:size(dataoi.data,2)
    se=std(dataoi.data(electr_oi_i,i,:))./sqrt(length(dataoi.data(electr_oi_i,i,:)));
    SE_upper(i)=mean(dataoi.data(electr_oi_i,i,:))+se;
    SE_lower(i)=mean(dataoi.data(electr_oi_i,i,:))-se;
end 
tempx=[time,fliplr(time)];
tempy=[SE_upper,fliplr(SE_lower)];
A=fill(tempx,tempy,'k');
A.EdgeColor='none';
A.FaceColor=colour;
A.FaceAlpha=.2;
%  make standard error 2
SE_upper=[];
SE_lower=[];
for i=1:size(dataoi.data,2)
    se=std(dataoi2.data(electr_oi_i,i,:))./sqrt(length(dataoi2.data(electr_oi_i,i,:)));
    SE_upper(i)=mean(dataoi2.data(electr_oi_i,i,:))+se;
    SE_lower(i)=mean(dataoi2.data(electr_oi_i,i,:))-se;
end 
tempx=[time,fliplr(time)];
tempy=[SE_upper,fliplr(SE_lower)];
A=fill(tempx,tempy,[1 .625 .25]);
A.EdgeColor='none';
A.FaceColor=colour2;
A.FaceAlpha=.2;
legend('SI-100','SI-80','Location','southeast')
% for i=1:size(dataoi.data,3)/10
%     plot(eeg_time,(dataoi.data(electr_oi_i,eeg_samples,i)),'Color',[.5 .5 .5]);
% end
y=ylim;
x=xlim;
ylim(y)
xlim(x)
xlabel('Time');
ylabel('Amplitude');
topo_time=[10 20 50 100 150];%[50 70 100 150 200];%[0:2];%ms
topo_time_interval=10;%topo_time(2)-topo_time(1);
topo_samples=topo_time+200;%topo_time*dt + 1*dt;
topo_loc=[-25 40 100 160 220];
max_clim=[0 0];

for i=[1:length(topo_samples)]
    temp=mean(mean(dataoi.data(:,topo_samples(i)-topo_time_interval:topo_samples(i)+topo_time_interval,:),3),2);
    max_clim(1)=min(max_clim(1),min(temp));
    max_clim(2)=max(max_clim(2),max(temp));
end

for i=1:length(topo_samples)
    %top topos
    subplot(3,5,i)
    topoplot(mean(mean(dataoi.data(:,topo_samples(i)-topo_time_interval:topo_samples(i)+topo_time_interval,:),3),2),EEG.chanlocs,'maplimits',max_clim./3);
    title(strcat(num2str(topo_time(i)),'ms'))
    %bottom topos
    subplot(3,5,10+i)
    topoplot(mean(mean(dataoi2.data(:,topo_samples(i)-topo_time_interval:topo_samples(i)+topo_time_interval,:),3),2),EEG.chanlocs,'maplimits',max_clim./3);
    a=title(strcat(num2str(topo_time(i)),'ms'));
    a.Color=colour2;
end

if ppt
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(1);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,850,550);%10,20,700,500
end





%% All Electrodes
subplarray=[nan,nan,nan,nan,1,60,2,nan,nan,nan,nan;nan,nan,nan,52,38,30,39,53,nan,nan,nan;nan,11,46,3,32,17,33,4,47,12,nan;nan,54,25,40,21,31,22,41,26,55,nan;nan,13,48,5,34,18,35,6,49,14,nan;29,56,27,42,23,61,24,43,28,57,nan;nan,15,50,7,36,19,37,8,51,16,nan;nan,nan,nan,58,44,62,45,59,nan,nan,nan;nan,nan,nan,nan,9,63,10,nan,nan,nan,nan;nan,nan,nan,nan,nan,20,nan,nan,nan,nan,nan];                
Data=tep1;
Data2=tep2;
figure('units','normalized','outerposition',[0 0 1 1])
clf
hold on
% time=[-1*dt+1:1.*dt];
time=EEG.times;
y_lims=[0 0];
colour=[.4 .4 .4];
for c=1:length(EEG.chanlocs)
    sub_i=find(subplarray'==c);
    sub=subplot(10,11,sub_i);     
    hold on
%   title(EEG.chanlocs(c).labels)
    plot(zeros(size([-10*10^6:1*10^6:10*10^6])),[-10*10^6:1*10^6:10*10^6],'Color',[.5 .5 .5])
    plot(time,zeros(size(time)),'Color',[.5 .5 .5])
    plot(time,mean(Data2.data(c,:,:),3));
    plot(time,mean(Data.data(c,:,:),3),'Color',colour2);
        %make standard error
        SE_upper=[];
        SE_lower=[];
        for i=1:size(Data.data,2)
            se=std(Data.data(c,i,:))./sqrt(length(Data.data(c,i,:)));
            SE_upper(i)=mean(Data.data(c,i,:))+se;
            SE_lower(i)=mean(Data.data(c,i,:))-se;
        end 
        tempx=[time,fliplr(time)];
        tempy=[SE_upper,fliplr(SE_lower)];
        A=fill(tempx,tempy,'k');
        A.EdgeColor='none';
        A.FaceColor=colour;
        A.FaceAlpha=.2;
        %make standard error2
        SE_upper=[];
        SE_lower=[];
        for i=1:size(Data.data,2)
            se=std(Data2.data(c,i,:))./sqrt(length(Data2.data(c,i,:)));
            SE_upper(i)=mean(Data2.data(c,i,:))+se;
            SE_lower(i)=mean(Data2.data(c,i,:))-se;
        end 
        tempx=[time,fliplr(time)];
        tempy=[SE_upper,fliplr(SE_lower)];
        A=fill(tempx,tempy,colour2);
        A.EdgeColor='none';
        A.FaceColor=colour2;
        A.FaceAlpha=.2;
    xlim([-50 300])
    y_lims(1)=min([y_lims(1),min(mean(Data.data(c,:,:),3))]);
    y_lims(2)=max([y_lims(2),max(mean(Data.data(c,:,:),3))]);
    set(gca,'visible','off')
    b=annotation('textbox','String',Data.chanlocs(c).labels);
    b.Position=sub.Position;
    b.Position(1)=b.Position(1)+.03;
%   b.FontSize=14;
    b.EdgeColor='none';
end

for c=1:length(EEG.chanlocs)
    sub_i=find(subplarray'==c);
    subplot(10,11,sub_i) 
    ylim(y_lims)
    % colour coding
    thirds=linspace(0,y_lims(2)-y_lims(1),4);
    this_diff=max(mean(Data.data(c,:,:),3))-min(mean(Data.data(c,:,:),3));
%     if this_diff<thirds(2)
%         plot(time,mean(Data.data(c,:,:),3),'Linewidth',2,'Color',[.344 .906 .344]);
%     elseif this_diff<thirds(3)
%         plot(time,mean(Data.data(c,:,:),3),'Linewidth',2,'Color',[1 .625 .25]);
%     else
%          plot(time,mean(Data.data(c,:,:),3),'Linewidth',2,'Color',[.75 0 0]);
%     end
    plot(time,mean(Data.data(c,:,:),3),'Linewidth',1,'Color',[.3 .3 .3]);
end        
if ppt
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(1);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,850,550);%10,20,700,500
end

close all