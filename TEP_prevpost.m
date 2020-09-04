
%% Let's see if TEPs change over the course of the session

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

use_trialnumber=1; %decide if we use the number of the trial or the timing of the trial (based on urevents) to establish where in the session it ccured
partic=4;
data_path='F:\Brown\TMS\Pilot\';
Partic=[2,4];
dt=1;
TESAICA=0; %if 0, runICA for blink, if 1 TESAICA automatic (not done here but pick which preprocessed file you want)
electr_oi='C3';

%% load data
if TESAICA
    name='TESA';
else
    name='run';
end
data=load(strcat(data_path,strcat('Beta0',num2str(partic),'_TEP_1k_',name,'ICA_filt100')));
EEG=data.EEG;

%% load urevents
raw_path='F:\Brown\TMS\Pilot\';
if partic==2
    filename_ext='1_4';
else
    filename_ext='4_10';
end
urEEG = pop_loadset('filename',strcat('BETA0',num2str(partic),'_merge_',filename_ext,'.set'),'filepath',raw_path);
urevents=urEEG.urevent;
clear urEEG

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
% print('-dpng','-r150',strcat('temp','.png'));
%  blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
%  Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
% Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,850,550);%10,20,700,500
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERP Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-- Basic--
%first, let's just plot all trials at electr_oi to see variability
plot_times=[-100 400];
figure
hold on
for trial =1:size(EEG.data,3)
    plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial))
end
title(electr_oi)
ylabel('Amplitude')
xlabel('Time')
print('-dpng','-r150',strcat('temp','.png'));
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500

%-- Avg--
% Now let's avgerage over binned trials as we go
%rebaseline?baseline at all? made almost no diff
base=[-200 0];
for trial=1:size(EEG.data,3)
    for elec=1:length(EEG.chanlocs)
        EEG.data(elec,:,trial)=EEG.data(elec,:,trial)-mean(EEG.data(elec,[find(EEG.times==base(1)):find(EEG.times==base(2))],trial));
    end
end


%avg over a few trials at a time
plot_times=[-100 400];
if partic==2
    avg_over=8;
elseif partic==4
    avg_over=30;%14;
end
figure
clf
hold on
c=parula;%hot cold
c_count=0;
d_count=0;
l_count=0;
this_data=[];
Line=[];
increment_colour_by=round(length(c)/((size(EEG.data,3)/avg_over)+1));
for trial =1:size(EEG.data,3)
    if d_count<=avg_over
        d_count=d_count+1;
        this_data(d_count,:)= EEG.data(electr_oi_i,[find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))],trial);       
    end
    
    if size(this_data,1)==avg_over | trial==size(EEG.data,3)        
        c_count=c_count+increment_colour_by;
        l_count=l_count+1;
        Line(l_count,:)=plot(EEG.times(find(EEG.times==plot_times(1)):find(EEG.times==plot_times(2))),mean(this_data),'Color',c(c_count,:),'Linewidth',2);
        this_data=[];
        d_count=0;
    end

end
legend([Line(1),Line(end)],{strcat('First ',num2str(avg_over), 'Trials'),strcat('Last ',num2str(avg_over),'(',num2str(size(EEG.data,3)-avg_over*(l_count-1)), ')Trials')})
title(electr_oi)
ylabel('Amplitude')
xlabel('Time')
print('-dpng','-r150',strcat('temp','.png'));
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's characterise our peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%let's take smaller ERPs for this
ERP = pop_epoch( EEG, { 'S  1' }, [-.1 .3], 'epochinfo', 'yes');

%% find times associated with each trial (=peak)
count=0;
timing=[];
for event=1:size(EEG.event,2)
    if ERP.event(event).type=='S  1'
        count=count+1;
        timing(count)=ERP.event(event).latency;        
    end
end
if length(timing)~=size(EEG.data,3)
    fprintf('\nNr of trials doesn''t match latencies found. \n')
end


%% First: get a guesstimate of where the main peaks are
smooth = pop_epoch( EEG, { 'S  1' }, [0 .25], 'epochinfo', 'yes');
smooth=tesa_filtbutter(smooth,5,30,4,'bandpass');
N1_time=[]; N1_amp=[];N1_interval=[];
P1_time=[]; P1_amp=[];P1_interval=[];
N2_time=[]; N2_amp=[];N2_interval=[];
P2_time=[]; P2_amp=[];P2_interval=[];
N3_time=[]; N3_amp=[];N3_interval=[];
P3_time=[]; P3_amp=[];P3_interval=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotcheck=0;
plotcheck_fig=figure;
if ppt
    %select a few random trials to plot
    ppt_fig=figure('units','normalized','outerposition',[0 0 1 1]);
    ppt_trials=randperm(size(EEG.data,3),9);
    ppt_count=0;
end
%first, find N2 (should be biggest)
for trial =1:size(ERP.data,3)
    figure(plotcheck_fig)
    clf
    if ppt & any(trial==ppt_trials)
        ppt_count=ppt_count+1;
        figure(ppt_fig)
        subplot(3,3,ppt_count)
        hold on
    end
    
    plot(smooth.times,smooth.data(electr_oi_i,:,trial),'Linewidth',1,'Color',[.6 .6 .6]);
    hold on
    plot(ERP.times,ERP.data(electr_oi_i,:,trial),'Linewidth',2,'Color','k');    
    
    %% N2
    look_here=[60 170];
    N2_time(trial)=smooth.times(find(smooth.data(electr_oi_i,:,trial)==min(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial))));    
    N2_amp(trial)=min(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial));
    N2_interval(trial,:)=[N2_time(trial)-25, N2_time(trial)+25];
    
    plot(N2_time(trial),N2_amp(trial),'o','Color',[1 0 .5])%,'MarkerFaceColor',[1 0 .5])
    
    look_here=N2_interval(trial,:);
    if partic==4 & trial==67
        look_here(1)=look_here(1)-20;
    end
    N2_time(trial)=ERP.times(find(ERP.data(electr_oi_i,:,trial)==min(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial))));    
    N2_amp(trial)=min(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial));
    plot(N2_time(trial),N2_amp(trial),'o','Color',[1 0 .5],'MarkerFaceColor',[1 0 .5])
    tx=text(N2_time(trial),N2_amp(trial)-1,'N2');
    tx.FontWeight='bold';
    tx.Color=[1 0 .5];
    
    %% P1
    if partic==2
        look_here=[max(50,N2_time(trial)-80) N2_time(trial)];
    elseif partic==4
        look_here=[max(20,N2_time(trial)-80) N2_time(trial)];
    end
    
    P1_time(trial)=smooth.times(find(smooth.data(electr_oi_i,:,trial)==max(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial))));    
    P1_amp(trial)=max(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial));
    P1_interval(trial,:)=[max(0,P1_time(trial)-25), P1_time(trial)+25];
    if partic==4 & (trial==11 |trial==37 | trial==43 | trial==68 | trial==90)
        P1_interval(trial,:)=[P1_time(trial)-15, P1_time(trial)+25];
    end
        
  
    plot(P1_time(trial),P1_amp(trial),'o','Color',[.344 .906 .344])%,'MarkerFaceColor',[.344 .906 .344])

    look_here=P1_interval(trial,:);
    if partic==4 & trial==90
        look_here(1)=28;
    end
    P1_time(trial)=ERP.times(find(ERP.data(electr_oi_i,:,trial)==max(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial))));    
    P1_amp(trial)=max(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial));
    plot(P1_time(trial),P1_amp(trial),'o','Color',[.344 .906 .344],'MarkerFaceColor',[.344 .906 .344])
    tx=text(P1_time(trial),P1_amp(trial)+1,'P1');
    tx.FontWeight='bold';
    tx.Color=[.344 .906 .344];
    
    %% N1
    look_here=[max(0,P1_time(trial)-40) P1_time(trial)];
    
    N1_time(trial)=smooth.times(find(smooth.data(electr_oi_i,:,trial)==min(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial))));    
    N1_amp(trial)=min(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial));
    N1_interval(trial,:)=[max(0,N1_time(trial)-15), N1_time(trial)+15];

    plot(N1_time(trial),N1_amp(trial),'o','Color',[1 .5 .75])%,'MarkerFaceColor',[1 .5 .75])

    look_here=N1_interval(trial,:);
    if partic==4 & trial==68
        look_here(2)=18;
    end
    N1_time(trial)=ERP.times(find(ERP.data(electr_oi_i,:,trial)==min(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial))));    
    N1_amp(trial)=min(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial));
    plot(N1_time(trial),N1_amp(trial),'o','Color',[1 .5 .75],'MarkerFaceColor',[1 .5 .75])
    tx=text(N1_time(trial),N1_amp(trial)-1,'N1');
    tx.FontWeight='bold';
    tx.Color=[1 .5 .75];
    
    %% P2
    look_here=[N2_time(trial) max(180,N2_time(trial)+50)];
    
    P2_time(trial)=smooth.times(find(smooth.data(electr_oi_i,:,trial)==max(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial))));    
    P2_amp(trial)=max(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial));
    P2_interval(trial,:)=[P2_time(trial)-20, P2_time(trial)+20];
    
    plot(P2_time(trial),P2_amp(trial),'o','Color',[.25 .75 .25])%,'MarkerFaceColor',[.25 .75 .25])

    look_here=P2_interval(trial,:);
    P2_time(trial)=ERP.times(find(ERP.data(electr_oi_i,:,trial)==max(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial))));    
    P2_amp(trial)=max(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial));
    plot(P2_time(trial),P2_amp(trial),'o','Color',[.25 .75 .25],'MarkerFaceColor',[.25 .75 .25])
    tx=text(P2_time(trial),P2_amp(trial)+1,'P2');
    tx.FontWeight='bold';
    tx.Color=[.25 .75 .25];
    
    %% N3
    look_here=[P2_time(trial) max(200,P2_time(trial)+50)];
    
    N3_time(trial)=smooth.times(find(smooth.data(electr_oi_i,:,trial)==min(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial))));    
    N3_amp(trial)=min(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial));
    N3_interval(trial,:)=[N3_time(trial)-15, N3_time(trial)+15];

    plot(N3_time(trial),N3_amp(trial),'o','Color',[.75 0 .375])%,'MarkerFaceColor',[.75 0 .375])

    look_here=N3_interval(trial,:);
    N3_time(trial)=ERP.times(find(ERP.data(electr_oi_i,:,trial)==min(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial))));    
    N3_amp(trial)=min(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial));
    plot(N3_time(trial),N3_amp(trial),'o','Color',[.75 0 .375],'MarkerFaceColor',[.75 0 .375])
    tx=text(N3_time(trial),N3_amp(trial)-1,'N3');
    tx.FontWeight='bold';
    tx.Color=[.75 0 .375];
    
    %% P3
    look_here=[N3_time(trial) min(smooth.times(end),N3_time(trial)+50)];
    
    P3_time(trial)=smooth.times(find(smooth.data(electr_oi_i,:,trial)==max(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial))));    
    P3_amp(trial)=max(smooth.data(electr_oi_i,[find(smooth.times==look_here(1)):find(smooth.times==look_here(2))],trial));
    P3_interval(trial,:)=[P3_time(trial)-20, P3_time(trial)+20];

    plot(P3_time(trial),P3_amp(trial),'o','Color',[0 .5 0])%,'MarkerFaceColor',[0 .5 0])

    look_here=P3_interval(trial,:);
    P3_time(trial)=ERP.times(find(ERP.data(electr_oi_i,:,trial)==max(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial))));    
    P3_amp(trial)=max(ERP.data(electr_oi_i,[find(ERP.times==look_here(1)):find(ERP.times==look_here(2))],trial));
    plot(P3_time(trial),P3_amp(trial),'o','Color',[0 .5 0],'MarkerFaceColor',[0 .5 0])
    tx=text(P3_time(trial),P3_amp(trial)+1,'P3');
    tx.FontWeight='bold';
    tx.Color=[0 .5 0];
    
    title(strcat('Trial',num2str(trial)))
    if plotcheck==1
        temp=input('h');
    end
end
if ppt
    close(plotcheck_fig)
    figure(ppt_fig)
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',0,0,950,540);%10,20,700,500
end

    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's characterise our peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  - Amplitude (max trough)
%  - Mean Amplitude (mean amp in 50ms window around max)
%  - Latency (of max)

% For N2 onlly, we''l also do:
%  - FWHM (full width at half maximum - time width at halfway point between max and 0)
%  - Slope pre peak
%  - Slope post peak
Peaks={'N1','P1','N2','P2','N3','P3'};
Mean_amp_time=[20,20,50,30,30,30];
Slope_time=[20];%window from peak, in both directions (N2 only)
for pick_peak_oi=1:length(Peaks)
    Peak=Peaks{pick_peak_oi};

    plotcheck=0; %plot each trial to check
    if plotcheck==1
        plotcheck_fig=figure;
    end

    final_fig=figure;
    title(strcat(Peak,'- all'))
    hold on
    xlabel('Trial')
    ylabel('Amplitude')

    mean_amp_time=Mean_amp_time(pick_peak_oi);
    reclass=[];
    reject=[];
    peak_time=[];
    peak_amp=[];
    mean_amp=[];
    fwhm=[];
    for trial =1:size(ERP.data,3)
        peak_time_temp=eval(strcat(Peak, '_time(trial)'));
        peak_amp_temp=eval(strcat(Peak, '_amp(trial)'));
        peak_i=find(ERP.times==peak_time_temp);

        if plotcheck==1
            figure(plotcheck_fig)
            clf
            hold on
            title(strcat('Trial',num2str(trial)))
            plot(ERP.times,ERP.data(electr_oi_i,:,trial),'Color',[.5 .5 .5])
            plot(peak_time_temp,peak_amp_temp,'ro')
            classify=input('Accept(Enter)/Correct(1)/Reject(0):');
            if classify== 1
                reclass=[reclass,trial];
                peak_amp_temp=input('Set New Trough Amp:');
                peak_time_temp=input('Set New Trough Time:');
                peak_i=find(ERP.times==peak_time_temp); 
                plot(peak_time_temp,peak_amp_temp,'ro')
            elseif classify==0
                reject=[reject,trial];
            end
        else
            if partic==2 & trial==69 & pick_peak_oi==3 %already did this manually
                peak_amp_temp=-6.3618;
                peak_time_temp=122;
                peak_i=find(ERP.times==peak_time_temp);           
            end
        end
        %mean amp
        mean_amp(trial)=mean(ERP.data(electr_oi_i,[peak_i-mean_amp_time/2:peak_i+mean_amp_time/2],trial));
        %fwhm
        if Peak=='N2'
            halfmax=peak_amp_temp/2;
            try 
                limit=.5;
                [sorted,sort_i]=sort(abs(ERP.data(electr_oi_i,:,trial)-halfmax));
                sort_i=sort_i(sorted<limit);
                sort_i_post=sort_i(sort_i-peak_i>0);
                sort_i_pre=sort_i(sort_i-peak_i<0);
                postpeak_i=min(sort_i_post);
                prepeak_i=max(sort_i_pre);
                if ERP.times(postpeak_i)-peak_time_temp >45 | peak_time_temp-ERP.times(prepeak_i)>45 | (partic==4 & trial==2)
                    breaktry
                end
                fwhm(trial)=ERP.times(postpeak_i)-ERP.times(prepeak_i);
            catch
                limit=1;
                [sorted,sort_i]=sort(abs(ERP.data(electr_oi_i,:,trial)-halfmax));
                sort_i=sort_i(sorted<limit);
                sort_i_post=sort_i(sort_i-peak_i>0);
                sort_i_pre=sort_i(sort_i-peak_i<0);
                postpeak_i=min(sort_i_post);
                prepeak_i=max(sort_i_pre);
                fwhm(trial)=ERP.times(postpeak_i)-ERP.times(prepeak_i);
            end
            
%         % also look at slope if N2
%         %slope based on fixed time window from peak
%         %pre peak slope
%         slope_data=ERP.data(electr_oi_i,[1+peak_i-Slope_time*dt:peak_i],trial);        
%         P_pre=polyfit(ERP.times(1+peak_i-Slope_time*dt:peak_i),slope_data,1);
%         pre_slope(trial)=P_pre(1);
%         %post peak slope
%         slope_data=ERP.data(electr_oi_i,[peak_i:peak_i+Slope_time*dt-1],trial);        
%         P_post=polyfit(ERP.times([peak_i:peak_i+Slope_time*dt-1]),slope_data,1);
%         post_slope(trial)=P_post(1);
        
        %slope based on time window from peak to fwhm
        % keep in miknd, this woln't giv e you a line between the two
        % points. the line will be whatever reduces teh distance to all
        % points bedtween tholse two points
        %pre peak slope
        slope_data=ERP.data(electr_oi_i,[prepeak_i:peak_i],trial);        
        P_pre=polyfit(ERP.times([prepeak_i:peak_i]),slope_data,1);
        pre_slope(trial)=P_pre(1);
        %post peak slope
        slope_data=ERP.data(electr_oi_i,[peak_i:postpeak_i],trial);        
        P_post=polyfit(ERP.times([peak_i:postpeak_i]),slope_data,1);
        post_slope(trial)=P_post(1);
        
        end


        peak_time(trial)=peak_time_temp;
        peak_amp(trial)=peak_amp_temp;
        if plotcheck==1
            figure(plotcheck_fig)
            %fwhm
            plot(ERP.times(postpeak_i),ERP.data(electr_oi_i,postpeak_i,trial),'bo')
            plot(ERP.times(prepeak_i),ERP.data(electr_oi_i,prepeak_i,trial),'bo')
            %slopes
%             %slope based on fixed time window from peak
%             xfit=linspace(peak_time(trial)-Slope_time+1,peak_time(trial),100);
%             yfit=polyval(P_pre,xfit);
%             plot(xfit,yfit,'-','Color','g','Linewidth',1)
%             
%             xfit=linspace(peak_time(trial),Slope_time+peak_time(trial)-1,100);
%             yfit=polyval(P_post,xfit);
%             plot(xfit,yfit,'-','Color','g','Linewidth',1)

            %slope based on time window from peak to fwhm
            xfit=linspace(ERP.times(prepeak_i),peak_time(trial),2);
            yfit=polyval(P_pre,xfit);
            plot(xfit,yfit,'-','Color','g','Linewidth',1)
            
            xfit=linspace(peak_time(trial),ERP.times(postpeak_i),100);
            yfit=polyval(P_post,xfit);
            plot(xfit,yfit,'-','Color','g','Linewidth',1)
            
            
            temp=input('FWHM(Enter)');
        end
        figure(final_fig)
        plot(ERP.times,ERP.data(electr_oi_i,:,trial),'Color',[.7 .7 .7],'Linewidth',.5)
        plot(peak_time_temp,peak_amp_temp,'o','Color',[0 0 .625])
    end
    if plotcheck
        close(plotcheck_fig)
    end
    figure(final_fig)
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500



    %% let's plot the relationships
    figure('units','normalized','outerposition', [1.5 -0.07 0.7 0.94]);
    if Peak=='N2'
        vals={'amp', 'mean','time','fwhm','preslope','postslope'};
    else
        vals={'amp', 'mean','time'};
    end
    c=parula(size(ERP.data,3));

    for val=1:length(vals)
        subplot(2,3,val)
        if val==1
            value=peak_amp;
            ylab='Peak Amplitude';
        elseif val==2
            value=mean_amp;
            ylab=('Mean Amplitude');
        elseif val==3
            value=peak_time;
            ylab='Latency';
        elseif val==4
            value=fwhm;
            ylab='FWHM';
        elseif val==5
            value=pre_slope;
            ylab='Slope';
        elseif val==6
            value=post_slope;
            ylab='Slope';
        end
        
        if use_trialnumber
            %line of best fit
            coeffs=polyfit([1:size(ERP.data,3)],value,1);
            xfit=linspace(1,size(ERP.data,3),1000);
            yfit=polyval(coeffs,xfit);
            plot(xfit,yfit,'-','Color',[.8 .8 .8],'Linewidth',3)
            %also corr but this obvs doesn't mean anything
            [r,p]=corrcoef([1:size(ERP.data,3)],value);
        else 
            coeffs=polyfit([timing],value,1);
            xfit=timing;
            yfit=polyval(coeffs,xfit);
            plot(xfit,yfit,'-','Color',[.8 .8 .8],'Linewidth',3)
            %also corr but this obvs doesn't mean anything
            [r,p]=corrcoef(timing,value);
        end
        hold on
        for trial=1:size(ERP.data,3)
            if use_trialnumber
               plot(trial,value(trial),'o','Color',c(trial,:),'MarkerFaceColor',c(trial,:))
            else
                plot(timing(trial),value(trial),'o','Color',c(trial,:),'MarkerFaceColor',c(trial,:))
            end
        end
        ylims=ylim;
        textpos=+(ylims(2)-ylims(1))/20;
        if r(1,2)<0
            textpos=textpos*(-1);
        end
        sig='';
        if p(1,2)<.05
            sig='*';
        end
        if use_trialnumber
           tx=text(xfit(end)-round(size(EEG.data,3)/5),double(yfit(end)+textpos),sprintf('r = %.2f%s',r(1,2),sig))
           xlim([0 size(EEG.data,3)+1])
        else
            tx=text(xfit(end)-round(xfit(end)/5),double(yfit(end)+textpos),sprintf('r = %.2f%s',r(1,2),sig))      
            xlim([timing(1) timing(end)])
        end
        tx.Color=[.8 .8 .8];
        tx.FontWeight='bold';
        ylabel(strcat(ylab))
        xlabel('Trial')
        title(strcat(Peak,'-',ylab))
        
    end
    print('-dpng','-r150',strcat('temp','.png'));
    blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
    Slide1 = Presentation.Slides.AddSlide(1,blankSlide);
    Image1 = Slide1.Shapes.AddPicture(strcat(cd,'/temp','.png'),'msoFalse','msoTrue',120,0,700,540);%10,20,700,500

end



%Presentation.SaveAs(strcat(cd,'\2Beta0',num2str(partic),'.ppt'))

close all



















