
supra=0 %if 0, default
all_times=[0 1 2 3 4 5 6 7 10 12 15 20]% 40 60 80 100 120 140 160 180 200 220 240 280 300 340 400];
all_times=[0,1,5,10,40,100];
% all_times=[20  100  140  180  220  280]
initial_times=[26,63,137];
figure
leg_name={};
line=[];
hold on
colour=autumn(length(all_times));
for time_i=1:length(all_times)
    if supra
        file_dir=strcat('C:\Users\ckohl\hnn_out\data\TEPpred_supra_Sarah_2_',num2str(all_times(time_i)),'ms/dpl.txt');
        onsets=[-195,-210];
    else
        file_dir=strcat('C:\Users\ckohl\hnn_out\data\TEPpred_default_Sarah_2_',num2str(all_times(time_i)),'ms/dpl.txt');
        onsets=[-145,-150];
    end
    %load
    dpl = readtable(file_dir);
    dpl = table2array(dpl);
    line(time_i)=plot(dpl(:,1),dpl(:,2), 'Color',colour(time_i,:), 'Linewidth',2);
    leg_name{end+1}=strcat(num2str(all_times(time_i)),'ms');
    plot([all_times(time_i),all_times(time_i)],onsets, 'Color',colour(time_i,:), 'Linewidth',2)
end
if supra
    file_dir=strcat('C:\Users\ckohl\hnn_out\data\supra_Sarah\dpl.txt');
else    
    file_dir=strcat('C:\Users\ckohl\hnn_out\data\default_Sarah\dpl.txt');
end
%load
dpl = readtable(file_dir);
dpl = table2array(dpl);
line(end+1)=plot(dpl(:,1),dpl(:,2), 'Color','k', 'Linewidth',2);
leg_name{end+1}='X';
legend(line,leg_name)
    
xlabel('Amplitude (nAm)')
ylabel('Time (ms)')
if supra
    ylim([-210 200])
else
    ylim([-150 100])
end
xlim([0 270])
