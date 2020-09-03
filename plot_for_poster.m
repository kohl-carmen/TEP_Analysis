%plot fot poster

sim_path='C:\Users\ckohl\hnn_out\data';
data_path='C:\Users\ckohl\Desktop\Current\TEP\mine_temp\';

sim_name='Sarah_defaultx3_plus10';
data_name='AVG_FOR_HNN_FROM_DIGITISER';

%load
sim=load(strcat(sim_path,'\',sim_name,'\dpl.txt'));

sim_time=sim(:,1);
sim=sim(:,2);

spk=load(strcat(sim_path,'\',sim_name,'\spk.txt'));


data=load(strcat(data_path,data_name,'.txt'));

data_time=data(:,1);
data=data(:,2);

%plot stuff
data_colour='k';
sim_colour='b';

clf
hold on
lines=[];
lines(1)=plot(data_time,data,'Linewidth',2,'Color',data_colour);
lines(2)=plot(sim_time,sim,'Linewidth',2,'Color',sim_colour);

% %tms lines
% ylims=ylim();
% plot([0,0],ylims,'--','Linewidth',3,'Color',[.5 .5 .5])
% plot([20,20],ylims,'--','Linewidth',3,'Color',[.5 .5 .5])
% plot([40,40],ylims,'--','Linewidth',3,'Color',[.5 .5 .5])

%inputs
p_colour=['r'];
d_colour='g';
tms_colour=[.5 .5 .5];
p=[37,57,77,147,167,187];
d=[74,94,114];
tms=[0,20,40];

plot(p,ylims(1),'o','Color',p_colour)
plot(d,ylims(2),'o','Color',d_colour)
plot(tms,ylims(2),'o','Color',tms_colour)

legend(lines,{'TEP','Simulation'})


cd('C:\Users\ckohl\Desktop\Current\TEP\')
print -depsc  spk_Temp1


%spik
fig1=figure(1);
clf
hold on
fig1.Renderer='Painters';
spk_times=spk(spk(:,2)<271,1);
spk_ids=spk(spk(:,2)<271,2);
hold on
cell_types = {'L2_basket', 'L2_pyramidal', 'L5_basket', 'L5_pyramidal'};
cell_colours={[.5 .75 1],[0 0 .75],[1 .5 0],[.75 0 0]};
cell_cutoff=[0 35 135 170 270]  ;  
for cell=1:length(cell_types)
    these_spikes=[];
    these_times=[];
    for id=1:length(spk_ids)
        if spk_ids(id)>cell_cutoff(cell) & spk_ids(id)<=cell_cutoff(cell+1)
            these_times=[these_times,spk_times(id)];
            these_spikes=[these_spikes,spk_ids(id)];
        end
    end
    plot(these_times,these_spikes,'.','Color',cell_colours{cell})
end

