% BrainSight: Get stimulation coordinates to use in mne later
% See Import_BrainSight.m for info
clear
file_dir='C:\Users\ckohl\Desktop\Current\Data\TMS\'
file_name='Exported Brainsight Data_MNI'
% mri_dir='C:\Users\ckohl\Documents\Virtual_Shared\Pilot\';
partic='02';
write_to_file=1;
output_dir=file_dir;

%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%% Get Data
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
fprintf('\n-----------------\n')
fprintf('Load Data')
fprintf('\n-----------------\n')
%First, we'll import the first column to see how many samples there are
opts = delimitedTextImportOptions("NumVariables", 33);
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["SampleName", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33"];
opts.SelectedVariableNames = "SampleName";
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
ExportedBrainsightData = readtable(strcat(file_dir,file_name), opts);
firstcolumn_str = table2array(ExportedBrainsightData);

samples_start_row=[];
samples_end_row=[];
targets_start_row=[];
targets_end_row=[];
target_names={};
% target_count=struct();
for row=1:size(firstcolumn_str,1)
    if ~isempty(targets_start_row) & isempty(targets_end_row)
        if firstcolumn_str{row}(1) ~='#'
            target_names{end+1}=firstcolumn_str{row};
%             target_count.(firstcolumn_str{row}(find(~isspace(firstcolumn_str{row}))))=0;
        end
    end
    if length(firstcolumn_str{row})> 7
        if firstcolumn_str{row}(1:8)=='# Target'
            targets_start_row=[targets_start_row,row+1];
        end
        if all(firstcolumn_str{row}(1:8)=='# Sample') & all(firstcolumn_str{row+1}(1:7)=='Sample ')
            samples_start_row=[samples_start_row,row+1];
            targets_end_row=[targets_end_row,row-1];
        end
        if firstcolumn_str{row}(1:6)=='Sample' & firstcolumn_str{row+1}(1)=='#'
            samples_end_row=[samples_end_row,row];
        end
    end
%     if row>=samples_start_row & (isempty(samples_end_row) | row==samples_end_row)
%         target_count.(AssocTarget{row-samples_start_row+1}(find(~isspace(AssocTarget{row-samples_start_row+1}))))=target_count.(AssocTarget{row-samples_start_row+1}(find(~isspace(AssocTarget{row-samples_start_row+1}))))+1;
%     end
end
if isempty(samples_start_row) | isempty(samples_end_row)
    fprintf('Target Detection failed: No Sample Start/End detected \n');
else
    fprintf('Target Detection successful \n');
    fprintf('%i Targets detected: \n',length(targets_start_row:targets_end_row));
    for t=1:length(target_names)
        fprintf('\t %s \n', target_names{t});
    end
end
if isempty(samples_start_row) | isempty(samples_end_row)
    fprintf('Sample Detection failed: No Sample Start/End detected \n')
elseif length(samples_start_row)>1 | length(samples_end_row)>1
    fprintf('Sample Detection failed: More than one Sample Start/End detected \n')
elseif samples_start_row >= samples_end_row
    fprintf('Sample Detection failed: Sample Start/End makes no sense \n')
else
    fprintf('Sample Detection successful \n')
    fprintf('%i Samples detected \n',samples_end_row-samples_start_row)
end

% get column headers
opts = delimitedTextImportOptions("NumVariables", 33);
opts.DataLines = [samples_start_row-1, samples_start_row-1];
opts.Delimiter = "\t";
ExportedBrainsightData = readtable(strcat(file_dir,file_name), opts);
column_headers = table2array(ExportedBrainsightData);
%chnage names so they can be valid variable names
for header=1:length(column_headers)
    column_headers{header}(column_headers{header}==' ')=[];
    column_headers{header}(column_headers{header}=='.')=[];
    column_headers{header}(column_headers{header}=='#')=[];
    column_headers{header}(column_headers{header}=='-')=[];    
end

% get data (loop through columns)
opts.VariableNames = column_headers;
for columns=1:length(column_headers)
    if columns>1 & columns <14
        opts.DataLines = [targets_start_row, targets_end_row];
        opts.SelectedVariableNames=column_headers{columns};
        tbls = readtable(strcat(file_dir,file_name), opts); 
        tbls = table2array(tbls);
        try
            eval(strcat('Target_',column_headers{columns+3}, ' =  double(tbls);'))
        catch
            for i=1:length(tbls)
                mat(i,1)=str2num(tbls{i});
            end
            eval(strcat('Target_',column_headers{columns+3}, ' =  mat;'))
        end
    end
    opts.DataLines = [samples_start_row, samples_end_row];
    opts.SelectedVariableNames=column_headers{columns};
    tbls = readtable(strcat(file_dir,file_name), opts); 
    tbls = table2array(tbls);
    if any(columns==[ 3,5:20,28:32])
        try
         eval(strcat(column_headers{columns}, ' =  double(tbls);'))
        catch
            for i=1:length(tbls)
                if isempty(str2num(tbls{i}))
                   mat(i,1)=nan;
                else
                    mat(i,1)=str2num(tbls{i});
                end
            end
            eval(strcat(column_headers{columns}, ' =  mat;'))
        end
    else
        eval(strcat(column_headers{columns}, ' =  tbls;'))
    end      
end

% Find samples which were part of triple-pulses 
Pulse_ID=nan(length(Time),1) ;%0= single, 1 2 3 triplets
Duration_since_last_pulse=nan(length(Time),1) ;
for sample=2:length(Time)

    Y=str2num(Date{sample}(1:4));
    M=str2num(Date{sample}(6:7));
    D=str2num(Date{sample}(9:10));
    
    h=str2num(Time{sample}(1:2));
    MI=str2num(Time{sample}(4:5));
    S=str2num(Time{sample}(7:8));
    MS=str2num(Time{sample}(10:end));
    
    sample_t = datetime(Y,M,D,h,MI,S,MS);
    
    Y=str2num(Date{sample-1}(1:4));
    M=str2num(Date{sample-1}(6:7));
    D=str2num(Date{sample-1}(9:10));
    
    h=str2num(Time{sample-1}(1:2));
    MI=str2num(Time{sample-1}(4:5));
    S=str2num(Time{sample-1}(7:8));
    MS=str2num(Time{sample-1}(10:end));
    
    prev_t = datetime(Y,M,D,h,MI,S,MS);
    
    Duration_since_last_pulse(sample)=milliseconds(sample_t-prev_t);
end
Temp=[Duration_since_last_pulse;999];
Duration_since_last_pulse_temp=[999;Duration_since_last_pulse];
Pulse_ID(Duration_since_last_pulse_temp>40 & Temp > 40)=0;   
Pulse_ID(Duration_since_last_pulse_temp>40 & Temp < 40)=1;
Pulse_ID(Duration_since_last_pulse_temp<40 & Temp < 40)=2;   
Pulse_ID(Duration_since_last_pulse_temp<40 & Temp > 40)=3;  
Pulse_ID=Pulse_ID(2:end);
if isnan(Pulse_ID(1))
    if Pulse_ID(2)==2 & Pulse_ID(3)==3
        Pulse_ID(1)=1;
    else
        Pulse_ID(1)=0;
    end
end
fprintf('Number of single pulses detected: %i \n',sum(Pulse_ID==0))
fprintf('Number of triplets detected: %i \n',sum(Pulse_ID==1))

% Get cont EMG in usable format
EMGData=[];
for trial=1:length(EMGData1)
    if EMGData1{trial}(1:6)~='(null)'
        EMGData(:,trial)=eval(strcat('[',EMGData1{trial},']'));
    end
end
EMGData=EMGData';

    
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------    
%% Define target
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
fprintf('\n-----------------\n')
fprintf('Define Targets')
fprintf('\n-----------------\n')
% % First, check if all triple-pulses had the same target
% targets_detected={};
% for trial=1:length(AssocTarget)
%     if Pulse_ID(trial)>0
%         targets_detected{end+1}=AssocTarget{trial};
%     end
% end
% First, check if all -pulses had the same target
%( this doesn't got through all the targtes listed, it goes through all the targets used)
targets_detected={};
for trial=1:length(AssocTarget)
    targets_detected{end+1}=AssocTarget{trial};
end
[targets_detected,temp,target_count]=unique(targets_detected);
if length(targets_detected)>1
    fprintf('More than one target detected during triple pulses:\n')
    for i=1:length(targets_detected)
        fprintf('\t(%d) %s (%d) \n',i,targets_detected{i},sum(target_count==i));
    end
    chosen_target=input('Pick target: ');
    target_oi=targets_detected{chosen_target};
else
    fprintf('Target:\n\t%s\n',targets_detected{1})
    target_oi=targets_detected{1};
end

for i=1:length(target_names)
    if length(target_names{i}(~isspace(target_names{i})))==length(target_oi(~isspace(target_oi)))
        if all(target_names{i}(~isspace(target_names{i}))==target_oi(~isspace(target_oi)))
            target_i=i;
        end
    end
end
T_X=Target_LocX(target_i);
T_Y=Target_LocY(target_i);
T_Z=Target_LocZ(target_i);

Target_vector=zeros(size(AssocTarget));
for trial=1:length(AssocTarget)
    if length(AssocTarget{trial})==length(target_oi)
        if AssocTarget{trial}== target_oi
            Target_vector(trial)=1;
        end
    end
end

cont=input('Press any key to continue');

target_coords=[T_X,T_Y,T_Z];
sample_coords=[LocX, LocY, LocZ];
sample_coords=sample_coords(Target_vector==1,:);

%Output file:
if write_to_file
    filename=strcat(output_dir,'BETA',num2str(partic),'_StimLocs_',target_oi(end-2:end),'.txt');
    output_file=fopen(filename,'a');
    fprintf(output_file,'%f\t%f\t%f\n',target_coords);
    for smp= 1:size(sample_coords,1)
        fprintf(output_file,'%f\t%f\t%f\n',sample_coords(smp,:));
    end   
    fprintf('\nSaved %d sample coords in %s\n',smp,filename)
end
