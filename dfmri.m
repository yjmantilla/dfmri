% conn_batch_*.m files provide examples for batch scripting in CONN
% For better documentation read conn_batch.m
% https://gist.github.com/A2ed/3114405e2eb216f53a6b has some errors ?
clc;
clear all;
clear BATCH;

data_path = 'Z:\dfmri\SAN\Rest';
MyFolderInfo = dir(data_path);
files = {MyFolderInfo.name}; % i hate matlab... 
file_names = {};
j = 1;
bytes = {MyFolderInfo.bytes};

for i = 1:length(files)
    if bytes{i}> 0
        file_names{j} = files{i};
        j = j + 1;
    end
end

conn_path = 'Y:/code/spm12';
spm_path = 'Y:/code/conn';

addpath(genpath(conn_path));
addpath(genpath(spm_path));

root_path = 'Z:/dfmri';
sitid = 'KKI';
tr = 2.5;
vol = 152;

%data_path = 

% we could create the list of filenames algorithmically --> actually this
% is fundamental for the workflow
nsubs = 3;
nses = 1;
san = {'_KKI_1779922_15.nii','_KKI_2903997_44.nii','_KKI_3699991_57.nii'};
san_ages = [11.04,11.95,12.98];

tdah = {'_KKI_1019436_2.nii','_KKI_1996183_21.nii','_KKI_3160561_50.nii'};
tdah_ages = [10.84,11.91,12.77];
san_path = 'SAN';
tdah_path = 'TDAH';

selected_files = san;
selected_path = san_path;

rest_path = 'Rest/swaurest';
anat_path = 'Anat/wc3canat';


name = 'dfmri.mat';
%Set file/dir to write
BATCH.filename = [root_path '/' name];

%Set number of subjects
BATCH.Setup.nsubjects = nsubs;
BATCH.Setup.nsessions = nses*ones(1,BATCH.Setup.nsubjects);

%Set TR value
BATCH.Setup.RT = tr.*ones(1,BATCH.Setup.nsubjects);

%Set acquisition type
BATCH.Setup.acquisitiontype = 1;

conn_batch(BATCH);

%%% NOTE BATCH.New for spatial pre-processing so it is not needed here %%%

%Structural
for sub = 1:BATCH.Setup.nsubjects
    % for default structural conn file:
    filename = 'Y:/code/conn/utils/surf/referenceT1.nii';
    
    % Anatomical
    %filename = [root_path '/' selected_path '/' anat_path selected_files{i}];
    BATCH.Setup.structurals{sub} = filename;
    clear filename;
end

conn_batch(BATCH);

%Functional

for sub = 1:BATCH.Setup.nsubjects
        for ses = 1:BATCH.Setup.nsessions(sub)
            filename = [root_path '/' selected_path '/' rest_path selected_files{sub}];
            BATCH.Setup.functionals{sub}{ses} = filename;
            clear filename;
        end
end

conn_batch(BATCH);
% Temporal Decomposition (Sliding Window)

window_length = 60;
overlap = 0.5;
onsets = 0:window_length*overlap:tr*vol; % should be length of something here?
nconditions = length(onsets);
BATCH.Setup.conditions.names=[{'rest'}, arrayfun(@(n)sprintf('rest x Time%d',n),1:nconditions,'uni',0)];
for sub = 1:BATCH.Setup.nsubjects %{ncondition}{nsub}{nses}
    for ses = 1:BATCH.Setup.nsessions(sub)
            BATCH.Setup.conditions.onsets{1}{sub}{ses}=[0];
            BATCH.Setup.conditions.durations{1}{sub}{ses}=[inf];
        for win = 1:length(onsets)
            BATCH.Setup.conditions.onsets{win+1}{sub}{ses}=[onsets(win)];
            BATCH.Setup.conditions.durations{win+1}{sub}{ses}=[window_length];
        end
    end

end

conn_batch(BATCH);

BATCH.Setup.analyses=[1,2,3,4];
BATCH.Setup.outputfiles=[1,1,1,1,1,1];
BATCH.Setup.done=1;
BATCH.Setup.overwrite='No';
conn_batch(BATCH);

BATCH.Preprocessing.done=0;
BATCH.Preprocessing.overwrite='No';
conn_batch(BATCH);

BATCH.Denoising.done = 1;
BATCH.Denoising.overwrite = 'No';
BATCH.Denoising.filter = [0.008,0.09];
BATCH.Denoising.detrending = 1;
BATCH.Denoising.regbp = 1;
BATCH.Denoising.despiking = 0;
BATCH.Denoising.confounds.names = {'White Matter','CSF','Effect of rest'};
conn_batch(BATCH);


% Crear todas las metricas de conectividad funcional

%conn
%conn('load',fullfile(root_path,name));
%conn gui_results
% do note that rest (in effect of rest) here could be any condition

% %Set the type of analyses to run
% BATCH.Setup.analyses=[1,2,3];
% 
% %Set voxel resolution to 2mm
% BATCH.Setup.voxelresolution = 1;
% 
% %Set BOLD signal units to percent signal change
% BATCH.Setup.analysisunits = 1;
% 
% 
% %Loop through files and associate runs and conditions
% for i = 1:18
%         for j = 1:6
%         filename = sprintf('/Users/Aiden/Imaging_data/EM_Meta/Connectivity/SpatialWM/Preprocessed_data/s0%d_run%d_filtered_func_data.nii', i, j);
%         BATCH.Setup.functionals{i}{j} = filename;
%         clear filename;
%     end
% end
% 
% %Loop through files and associate structurals
% for i = 1:18
%     filename = sprintf('/Users/Aiden/Imaging_data/EM_Meta/Connectivity/SpatialWM/Preprocessed_data/mprage_brain_s0%d.nii', i);
%     BATCH.Setup.functionals{i} = filename;
%     clear filename;
% end
% 
% %Create lists to hold condition names and cue type
% conditionlist = ['correct' ; 'INC    '];
% typelist = ['Cue  ' ; 'Probe'];
% 
% %Specify conditions. These are indexed as 1, 2, 3, 4
% BATCH.Setup.conditions.names = {'cue_corr', 'cue_inc', 'probe_corr', 'probe_inc'};
% 
% %Loop through onset files to insert based on condition
% 
% %Loop through subjects
% for i = 1:18
%     
%     %Loop through runs
%     for j = 1:6
%         
%         %Loop through conditions
%         for k = 1:4
%             
%             %read onset vector from text file
%             filename = sprintf('/Users/Aiden/Imaging_data/EM_Meta/Connectivity/SpatialWM/conn_onsets/no_condition/%d_%d_%d', i, j, k);
%             f = dlmread(filename);
%             
%             %write onset vector into batch structure
%             BATCH.Setup.conditions.onsets{i}{j}{k} = [f];
%             
%             %write duration time
%             BATCH.Setup.conditions.durations{i}{j}{k} = [3];
%             
%             clear filename;
%             
%         end
%     end
% end
% 
% conn_batch(BATCH);