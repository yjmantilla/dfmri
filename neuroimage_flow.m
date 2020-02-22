% conn_batch_*.m files provide examples for batch scripting in CONN
% For better documentation read conn_batch.m
% https://gist.github.com/A2ed/3114405e2eb216f53a6b has some errors ?
clc;
clear all;
clear BATCH;

% Library Importing
conn_path = 'Y:/code/spm12';
spm_path = 'Y:/code/conn';

addpath(genpath(conn_path));
addpath(genpath(spm_path));


% Data where files are located
data_path = 'Z:\dfmri\SAN\Rest\';

% Generate filelist
MyFolderInfo = dir(data_path);
files = {MyFolderInfo.name}; % i hate matlab..., this is to be able to index the names
file_names = {};
j = 1;

bytes = {MyFolderInfo.bytes}; % for checking if any file is empty...

for i = 1:length(files)
    if bytes{i}> 0 
        file_names{j} = files{i};
        j = j + 1;
    end
end


% Path to save the project
root_path = 'Z:/dfmri';

% Dataset characteristics
tr = 2.5;
vol = 152;
nsubs = length(file_names); % assume one file per subject
nses = 1; % assume only one session

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
    
    % Using Anatomical Files
    %filename = [anat_path selected_files{i}];
    BATCH.Setup.structurals{sub} = filename;
    clear filename;
end

conn_batch(BATCH);

%Functional

for sub = 1:BATCH.Setup.nsubjects
        for ses = 1:BATCH.Setup.nsessions(sub)
            filename = [data_path file_names{sub}];
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
        for win = 1:length(onsets) % for each window create a new condition
                                   % win + 1 since 1st condition is continuous rest 
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

