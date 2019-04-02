function S = sr3_subj(subj, fig, block, trial)
%% function S = sr3_subj(subj, fig, block, trial)
% Subject routine for SequenceRepetition experiment 3
%
% Example calls:
%               S=sr3_subj('s01');         % which subj
%               S=sr3_subj('s01',0);       % which subj, with/without plot
%               S=sr3_subj('s01',0,3);     % which subj, with/without plot, which block, all trials
%               S=sr3_subj('s09',1,7,5);   % which subj, with/without plot, which block, and which trial in this block
%
%%
if nargin<2
    fig=0; %don't produce a plot
end;
%%
pathToData='/Users/gariani/Documents/data/SequenceRepetition/sr3'; %path to data
datafilename=fullfile(pathToData,sprintf('sr3_%s.dat',subj)); %input
outfilename=fullfile(pathToData,sprintf('sr3_%s.mat',subj)); %output
%%
S=[]; %preallocate an empty output structure
D=dload(datafilename); %load dataset for this subj
if (nargin<3)
    trials=1:numel(D.TN); %analyze all trials for all blocks
    block=-1; %initialize block count variable
else
    if (nargin<4)
        trials=find(D.BN==block & D.TN==1):find(D.BN==block & D.TN==numel(unique(D.TN))); %analyze all trials for this block
    else
        trials=find(D.BN==block & D.TN==trial); %analyze this trial of this block
    end;
end;
%%
for t=trials
    
    if ~exist('MOV','var') || (D.BN(t)~=block)
        block=D.BN(t); %update block number within the loop
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% specific to sr3! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(subj, 's09') && block == 8 || (strcmp(subj, 's40') && block == 2) % weird MOV files, do not match with .dat files
        else
            MOV=movload(fullfile(pathToData,sprintf('sr3_%s_%02d.mov',subj,block))); %load MOV file for this block
        end
    end;
    
    fprintf(1,'\nsubject: %s   block: %02d   trial: %02d\n\n',subj,block,D.TN(t));
    fig_name=sprintf('sr3_%s_b%02d_t%02d',subj,block,D.TN(t));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% specific to sr3! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(subj, 's09') && block == 8 || (strcmp(subj, 's40') && block == 2) % weird MOV files, do not match with .dat files
        T = getrow(D, t);
        fNames = fieldnames(T);
        resp = [];
        for i = 1:length(fNames)
            if length(fNames{i}) > 8
                if strcmp(fNames{i}(1:8), 'response')
                    eval(['resp=[resp, T.', fNames{i}, '];']);
                end
            end
        end
        T.numPress = sum(resp>0); % add info about number of presses
    else
        T=sr3_trial(MOV{D.TN(t)},getrow(D,t),fig,fig_name);
    end
    
    S=addstruct(S,T,'row','force');
end
%%
if (nargin<3)
    save(outfilename,'-struct','S');
end;
