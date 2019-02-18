function [varargout] = sr3_makeTargetFiles(varargin)
% function [varargout] = sr3_makeTargetFiles(varagin)
% Creates .tgt files for subject(s) s and block(s) b
%
% example calls:
%               [G] = sr3_makeTargetFiles([99], [1:12]);
%               [G] = sr3_makeTargetFiles([1:20], [1:12]);
%
% --
% gariani@uwo.ca - 2018.04.05

%% make target files for all subjects and blocks per subject
G = struct();
for s = varargin{1}
    S = struct();
    for b = varargin{2}
        fprintf(1, '\nsubj: %d   block: %d\n', s, b);
        [B] = sr3_target(s, b); % B=block (all trials)
        %fprintf(1, '\nGo-Go Repetitions: %02d trials (%03d%%)   Nogo-Go Repetitions: %02d trials (%02d%%)',     B.GoGoRep,      round(B.GoGoRep_per*100),       B.NogoGoRep,    round(B.NogoGoRep_per*100));
        %fprintf(1, '       Go-Go Switch: %02d trials (%03d%%)        Nogo-Go Switch: %02d trials (%02d%%)\n',   B.GoGoSwitch,   round(B.GoGoSwitch_per*100),    B.NogoGoSwitch, round(B.NogoGoSwitch_per*100));
        S = addstruct(S, B); % S=subject (all blocks)
    end
    fprintf(1, '\nGo-Go Repetitions: %03d trials (%02d%%)   Nogo-Go Repetitions: %02d trials (%02d%%)',     sum(S.GoGoRep),      round(mean(S.GoGoRep_per*100)),       sum(S.NogoGoRep),    round(mean(S.NogoGoRep_per*100)));
    fprintf(1, '       Go-Go Switch: %03d trials (%02d%%)        Nogo-Go Switch: %02d trials (%02d%%)\n',   sum(S.GoGoSwitch),   round(mean(S.GoGoSwitch_per*100)),    sum(S.NogoGoSwitch), round(mean(S.NogoGoSwitch_per*100)));
    G = addstruct(G, S); % G=group (all subjects)
end
%fprintf(1, '\nGo-Go Repetitions: %04d trials (%02d%%)   Nogo-Go Repetitions: %03d trials (%02d%%)',     sum(G.GoGoRep),      round(mean(G.GoGoRep_per*100)),       sum(G.NogoGoRep),    round(mean(G.NogoGoRep_per*100)));
%fprintf(1, '       Go-Go Switch: %04d trials (%02d%%)        Nogo-Go Switch: %03d trials (%02d%%)\n',   sum(G.GoGoSwitch),   round(mean(G.GoGoSwitch_per*100)),    sum(G.NogoGoSwitch), round(mean(G.NogoGoSwitch_per*100)));
varargout{1} = G;
end

function [varargout] = sr3_target(s,b)
% function [varargout] = sr3_target(s,b)
% This function generates .tgt files, one per block
%
% inputs: vector of subject numbers (s), vector of block numbers (b)
% output: saved filename (fn), block structure (B)

%% define target folder
targetFolder = '../../../../robotcode/projects/SequenceRepetition/sr3/target'; %save tgt files in the right relative path
if ~exist(targetFolder,'dir'); mkdir(targetFolder); end %create target folder if it doesn't already exist

%% experimental details
%allseq = perms([1,2,3,5]); % all possible sequences with each of the four fingers used only once
allseq = permn([1,2,3,4,5], 4); % all possible sequences taken 4 at a time, with repetitions
allseq_noRun = zeros(1,size(allseq,2)); % remove "runs" (e.g. 1,2,3, ... 3,2,1, ...)
rowCount = 1;
for r = 1:size(allseq, 1)
    
    isRun = zeros(1,size(allseq, 2) - 2);
    for c = 1:size(allseq, 2) - 2
        if allseq(r, c)==(allseq(r, c+1) - 1) && allseq(r, c)==(allseq(r, c+2) - 2) % flag descending runs
            isRun(1, c) = 1;
        elseif allseq(r, c)==(allseq(r, c+1) + 1) && allseq(r,c)==(allseq(r, c+2) + 2) % flag ascending runs
            isRun(1, c) = 1;
        elseif any(allseq(r, c)==(allseq(r, c+1:end))) || any(allseq(r, c+1)==(allseq(r, c+2:end))) % flag same-number repetitions
            isRun(1, c) = 1;
        else
            isRun(1, c) = 0;
        end
    end
    
    if all(isRun==0)
        allseq_noRun(rowCount,:) = allseq(r,:); % pool of similarly difficult sequences to pick from
        rowCount = rowCount + 1;
    else
        continue
    end
    
end
randind = randperm( size(allseq_noRun, 1)); % randomize pool of sequences
allseq_noRun = allseq_noRun(randind, :); % select random sequences from pool
%%%
nSeq = 4; % how many sequences do you want to pick?
%%%
seq = allseq_noRun(1:nSeq, :); % pick nSeq number of sequences
seqNum = randind(1:nSeq)'; % which sequences from pool?
maxTrials = 50; % roughly how many trials per block?
trials = floor(maxTrials/nSeq) * nSeq; % exactly how many trials per block
p_rep = 0.50; % probability of same sequence repetition
p_switch = (1 - p_rep) / (nSeq-1); % probability of switch to each other sequence
T = eye(nSeq); T(T==1) = p_rep; T(T==0) = p_switch; % transition matrix
p_nogo = 0.30; % nogo probability
nogo_trials = round(p_nogo * trials); % proportion of nogo trials

%% fill in dataframe structure B for this block
%%%
exeType = ones(trials, 1); exeType(1:nogo_trials) = 0; exeType = exeType(randperm(trials)); % 1=go; 0=nogo;
%%%
p1 = ones(1,nSeq) * (1/nSeq); % initial probablity (the same for all sequences)
z = pickRandSeq(p1); % initialize the first random pick
for t = 1 : trials
    B.SN(t,1) = s; % subject number
    B.seqNum(t,1) = seqNum(find(z,1)); % which sequence?
    
    if (t > 1) && (B.seqNum(t, 1) == B.seqNum(t-1, 1)) % check if this is a repetition or not
        B.isRep(t,1) = 1; % is this a repetition? (0=no; 1=yes)
        B.repNum(t,1) = B.repNum(t-1, 1) + 1; % what repetition number is this?
        B.exeType(t,1) = exeType(t); % is this a go or nogo trial? (1=go; 0=nogo)
        
        if B.exeType(t-1, 1) == 1 % check if this trial was preceded by a go
            % repType explanation:
            % the first number corresponds to whether the current trial n is a repetition (1) or not (0),
            % the second number corresponds to whether the n-1 trial was a go (1) or a nogo (0),
            % the third number corresponds to whether the current trial n is a go (1) or a nogo (0)
            B.repType(t,1) = str2double( num2str([B.isRep(t,1), B.exeType(t-1,1), B.exeType(t,1)], '%d%d%d')); % repetition(1), preceded by a go(1), current go(1) or nogo(0)
            %B.repType{t,1} = num2str([B.isRep(t,1), B.exeType(t-1,1), B.exeType(t,1)], '%d%d%d'); % repetition(1), preceded by a go(1), current go(1) or nogo(0)
        else % previous trial was a nogo
            B.repType(t,1) = str2double( num2str([B.isRep(t,1), B.exeType(t-1,1), B.exeType(t,1)], '%d%d%d')); % repetition(1), preceded by a nogo(0), current go(1) or nogo(0)
            %B.repType{t,1} = num2str([B.isRep(t,1), B.exeType(t-1,1), B.exeType(t,1)], '%d%d%d'); % repetition(1), preceded by a nogo(0), current go(1) or nogo(0)
        end
        
    else % switch trial
        B.isRep(t,1) = 0; % is this a repetition? (0=no; 1=yes)
        B.repNum(t,1) = 0; % what repetition number is this?
        B.exeType(t,1) = exeType(t); % is this a go or nogo trial? (1=go; 0=nogo)
        
        if (t > 1) && (B.exeType(t-1, 1) == 1) % check if this trial was preceded by a go
            % repType explanation:
            % the first number corresponds to whether the current trial n is a repetition (1) or not (0),
            % the second number corresponds to whether the n-1 trial was a go (1) or a nogo (0),
            % the third number corresponds to whether the current trial n is a go (1) or a nogo (0)
            B.repType(t,1) = str2double( num2str([B.isRep(t,1), B.exeType(t-1,1), B.exeType(t,1)], '%d%d%d')); % switch(0), preceded by a go(1), current go(1) or nogo(0);
            %B.repType{t,1} = num2str([B.isRep(t,1), B.exeType(t-1,1), B.exeType(t,1)], '%d%d%d'); % switch(0), preceded by a go(1), current go(1) or nogo(0);
        else % previous trial was a nogo
            if t == 1 % the first trial is by definition always preceded by a nogo
                B.repType(t,1) = str2double( num2str([B.isRep(t,1), 0, B.exeType(t,1)], '%d%d%d')); % switch(0), preceded by a nogo(0), current go(1) or nogo(0);
                %B.repType{t,1} = num2str([B.isRep(t,1), 0, B.exeType(t,1)], '%d%d%d'); % switch(0), preceded by a nogo(0), current go(1) or nogo(0);
            else
                B.repType(t,1) = str2double( num2str([B.isRep(t,1), B.exeType(t-1,1), B.exeType(t,1)], '%d%d%d')); % switch(0), preceded by a nogo(0), current go(1) or nogo(0);
                %B.repType{t,1} = num2str([B.isRep(t,1), B.exeType(t-1,1), B.exeType(t,1)], '%d%d%d'); % switch(0), preceded by a nogo(0), current go(1) or nogo(0);
            end
        end
        
    end
    
    B.prepTime(t,1) = 2500; % fixed preparation time (in ms)
    B.iti(t,1) = 500; % fixed inter-trial-interval duration (in ms)
    B.feedback(t,1) = 1; % give feedback on go trials only (yes/no)
    B.hand(t,1) = 2; % right hand (left=1/right=2)
    B.cuePress{t,1} = num2str(seq(find(z,1), :), '%d%d%d%d'); % which cue press?
    B.cueMask{t,1} = char('****'); % mask the cue press after go/nogo
    B.press1{t,1} = B.cuePress{t,1}(1); % 1st finger press
    B.press2{t,1} = B.cuePress{t,1}(2); % 2nd finger press
    B.press3{t,1} = B.cuePress{t,1}(3); % 3rd finger press
    B.press4{t,1} = B.cuePress{t,1}(4); % 4th finger press
    
    p = T * z; % update probabilities
    z = pickRandSeq(p); % update pick
end

%% save structure B as a target file (.tgt)
filename = fullfile(targetFolder, sprintf('sr3_s%02d_b%02d.tgt', s, b));
dsave(filename, B);

%% check how many repetitions per each condition of interest
B.GoGoRep = sum(B.repType == 111);
B.NogoGoRep = sum(B.repType == 101);
B.GoGoSwitch = sum(B.repType == 011);
B.NogoGoSwitch = sum(B.repType == 001);

B.GoGoRep_per = sum(B.repType == 111) / trials;
B.NogoGoRep_per = sum(B.repType == 101) / trials;
B.GoGoSwitch_per = sum(B.repType == 011) / trials;
B.NogoGoSwitch_per = sum(B.repType == 001) / trials;

%% return output data structure B
B.fn = filename;
varargout{1} = B;
end

function [z] = pickRandSeq(p)
% function [z] = pickRandSeq(p)
%
% input: p, probability vector of length nSeq (probability of each sequence being selected)
% output: z, indicator variable of length nSeq (which sequence has been selected: 1=yes, 0=no)

%% make p a column vector
p = p(:);

%% initialize indicator variable z
z = zeros(numel(p), 1);

%% define cumulative probability
cump = cumsum(p);

%% pick random number from continuous uniform distribution with range [0 1]
r = unifrnd(0,1);

%% check in which class the number falls
i = find(r <= cump, 1);

%% select and return the appropriate class
z(i) = 1;
end