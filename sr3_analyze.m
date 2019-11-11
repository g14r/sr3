function [varargout] = sr3_analyze(what, varargin)
%% function [varargout] = sr3_analyze(what, varargin)
% SequenceRepetition experiment 3: analysis of behavioral data
%
% usage|example calls:
%
%                       sr3_analyze('all_subj');                                %pre-analysis: run the subject routine for all_subj
%                       sr3_analyze('all_subj', {'s01'});                       %pre-analysis: run the subject routine for selected subjects
%                       [all_data] = sr3_analyze('all_data');                   %pre-analysis: create .mat file containing data from subjects in subj
%
%                       [D] = sr3_analyze('repEff_ET');                         %group          analysis of repetition effect on execution time (ET)
%                       [D] = sr3_analyze('repEff_ET', 's01');                  %single subject analysis of repetition effect on execution time (ET)
%
%                       [D] = sr3_analyze('N-1_ET');                            %group          analysis of N-1 trial go/no-go on repetition effect for ET
%                       [D] = sr3_analyze('N-1_ET', 's01');                     %single subject analysis of N-1 trial go/no-go on repetition effect for ET
%
%                       [D] = sr3_analyze('repEff_RT');                         %group          analysis of repetition effect on execution time (RT)
%                       [D] = sr3_analyze('repEff_RT', 's01');                  %single subject analysis of repetition effect on execution time (RT)
%
%                       [D] = sr3_analyze('N-1_RT');                            %group          analysis of N-1 trial (go/no-go) on repetition effect for RT
%                       [D] = sr3_analyze('N-1_RT', 's01');                     %single subject analysis of N-1 trial (go/no-go) on repetition effect for RT
%
%                       [D] = sr3_analyze('firstFinger_ET');                    %group          analysis of same first finger repetition (different sequence) influence on ET (vs same sequence)
%                       [D] = sr3_analyze('firstFinger_ET', 's01');             %single subject analysis of same first finger repetition (different sequence) influence on ET (vs same sequence)
%
%                       [D] = sr3_analyze('firstFinger_RT');                    %group          analysis of same first finger repetition (different sequence) influence on RT (vs same sequence)
%                       [D] = sr3_analyze('firstFinger_RT', 's01');             %single subject analysis of same first finger repetition (different sequence) influence on RT (vs same sequence)
%
%                       [D] = sr3_analyze('slf_RT');                            %group          analysis of same last and first finger of next seq (different sequence) influence on RT
%                       [D] = sr3_analyze('slf_RT', 's01');                     %single subject analysis of same last and first finger of next seq (different sequence) influence on RT
% 
%                       [D] = sr3_analyze('FAp_RT');                            %group          analysis of False Alarm probability influence on RT (no-go trials only)
%                       [D] = sr3_analyze('FAp_RT', 's01');                     %single subject analysis of False Alarm probability influence on RT (no-go trials only)
%
%                       [D] = sr3_analyze('ET_RT');                             %group          analysis of relationship between ETs and RTs (correlation)
%                       [D] = sr3_analyze('ET_RT', 's01');                      %single subject analysis of relationship between ETs and RTs (correlation)
%
%                       [D] = sr3_analyze('first2IPIs_RT');                     %group          analysis of relationship between first 2 IPIs and RTs (correlation)
%                       [D] = sr3_analyze('first2IPIs_RT', 's01');              %single subject analysis of relationship between first 2 IPIs and RTs (correlation)
%
%                       [D] = sr3_analyze('N-1_IPI');                           %group          analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs
%                       [D] = sr3_analyze('N-1_IPI', 's01');                    %single subject analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs
%
%                       [D] = sr3_analyze('repEff_N-1_smt');                    %group          analysis of N-1 trial go/no-go on repetition effect, only sequences with same middle transition (smt)
%                       [D] = sr3_analyze('repEff_N-1_smt', 's01');             %single subject analysis of N-1 trial go/no-go on repetition effect, only sequences with same middle transition (smt)
%
%                       [D] = sr3_analyze('N-1_IPI_smt');                       %group          analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, only sequences with same middle transition (smt)
%                       [D] = sr3_analyze('N-1_IPI_smt', 's01');                %single subject analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, only sequences with same middle transition (smt)
%
%                       [D] = sr3_analyze('sft_ET');                            %group          analysis of same finger transition different sequences on ET (sft_ET)
%                       [D] = sr3_analyze('sft_ET', 's01');                     %single subject analysis of same finger transition different sequences on ET (sft_ET)
%
%                       [D] = sr3_analyze('sft_RT');                            %group          analysis of same finger transition different sequences on RT (sft_RT)
%                       [D] = sr3_analyze('sft_RT', 's01');                     %single subject analysis of same finger transition different sequences on RT (sft_RT)
%
%                       [D] = sr3_analyze('N-1_IPI_sft');                       %group          analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, same transition(s) in different sequences (sft)
%                       [D] = sr3_analyze('N-1_IPI_sft', 's01');                %single subject analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, same transition(s) in different sequences (sft)
%
%                       [D] = sr3_analyze('repEff_N-1_smt_fls');                %group          analysis of N-1 trial go/no-go on repetition effect, same middle transition in sequences with first and last press swapped (smt_fls)
%                       [D] = sr3_analyze('repEff_N-1_smt_fls', 's01');         %single subject analysis of N-1 trial go/no-go on repetition effect, same middle transition in sequences with first and last press swapped (smt_fls)
%
%                       [D] = sr3_analyze('N-1_IPI_smt_fls');                   %group          analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, same middle transition in sequences with first and last press swapped (smt_fls)
%                       [D] = sr3_analyze('N-1_IPI_smt_fls', 's01');            %single subject analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, same middle transition in sequences with first and last press swapped (smt_fls)
%
%                       [D] = sr3_analyze('repEff_scatter');                    %group          analysis of relationship between ET, RT and repetition effect
%                       [D] = sr3_analyze('repEff_scatter', 's01');             %single subject analysis of relationship between ET, RT and repetition effect
%
% --
% gariani@uwo.ca - 2018.05.18


%% paths
pathToData = '../../../../data/SequenceRepetition/sr3';
pathToAnalyze = '../../../../data/SequenceRepetition/sr3/analyze';
if ~exist(pathToAnalyze, 'dir'); mkdir(pathToAnalyze); end % if it doesn't exist already, create analyze folder

%% globals

% subjects
% incomplete: 's17', 's22', 's27', 's34', 's39'
% high-error: 's01', 's08', 's31'
subj = {
    's01', 's02', 's03', 's04', 's05', 's06', 's07', 's08', 's09', 's10', ...
    's11', 's12', 's13', 's14', 's15', 's16',        's18', 's19', 's20', ...
    's21',        's23', 's24', 's25', 's26',        's28', 's29', 's30', ...
    's31', 's32', 's33',        's35', 's36', 's37', 's38',        's40', ...
    's41', 's42', 's43', 's44', 's45', 's46', 's47', 's48', 's49', 's50'
    };

ns = numel(subj);
subvec = zeros(1,ns);
for i = 1:ns; subvec(1,i) = str2double(subj{i}(2:3)); end

% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
blue = [49,130,189]/255;
lightblue = [158,202,225]/255;
red = [222,45,38]/255;
lightred = [252,146,114]/255;
green = [49,163,84]/255;
lightgreen = [161,217,155]/255;
orange = [253,141,60]/255;
yellow = [254,196,79]/255;
lightyellow = [255,237,160]/255;
purple = [117,107,177]/255;
lightpurple = [188,189,220]/255;
darkgray = [50,50,50]/255;
gray2 = [100,100,100]/255;
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;
white = [255,255,255]/255;

% plot defaults
fs = 20; %default font size for all figures
lw = 4; %4; %default line width for all figures
ms = 10; %12; %default marker size for all figures
rs = 0; %0:5; %[0,1,2,3,4,5,6]; %default repetition number subset (full range = 0:10; 0 = full set)
maxRep = 3; % default ceiling level for repetition number (0 = no ceiling)

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,darkgray,gray,gray2,lightgray,green,lightgreen,black,silver,white,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
isrepsty = style.custom({lightgray, darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
lightgraysty = style.custom({lightgray}, 'markertype','none', 'linewidth',lw, 'errorwidth',lw, 'errorcap',0, 'linestyle','--');
darkgraysty = style.custom({darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
blacksty = style.custom({black}, 'markertype','none', 'linewidth',lw, 'linestyle','--','errorbars','plusminus', 'errorwidth',lw, 'errorcap',0);%, 'linestyle','none');
graysty = style.custom({lightgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
nm1sty = style.custom({cbs_pink, cbs_blue}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
sffsty = style.custom({lightgray, gray, darkgray}, 'markersize',ms, 'linewidth',lw);
boxplotsty = style.custom({darkgray}, 'markertype','none', 'linewidth',lw);
isrepstybox = style.custom({gray, black}, 'markersize',lw, 'linewidth',lw);

% legends
isrepleg = {'Switch', 'Repetition'};
nm1leg = {'N-1 nogo', 'N-1 go'};

%% types of analysis
switch (what)
    case 'all_subj' % pre-analysis: run the subject routine for all_subj
        if nargin>1; subj = varargin{1}; end
        for s = 1:numel(subj)
            sr3_subj(subj{s}, 0); % run sr3_subj.m routine (without plot)
        end
        
    case 'all_data' % pre-analysis: create .mat file containing data from all subjects
        all_data = [];
        for s = 1:ns
            fprintf('\n%s\n\n', subj{s});
            D = load(fullfile(pathToData, sprintf('sr3_%s.mat', subj{s}))); % load data structure for each subject
            %-------------------------------------------------------------------------------------------------------------------------------------
            % add IPI info
            D.IPI = diff([D.pressTime1, D.pressTime2, D.pressTime3, D.pressTime4], 1, 2);
            D.IPI_1 = D.IPI(:, 1); D.IPI_2 = D.IPI(:, 2); D.IPI_3 = D.IPI(:, 3);
            
            %-------------------------------------------------------------------------------------------------------------------------------------
            % add N-1 info
            D.nm1 = zeros(numel(D.repType), 1); % n-1 exeType flag: 2=go, 1=nogo, 0=current nogo
            D.nm1(D.repType==111) = 2;        % seq rep,    n-1 go,       go
            D.nm1(D.repType==011) = 2;        % seq switch, n-1 go,       go
            D.nm1(D.repType==101) = 1;        % seq rep,    n-1 nogo,     go
            D.nm1(D.repType==001) = 1;        % seq switch, n-1 nogo,     go
            
            B = [];
            for b = 1:max(D.BN)
                T = getrow(D, D.BN==b);
                %-------------------------------------------------------------------------------------------------------------------------------------
                % detect same first finger (sff) in different sequences
                T.sff = zeros(numel(T.TN), 1);
                for t = 2:numel(T.sff)
                    if (T.isRep(t) == 1)
                        % repetition
                        T.sff(t) = 2;
                    elseif (T.isRep(t) == 0) && (T.press1(t, 1) == T.press1(t-1, 1))
                        % switch, same first finger (sff)
                        T.sff(t) = 1;
                    elseif (T.isRep(t) == 0)
                        % switch
                        T.sff(t) = 0;
                    else % not categorized
                        error('This transition does not fall into any of the defined categories!');
                    end
                end
                %-------------------------------------------------------------------------------------------------------------------------------------
                % detect same middle transition (smt) in different sequences
                T.smt = zeros(numel(T.TN), 1);
                for t = 2:numel(T.smt)
                    if (T.isRep(t) == 1)
                        % repetition
                        T.smt(t) = 2;
                    elseif (T.isRep(t) == 0) && (T.press2(t, 1) == T.press2(t-1, 1)) && (T.press3(t, 1) == T.press3(t-1, 1))
                        % switch, same middle transition (smt)
                        T.smt(t) = 1;
                    elseif (T.isRep(t) == 0)
                        % switch
                        T.smt(t) = 0;
                    else % not categorized
                        error('This transition does not fall into any of the defined categories!');
                    end
                end
                %-------------------------------------------------------------------------------------------------------------------------------------
                % detect same middle transition (smt_fls) in different sequences
                % with the first and last numbers swapped
                T.smt_fls = zeros(numel(T.TN), 1);
                for t = 2:numel(T.smt_fls)
                    if (T.isRep(t) == 1)
                        % repetition
                        T.smt_fls(t) = 2;
                    elseif (T.isRep(t) == 0) && (T.press2(t, 1) == T.press2(t-1, 1))     &&     (T.press3(t, 1) == T.press3(t-1, 1))    &&     (T.press1(t, 1) == T.press4(t-1, 1))     &&     (T.press4(t, 1) == T.press1(t-1, 1))
                        % switch, same middle transition with the first and last numbers swapped (smt_fls)
                        T.smt_fls(t) = 1;
                    elseif (T.isRep(t) == 0)
                        % switch
                        T.smt_fls(t) = 0;
                    else
                        % not categorized
                        error('This transition does not fall into any of the defined categories!');
                    end
                end
                %-------------------------------------------------------------------------------------------------------------------------------------
                % detect same last finger and first of next seq (slf)
                T.slf = zeros(numel(T.TN), 1);
                for t = 2:numel(T.slf)
                    if (T.isRep(t) == 1)
                        % repetition
                        T.slf(t) = 2;
                    elseif (T.press1(t, 1) == T.press4(t-1, 1))
                        % Switch slf
                        T.slf(t) = 1;
                    else
                        % Switch not slf
                        T.slf(t) = 0;
                    end
                end
                %-------------------------------------------------------------------------------------------------------------------------------------
                % detect same finger transition (sft) in different sequences
                T.sft = zeros(numel(T.TN), 1);
                for t = 2:numel(T.sft)
                    % build transition matrix for consecutive sequences
                    seq1 = num2str(T.cuePress(t));
                    tm1 = zeros(5);
                    for i = 2:(numel(seq1))
                        tm1(str2double(seq1(i)), str2double(seq1(i-1))) = 1;
                    end
                    seq2 = num2str(T.cuePress(t-1));
                    tm2 = zeros(5);
                    for i = 2:(numel(seq2))
                        tm2(str2double(seq2(i)), str2double(seq2(i-1))) = 1;
                    end
                    %
                    if (T.isRep(t) == 1)
                        % repetition
                        T.sft(t) = 123;
                    elseif (T.isRep(t) == 0) && (T.press2(t, 1) == T.press2(t-1, 1)) && (T.press3(t, 1) == T.press3(t-1, 1)) && (T.press4(t, 1) == T.press4(t-1, 1))
                        % switch, same scond and third transitions
                        T.sft(t) = 23;
                    elseif (T.isRep(t) == 0) && (T.press1(t, 1) == T.press1(t-1, 1)) && (T.press2(t, 1) == T.press2(t-1, 1)) && (T.press3(t, 1) == T.press3(t-1, 1))
                        % switch, same first and second transitions
                        T.sft(t) = 12;
                    elseif (T.isRep(t) == 0) && (T.press3(t, 1) == T.press3(t-1, 1)) && (T.press4(t, 1) == T.press4(t-1, 1))
                        % switch, same third transition
                        T.sft(t) = 3;
                    elseif (T.isRep(t) == 0) && (T.press2(t, 1) == T.press2(t-1, 1)) && (T.press3(t, 1) == T.press3(t-1, 1))
                        % switch, same second transition
                        T.sft(t) = 2;
                    elseif (T.isRep(t) == 0) && (T.press1(t, 1) == T.press1(t-1, 1)) && (T.press2(t, 1) == T.press2(t-1, 1))
                        % switch, same first transition
                        T.sft(t) = 1;
                    elseif (T.isRep(t) == 0) && ~any(any(tm1==1 & tm1==tm2))
                        % switch, all transitions are different
                        T.sft(t) = 0;
                    elseif (T.isRep(t) == 0)
                        % switch, some transitions may be shared
                        T.sft(t) = -1;
                    else
                        % not categorized
                        error('This transition does not fall into any of the defined categories!');
                    end
                end
                B = addstruct(B, T);
            end
            all_data = addstruct(all_data, B); % append data structures from each subject
        end
        save( fullfile( pathToAnalyze, 'sr3_all_data.mat'), '-struct', 'all_data'); % save all_data.mat file
        varargout = {all_data};
        
    case 'repEff_ET' % analysis of repetition effect on ET
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Repetition effect - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Repetition effect - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %%% figure 1 %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep eff
        T = tapply(D, {'SN', 'isRep'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.isRep, T.SN, T.normET, 'length');
        
        subplot(2,2,1); title('Repetition effect'); hold on;
        plt.box(T.isRep, T.normET, 'style',boxplotsty);
        
        T.isRep(T.isRep == 1) = 2; T.isRep(T.isRep == 0) = 1;
        hold on;
        plt.line(T.isRep, T.normET, 'plotfcn','mean', 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','off');
        hold on;
        plt.line(T.isRep, T.normET, 'plotfcn','mean', 'errorbars','plusminus', 'style',blacksty, 'leg','off');
        xticklabels({'Switch', 'Repetition'}); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square;
        ylim([640 780]);
        
        % stats
        [T.t1,T.p1]=ttest(T.ET(T.isRep==1), T.ET(T.isRep==2), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep num
        if maxRep > 0; D.repNum(D.repNum >= maxRep) = maxRep; end % put a ceiling to nReps
        T = tapply(D, {'SN', 'repNum'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);% & ismember(D.repType,[111,11]));
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.repNum, T.SN, T.normET, 'length');
        
        % stats
        %ttest(T.ET(T.repNum==2), T.ET(T.repNum==4), 2, 'paired');
        %ttest(T.ET(T.repNum==2), T.ET(T.repNum==3), 2, 'paired');
        %ttest(T.ET(T.repNum==3), T.ET(T.repNum==4), 2, 'paired');
        
        subplot(2,2,2); title('Repetition number');
        %[~,~] = plt.line(T.repNum, T.normET, 'split',[T.SN], 'errorbars','shade', 'style',lightgraysty, 'leg','off',  'leglocation','northeast');
        %hold on;
        [~,~] = plt.line(T.repNum, T.normET, 'plotfcn','mean', 'errorbars','shade', 'style',darkgraysty, 'leg',{'Group median'},  'leglocation','northeast');
        xlabel('Repetition number'); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square;
        ylim([640 780]);
        %plt.match('y')
        labels = cell(1, max(D.repNum)+1); labels{1, 1} = 'Switch';
        for l = 1 : max(D.repNum); labels{1, l+1} = num2str(l); end
        labels{1, end} = sprintf('%s+', labels{1, end}); xticklabels(labels);
        
        % stats
        [T.t2,T.p2]=ttest(T.ET(T.repNum==0), T.ET(T.repNum==1), 2, 'paired');
        [T.t3,T.p3]=ttest(T.ET(T.repNum==1), T.ET(T.repNum==2), 2, 'paired');
        [T.t4,T.p4]=ttest(T.ET(T.repNum==2), T.ET(T.repNum==3), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        nq = 11; % number of quantiles
        [D.I,D.M,~] = split_data(D.ET, 'split',[D.SN D.isRep], 'numquant',nq, 'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        T = tapply(D, {'SN', 'I', 'isRep'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        subplot(2,2,3); title(''); hold on;
        plt.line(T.I, T.ET, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xticklabels(linspace(0,100,nq));
        xlabel('ET percentile (%)'); ylabel('ET (ms)');  set(gca,'fontsize',fs); axis square;
        xlim([0.5 nq+0.5]);
        ylim([420 1280]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        [D.I,D.M,~] = split_data(D.ET, 'split',[D.SN D.isRep], 'numquant',nq, 'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        T = tapply(D, {'SN', 'I'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.ET, 'nanmedian', 'name', 'ETs', 'subset',D.isRep==0}, ...
            {D.ET, 'nanmedian', 'name', 'ETr', 'subset',D.isRep==1}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ETs', 'ETr', 'ET'}, 'sub');
        
        subplot(2,2,4); title(''); hold on;
        plt.line(T.I, (nanplus(T.ETs,-T.ETr)./T.ET)*100, 'style',darkgraysty);
        xticklabels(linspace(0,100,nq));
        xlabel('ET percentile (%)'); ylabel('Repetition difference (% of ET)');  set(gca,'fontsize',fs); axis square;
        xlim([0.5 nq+0.5]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--');
        
        % stats
        ttest((T.ETs(T.I==1)-T.ETr(T.I==1)), 0, 2, 'onesample');
        
        % open figure
        if nargin>1; figure('Name',sprintf('Repetition effect - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Repetition effect - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        
        %%% figure 2 -------------------------------------------------------------------------------------------------------------------------------------
        % Train
        T = tapply(D, {'SN', 'BN', 'isRep'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.isRep, T.SN, T.normET, 'length');
        
        subplot(2,2,1); title('Training');
        [~,~] = plt.line(T.BN, T.normET, 'split',[T.isRep], 'errorbars','shade', 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        xlabel('Block number'); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square; ylim([480 1020]); %ylim([600 1000]);
        xlim([0 25]);
        
        % stats
        ttest(T.ET(T.BN==1 & T.isRep==0), T.ET(T.BN==1 & T.isRep==1), 2, 'paired');
        ttest(T.ET(T.BN==12 & T.isRep==0), T.ET(T.BN==12 & T.isRep==1), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Diff
        T = tapply(D, {'SN', 'BN'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.ET, 'nanmedian', 'name', 'ETswc', 'subset', D.isRep==0}, ...
            {D.ET, 'nanmedian', 'name', 'ETrep', 'subset', D.isRep==1}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET', 'ETswc', 'ETrep'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.isRep, T.SN, T.normET, 'length');
        
        subplot(2,2,2); title('Difference');
        %[~,~] = plt.line(T.BN, (T.normETswc-T.normETrep), 'split',[T.SN], 'errorbars','shade', 'style',lightgraysty, 'leg','off', 'leglocation','southeast');
        %hold on;
        
        %[~,~] = plt.line(T.BN, ((T.ETswc-T.ETrep)./T.ET)*100, 'split',[], 'errorbars','shade', 'style',graysty, 'leg','skip');
        [~,~] = plt.line(T.BN, (nanplus(T.ETswc,-T.ETrep)./T.ET)*100, 'split',[], 'errorbars','shade', 'style',graysty, 'leg','skip');
        hold on;
        [~,~] = plt.scatter(T.BN, (nanplus(T.ETswc,-T.ETrep)./T.ET)*100, 'style',blacksty, 'leg','skip');
        xlabel('Block number'); ylabel('Repetition difference (% of ET)'); set(gca,'fontsize',fs); axis square; ylim([-1.5 7.5]); %ylim([-15 65]); %ylim([-200 200]);
        xlim([0 25]);
        drawline(0,'dir','horz','linestyle','--','color',black);
        
        % stats
        delta = ((T.ETswc-T.ETrep)./T.ET)*100;
        ttest(delta(T.BN==1), delta(T.BN==12), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % create summary tables for ACC
        T = tapply(D,{'SN', 'isRep'},...
            {(1-D.isError)*100,'nanmean', 'name','ACC'}, ...
            'subset',D.timingError==0 & D.exeType==1);
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.day, T.BN, T.ACC, 'length');
        
        subplot(2,2,3); hold on;
        plt.box(T.isRep, T.normACC, 'style',boxplotsty);
        T.isRep(T.isRep == 1) = 2; T.isRep(T.isRep == 0) = 1;
        hold on;
        plt.line(T.isRep, T.normACC, 'plotfcn','mean', 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','off');
        hold on;
        plt.line(T.isRep, T.normACC, 'plotfcn','mean', 'errorbars','plusminus', 'style',blacksty, 'leg','off');
        xticklabels({'Switch', 'Repetition'}); ylabel('Accuracy (%)'); set(gca,'fontsize',fs); axis square;
        ylim([73 97]);
        
        % stats
        ttest(T.ACC(T.isRep==1), T.ACC(T.isRep==2), 2, 'paired');
        
        % how many false starts (timing errors)? --> 3.48% in total
        %fs = (sum(D.timingError)/numel(D.timingError))*100;
        
        T = tapply(D,{'SN', 'isRep'},...
            {(D.timingError)*100,'nanmean', 'name','tACC'});
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'tACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.day, T.BN, T.ACC, 'length');
        
        subplot(2,2,3); hold on;
        plt.box(T.isRep, T.normtACC, 'style',boxplotsty);
        T.isRep(T.isRep == 1) = 2; T.isRep(T.isRep == 0) = 1;
        hold on;
        plt.line(T.isRep, T.normtACC, 'plotfcn','mean', 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','off');
        hold on;
        plt.line(T.isRep, T.normtACC, 'plotfcn','mean', 'errorbars','plusminus', 'style',blacksty, 'leg','off');
        xticklabels({'Switch', 'Repetition'}); ylabel('False start (%)'); set(gca,'fontsize',fs); axis square;
        %ylim([73 97]);
        
        % stats
        ttest(T.tACC(T.isRep==1), T.tACC(T.isRep==2), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep num
        if maxRep > 0; D.repNum(D.repNum >= maxRep) = maxRep; end % put a ceiling to nReps
        T = tapply(D, {'SN', 'repNum'}, ...
            {(1-D.isError)*100,'nanmean', 'name','ACC'}, ...
            'subset',D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.repNum, T.SN, T.ACC, 'length');
        
        subplot(2,2,4);
        plt.line(T.repNum, T.normACC, 'errorbars','shade', 'style',darkgraysty, 'leglocation','northeast');
        xlabel('Repetition number'); ylabel('ACC (%)'); set(gca,'fontsize',fs); axis square;
        ylim([73 97]);
        labels = cell(1, max(D.repNum)+1); labels{1, 1} = 'Switch';
        for l = 1 : max(D.repNum); labels{1, l+1} = num2str(l); end
        labels{1, end} = sprintf('%s+', labels{1, end}); xticklabels(labels);
        
        % stats
        ttest(T.ACC(T.repNum==0), T.ACC(T.repNum==1), 2, 'paired');
        ttest(T.ACC(T.repNum==1), T.ACC(T.repNum==2), 2, 'paired');
        ttest(T.ACC(T.repNum==2), T.ACC(T.repNum==3), 2, 'paired');
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'N-1_ET' % analysis of N-1 trial (go/no-go) on repetition effect for ET
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('N-1 on repetition effect - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('N-1 on repetition effect - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %         %-------------------------------------------------------------------------------------------------------------------------------------
        %         % n-1 no-go
        %         T = tapply(D, {'SN', 'BN', 'repType'}, ...
        %             {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.nm1==1}, ... %{D.ET, 'nanmedian', 'name', 'ET', 'subset', D.repType==1 | D.repType==11}, ...
        %             'subset', D.isError==0 & D.timingError==0);
        %
        %         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'ET'}, 'sub');
        %
        %         % make sure that you have one value per subject for each condition
        %         % pivottable(T.repType, T.SN, T.normET, 'length');
        %
        %         subplot(2,2,1); title('N-1 no-go trials');
        %         [~,~] = plt.line(T.BN, T.normET, 'split',T.repType, 'errorbars','shade', 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        %         xlabel('Block number'); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square;
        %
        %         %-------------------------------------------------------------------------------------------------------------------------------------
        %         % n-1 go
        %         T = tapply(D, {'SN', 'BN', 'repType'}, ...
        %             {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.nm1==2}, ...
        %             'subset', D.isError==0 & D.timingError==0);% & D.repNum<=1);
        %
        %         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'ET'}, 'sub');
        %
        %         % make sure that you have one value per subject for each condition
        %         % pivottable(T.repType, T.SN, T.normET, 'length');
        %
        %         subplot(2,2,2); title('N-1 go trials');
        %         [~,~] = plt.line(T.BN, T.normET, 'split',T.repType, 'errorbars','shade', 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        %         xlabel('Block number'); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square; plt.match('y');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % box plot
        T = tapply(D, {'SN', 'isRep', 'nm1'}, ...
            {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0 & D.repNum<2);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.repType, T.SN, T.normET, 'length');
        
        subplot(2,2,1); title('ANOVA'); hold on;
        plt.box(T.nm1, T.normET, 'split',T.isRep, 'style',isrepstybox, 'leg',isrepleg, 'leglocation','northeast');
        hold on;
        plt.line([T.nm1 T.isRep], T.normET, 'plotfcn','mean', 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','skip');
        hold on;
        plt.line([T.nm1 T.isRep], T.normET, 'plotfcn','mean', 'errorbars','plusminus', 'style',blacksty);
        xlabel('Previous trial'); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square; %ylim([635 835]);
        xt = xticks; xticks([sum(xt(1:2))/2 sum(xt(3:4))/2]); xticklabels({'No-go', 'Go'});
        
        % stats
        %ttest(T.ET(T.nm1==1 & T.isRep==0), T.ET(T.nm1==1 & T.isRep==1), 2, 'paired');
        %ttest(T.ET(T.nm1==2 & T.isRep==0), T.ET(T.nm1==2 & T.isRep==1), 2, 'paired');
        %T.ANOVA = anovaMixed(T.ET, T.SN,'within', [T.nm1, T.isRep], {'n-1','isRep'});
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        nq = 11; % number of quantiles
        [D.I,D.M,~] = split_data(D.ET, 'split',[D.SN D.isRep D.nm1], 'numquant',nq, 'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        T = tapply(D, {'SN', 'I', 'isRep', 'nm1'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        subplot(2,2,2); title('Previous trial No-go (left), or Go (right)'); hold on;
        plt.line([T.nm1 T.I], T.normET, 'split',[T.isRep], 'style',isrepsty, 'leg',isrepleg, 'leglocation','northwest');
        xticklabels(repmat(linspace(0,100,nq),1,2));
        xlabel('ET percentile (%)'); ylabel('ET (ms)');  set(gca,'fontsize',fs); %axis square;
        ylim([420 1280]);
        xt = xticks; drawline([xt(6),xt(17)], 'dir','vert', 'linestyle','--');
        xlim([xt(1)-0.5 xt(end)+0.5]);
        
        ttest( T.ET(T.nm1==1 & T.I==11 & T.isRep==0) - T.ET(T.nm1==1 & T.I==11 & T.isRep==1), 0, 2, 'onesample');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % box plot
        T = tapply(D, {'SN', 'isRep', 'nm1'}, ...
            {(1-D.isError)*100,'nanmean', 'name','ACC'}, ...
            'subset',D.timingError==0 & D.exeType==1 & D.repNum<2);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ACC'}, 'sub');
        
        subplot(2,2,3);
        plt.box(T.nm1, T.normACC, 'split',T.isRep, 'style',isrepstybox, 'leg',isrepleg, 'leglocation','northeast');
        hold on;
        plt.line([T.nm1 T.isRep], T.normACC, 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','skip');
        hold on;
        plt.line([T.nm1 T.isRep], T.normACC, 'errorbars','plusminus', 'style',blacksty);
        xlabel('Previous trial'); ylabel('Accuracy (%)'); set(gca,'fontsize',fs); axis square; ylim([69 101]);
        xt = xticks; xticks([sum(xt(1:2))/2 sum(xt(3:4))/2]); xticklabels({'No-go', 'Go'});
        
        % stats
        %ttest(T.ACC(T.nm1==1 & T.isRep==0), T.ACC(T.nm1==1 & T.isRep==1), 2, 'paired');
        %ttest(T.ACC(T.nm1==2 & T.isRep==0), T.ACC(T.nm1==2 & T.isRep==1), 2, 'paired');
        %T.ANOVA = anovaMixed(T.ACC, T.SN,'within', [T.nm1, T.isRep], {'n-1','isRep'});
        
        T = tapply(D, {'SN', 'isRep', 'nm1'}, ...
            {(D.timingError)*100,'nanmean', 'name','tACC'}, ...
            'subset', D.exeType==1 & D.repNum<2);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'tACC'}, 'sub');
        
        subplot(2,2,3);
        plt.box(T.nm1, T.normtACC, 'split',T.isRep, 'style',isrepstybox, 'leg',isrepleg, 'leglocation','northeast');
        hold on;
        plt.line([T.nm1 T.isRep], T.normtACC, 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','skip');
        hold on;
        plt.line([T.nm1 T.isRep], T.normtACC, 'errorbars','plusminus', 'style',blacksty);
        xlabel('Previous trial'); ylabel('False start (%)'); set(gca,'fontsize',fs); axis square; %ylim([69 101]);
        xt = xticks; xticks([sum(xt(1:2))/2 sum(xt(3:4))/2]); xticklabels({'No-go', 'Go'});
        
        % stats
        ttest(T.tACC(T.nm1==1 & T.isRep==0), T.tACC(T.nm1==1 & T.isRep==1), 2, 'paired');
        ttest(T.tACC(T.nm1==2 & T.isRep==0), T.tACC(T.nm1==2 & T.isRep==1), 2, 'paired');
        T.ANOVA = anovaMixed(T.tACC, T.SN,'within', [T.nm1, T.isRep], {'n-1','isRep'});
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        [D.I,D.M,~] = split_data(D.ET, 'split',[D.SN D.isRep D.nm1], 'numquant',nq, 'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        T = tapply(D, {'SN', 'I', 'nm1'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.ET, 'nanmedian', 'name', 'ETs', 'subset',D.isRep==0}, ...
            {D.ET, 'nanmedian', 'name', 'ETr', 'subset',D.isRep==1}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ETs', 'ETr', 'ET'}, 'sub');
        
        subplot(2,2,4); title(''); hold on;
        [x,y]=plt.line(T.I, ((T.ETs-T.ETr)./T.ET)*100, 'split',T.nm1, 'style',isrepsty, 'leg',nm1leg, 'leglocation','northwest');
        xticklabels(linspace(0,100,nq));
        xlabel('ET percentile (%)'); ylabel('Repetition difference (% of ET)');  set(gca,'fontsize',fs); axis square;
        xt = xticks; xlim([xt(1)-0.5 xt(end)+0.5]); ylim([-5 9]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--'); drawline(0, 'dir','horz', 'linestyle',':');
        
        ttest(((T.ETs(T.nm1==2 & T.I==1)-T.ETr(T.nm1==2 & T.I==1))./T.ET(T.nm1==2 & T.I==1))*100, ((T.ETs(T.nm1==2 & T.I==11)-T.ETr(T.nm1==2 & T.I==11))./T.ET(T.nm1==2 & T.I==11))*100, 2, 'paired')
        
        % stats
        T = tapply(D, {'SN', 'nm1'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.ET, 'nanmedian', 'name', 'ETs', 'subset',D.isRep==0}, ...
            {D.ET, 'nanmedian', 'name', 'ETr', 'subset',D.isRep==1}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        %[t1,p1]=ttest(((T.ETs(T.nm1==1)-T.ETr(T.nm1==1))./T.ET(T.nm1==1))*100, 0, 2, 'onesample');
        %[t2,p2]=ttest(((T.ETs(T.nm1==2)-T.ETr(T.nm1==2))./T.ET(T.nm1==2))*100, 0, 2, 'onesample');
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'repEff_RT' % analysis of repetition effect on execution time (RT)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('RT analysis - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('RT analysis - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep eff
        T = tapply(D, {'SN', 'isRep'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.isRep, T.SN, T.normRT, 'length');
        
        subplot(2,2,1); title('Repetition effect'); hold on;
        plt.box(T.isRep, T.normRT, 'style',boxplotsty);
        
        T.isRep(T.isRep == 1) = 2; T.isRep(T.isRep == 0) = 1;
        hold on;
        plt.line(T.isRep, T.normRT, 'plotfcn','mean', 'split',[T.SN], 'errorbars','', 'style',lightgraysty, 'leg','off');
        hold on;
        plt.line(T.isRep, T.normRT, 'plotfcn','mean', 'errorbars','plusminus', 'style',blacksty, 'leg','off');
        xticklabels({'Switch', 'Repetition'}); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square; ylim([410 490]);
        
        % stats
        [T.t1,T.p1]=ttest(T.RT(T.isRep==1), T.RT(T.isRep==2), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep num
        if maxRep > 0; D.repNum(D.repNum >= maxRep) = maxRep; end % put a ceiling to nReps
        T = tapply(D, {'SN', 'repNum'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.repNum, T.SN, T.normRT, 'length');
        
        % stats
        [T.t2,T.p2]=ttest(T.RT(T.repNum==0), T.RT(T.repNum==1), 2, 'paired');
        [T.t3,T.p3]=ttest(T.RT(T.repNum==1), T.RT(T.repNum==2), 2, 'paired');
        [T.t4,T.p4]=ttest(T.RT(T.repNum==2), T.RT(T.repNum==3), 2, 'paired');
        
        
        subplot(2,2,2); title('Repetition number');
        %[~,~] = plt.line(T.repNum, T.normRT, 'split',[T.SN], 'errorbars','shade', 'style',lightgraysty, 'leg','off',  'leglocation','northeast');
        %hold on;
        [~,~] = plt.line(T.repNum, T.normRT, 'plotfcn','mean', 'errorbars','shade', 'style',darkgraysty, 'leg',{'Group median'},  'leglocation','northeast');
        xlabel('Repetition number'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square; ylim([410 490]);
        labels = cell(1, max(D.repNum)+1); labels{1, 1} = 'Switch';
        for l = 1 : max(D.repNum); labels{1, l+1} = num2str(l); end
        labels{1, end} = sprintf('%s+', labels{1, end}); xticklabels(labels);
        %plt.match('y');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        nq = 11; % number of quantiles
        [D.I,D.M,~] = split_data(D.RT, 'split',[D.SN D.isRep], 'numquant',nq, 'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        T = tapply(D, {'SN', 'I', 'isRep'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        subplot(2,2,3); title(''); hold on;
        plt.line(T.I, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xticklabels(linspace(0,100,nq));
        xlabel('RT percentile (%)'); ylabel('RT (ms)');  set(gca,'fontsize',fs); axis square;
        xlim([0.5 nq+0.5]);
        ylim([300 800]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        [D.I,D.M,~] = split_data(D.RT, 'split',[D.SN D.isRep], 'numquant',nq, 'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        T = tapply(D, {'SN', 'I'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            {D.RT, 'nanmedian', 'name', 'RTs', 'subset',D.isRep==0}, ...
            {D.RT, 'nanmedian', 'name', 'RTr', 'subset',D.isRep==1}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RTs', 'RTr', 'RT'}, 'sub');
        
        subplot(2,2,4); title(''); hold on;
        plt.line(T.I, (nanplus(T.normRTs,-T.normRTr)./T.normRT)*100, 'style',darkgraysty);
        xticklabels(linspace(0,100,nq));
        xlabel('RT percentile (%)'); ylabel('Repetition difference (% of RT)');  set(gca,'fontsize',fs); axis square;
        xlim([0.5 nq+0.5]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--');
        
        %         %-------------------------------------------------------------------------------------------------------------------------------------
        %         % Train
        %         T = tapply(D, {'SN', 'BN', 'isRep'}, ...
        %             {D.RT, 'nanmedian', 'name', 'RT'}, ...
        %             'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        %
        %         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'RT'}, 'sub');
        %
        %         % make sure that you have one value per subject for each condition
        %         % pivottable(T.isRep, T.SN, T.normRT, 'length');
        %
        %         subplot(2,2,3); title('Training');
        %         [~,~] = plt.line(T.BN, T.normRT, 'split',[T.isRep], 'errorbars','shade', 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        %         xlabel('Block number'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square; ylim([400 600]);
        %
        %         %-------------------------------------------------------------------------------------------------------------------------------------
        %         % Diff
        %         T = tapply(D, {'SN', 'BN'}, ...
        %             {D.RT, 'nanmedian', 'name', 'RTswc', 'subset', D.isRep==0}, ...
        %             {D.RT, 'nanmedian', 'name', 'RTrep', 'subset', D.isRep==1}, ...
        %             'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        %
        %         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'RTswc', 'RTrep'}, 'sub');
        %
        %         % make sure that you have one value per subject for each condition
        %         % pivottable(T.isRep, T.SN, T.normRT, 'length');
        %
        %         subplot(2,2,4); title('Difference');
        %         %[~,~] = plt.line(T.BN, (T.normRTswc-T.normRTrep), 'split',[T.SN], 'errorbars','shade', 'style',lightgraysty, 'leg','off', 'leglocation','southeast');
        %         %hold on;
        %         [~,~] = plt.line(T.BN, (T.normRTswc-T.normRTrep), 'split',[], 'errorbars','shade', 'style',blacksty, 'leg',{'Group median'},  'leglocation','southeast');
        %         xlabel('Block number'); ylabel('Switch - Repetition RT (ms)'); set(gca,'fontsize',fs); axis square; ylim([-200 200]);
        %         drawline(0,'dir','horz','linestyle','--','color',black);
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'N-1_RT' % analysis of N-1 trial (go/no-go) on repetition effect for RT
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        %         % open figure
        if nargin>1; figure('Name',sprintf('N-1 on repetition effect for RT - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('N-1 on repetition effect for RT - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %         %-------------------------------------------------------------------------------------------------------------------------------------
        %         % N-1 no-go trials (Swc trials)
        %         T = tapply(D, {'SN', 'BN', 'repType'}, ...
        %             {D.RT, 'nanmedian', 'name', 'RT', 'subset', D.nm1==1}, ... , D.repType==1 | D.repType==11}, ...
        %             'subset', D.isError==0 & D.timingError==0);
        %
        %         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'RT'}, 'sub');
        %
        %         % make sure that you have one value per subject for each condition
        %         % pivottable(T.repType, T.SN, T.normRT, 'length');
        %
        %         subplot(2,2,1); title('N-1 no-go trials');
        %         [~,~] = plt.line(T.BN, T.normRT, 'split',[T.repType], 'errorbars','shade', 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        %         xlabel('Block number'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square; plt.match('y');
        %
        %         %-------------------------------------------------------------------------------------------------------------------------------------
        %         % N-1 go trials (Rep trials)
        %         T = tapply(D, {'SN', 'BN', 'repType'}, ...
        %             {D.RT, 'nanmedian', 'name', 'RT', 'subset', D.nm1==2}, ... , D.repType==101 | D.repType==111}, ...
        %             'subset', D.isError==0 & D.timingError==0);% & D.repNum<=1);
        %
        %         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'RT'}, 'sub');
        %
        %         % make sure that you have one value per subject for each condition
        %         % pivottable(T.repType, T.SN, T.normRT, 'length');
        %
        %         subplot(2,2,2); title('N-1 go trials');
        %         [~,~] = plt.line(T.BN, T.normRT, 'split',[T.repType], 'errorbars','shade', 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        %         xlabel('Block number'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square;
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % box plot
        T = tapply(D, {'SN', 'isRep', 'nm1'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.repType, T.SN, T.normRT, 'length');
        
        %         subplot(2,2,3); title('ANOVA');
        %         [x, ~, ~] = plt.bar(T.nm1, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        %         xlabel('N-1 execution type'); ylabel('Sequence execution time (ms)'); set(gca,'fontsize',fs); axis square;
        %         ylim([700 900]); xticklabels({'no-go', 'go'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2]);
        
        subplot(2,2,1); title('ANOVA'); hold on;
        %[x, ~, ~] = plt.bar(T.nm1, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        %xlabel(''); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square;
        %ylim([400 600]); xticklabels(nm1leg); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2]);
        plt.box(T.nm1, T.normRT, 'split',T.isRep, 'style',isrepstybox, 'leg',isrepleg, 'leglocation','northeast');
        hold on;
        plt.line([T.nm1 T.isRep], T.normRT, 'plotfcn','mean', 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','skip');
        hold on;
        plt.line([T.nm1 T.isRep], T.normRT, 'plotfcn','mean', 'errorbars','plusminus', 'style',blacksty);
        xlabel('Previous trial'); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square; %ylim([340 580]);
        xt = xticks; xticks([sum(xt(1:2))/2 sum(xt(3:4))/2]); xticklabels({'No-go', 'Go'});
        
        % stats
        ttest(T.RT(T.nm1==1 & T.isRep==0), T.RT(T.nm1==1 & T.isRep==1), 2, 'paired');
        ttest(T.RT(T.nm1==2 & T.isRep==0), T.RT(T.nm1==2 & T.isRep==1), 2, 'paired');
        T.ANOVA = anovaMixed(T.RT, T.SN,'within', [T.nm1, T.isRep], {'n-1','isRep'});
        
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        nq = 11; % number of quantiles
        [D.I,D.M,~] = split_data(D.RT, 'split',[D.SN D.isRep D.nm1], 'numquant',nq, 'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        T = tapply(D, {'SN', 'I', 'isRep', 'nm1'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        subplot(2,2,2); title('Previous trial No-go (left), or Go (right)'); hold on;
        plt.line([T.nm1 T.I], T.normRT, 'split',[T.isRep], 'style',isrepsty, 'leg',isrepleg, 'leglocation','northwest');
        xticklabels(repmat(linspace(0,100,nq),1,2)); xlabel('RT percentile (%)'); ylabel('RT (ms)');  set(gca,'fontsize',fs); %axis square;
        ylim([300 800]);
        xt = xticks; drawline([xt(6),xt(17)], 'dir','vert', 'linestyle','--');
        xlim([xt(1)-0.5 xt(end)+0.5]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        %nq = 11; % number of quantiles
        [D.I,D.M,~] = split_data(D.RT, 'split',[D.SN D.isRep D.nm1], 'numquant',nq, 'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        T = tapply(D, {'SN', 'I', 'nm1'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            {D.RT, 'nanmedian', 'name', 'RTs', 'subset',D.isRep==0}, ...
            {D.RT, 'nanmedian', 'name', 'RTr', 'subset',D.isRep==1}, ...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RTs', 'RTr', 'RT'}, 'sub');
        
        subplot(2,2,3:4); title(''); hold on;
        plt.line(T.I, (nanplus(T.normRTs,-T.normRTr)./T.normRT)*100, 'split',T.nm1, 'style',isrepsty, 'leg',nm1leg, 'leglocation','northwest');
        xticklabels(linspace(0,100,nq)); xlabel('RT percentile (%)'); ylabel('Repetition difference (% of RT)');  set(gca,'fontsize',fs); axis square;
        xt = xticks; xlim([xt(1)-0.5 xt(end)+0.5]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--'); drawline(0, 'dir','horz', 'linestyle',':');
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'firstFinger_ET' % analysis of same first finger repetition (different sequence) influence on ET (vs same sequence)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Same first finger on ET - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Same first finger on ET - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % create summary table
        T = tapply(D, {'SN', 'sff'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.sff, T.SN, T.normET, 'length');
        
        subplot(2,2,1:4); title('First finger influence on ET'); hold on;
        [x, y, e] = plt.bar(T.sff, T.normET, 'split',T.sff, 'style',sffsty, 'capwidth',0.1, 'leg','off');
        xticklabels({'Swc', 'Sff', 'Rep'}); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square; ylim([600 800]);
        
        % stats
        [t,p] = ttest(T.ET(T.sff==0), T.ET(T.sff==1), 2, 'paired');
        text(x(1)-.3, max(y)+sum(e), sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',fs);
        [t,p] = ttest(T.ET(T.sff==0), T.ET(T.sff==2), 2, 'paired');
        text(x(2)-.3, max(y)+sum(e)*1.5, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',fs);
        [t,p] = ttest(T.ET(T.sff==1), T.ET(T.sff==2), 2, 'paired');
        text(x(3)-.3, max(y)+sum(e), sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',fs);
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'firstFinger_RT' % analysis of same first finger repetition (different sequence) influence on RT (vs same sequence)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Same first finger on RT - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Same first finger on RT - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % detect same first finger (sff) different sequences
        D.sff = zeros(numel(D.TN), 1);
        for t = 2:numel(D.TN)
            if (D.press1(t, 1) == D.press1(t-1, 1)) && (D.isRep(t) == 0)
                D.sff(t) = 1;
            elseif (D.press1(t, 1) == D.press1(t-1, 1)) && (D.isRep(t) == 1)
                D.sff(t) = 2;
            else
                D.sff(t) = 0;
            end
        end
        
        % create summary table
        T = tapply(D, {'SN', 'sff'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.sff, T.SN, T.normRT, 'length');
        
        subplot(2,2,1:4); title('First finger influence on RT'); hold on;
        [x, y, e] = plt.bar(T.sff, T.normRT, 'split',T.sff, 'style',sffsty, 'capwidth',0.1, 'leg','off');
        xticklabels({'Swc', 'Sff', 'Rep'}); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square; ylim([400 500]);
        
        % stats
        [t,p] = ttest(T.RT(T.sff==0), T.RT(T.sff==1), 2, 'paired');
        text(x(1)-.3, max(y)+sum(e), sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',fs);
        [t,p] = ttest(T.RT(T.sff==0), T.RT(T.sff==2), 2, 'paired');
        text(x(2)-.3, max(y)+sum(e)*1.5, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',fs);
        [t,p] = ttest(T.RT(T.sff==1), T.RT(T.sff==2), 2, 'paired');
        text(x(3)-.3, max(y)+sum(e), sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',fs);
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
    
    case 'slf_RT' % analysis of same last and first finger of next seq (different sequence) influence on RT
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('SLF on RT - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('SLF on RT - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % create summary table
        T = tapply(D, {'SN', 'slf'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.slf, T.SN, T.normRT, 'length');
        
        subplot(2,2,1:4); title('First finger influence on RT'); hold on;
        [x, y, e] = plt.bar(T.slf, T.normRT, 'style',sffsty, 'capwidth',0.1, 'leg','off');
        xticklabels({'not SLF', 'SLF', 'Rep'}); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square; ylim([400 500]);
        
        % stats
        [t,p] = ttest(T.RT(T.slf==0), T.RT(T.slf==1), 2, 'paired');
        text(x(1)-.3, max(y)+sum(e), sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',fs);
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'FAp_RT' % analysis of False Alarm probability influence on RT (no-go trials only)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('False alarm RT analysis - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('False alarm RT analysis - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % flag False Alarms (FAs)
        D.FA = zeros(numel(D.TN), 1);
        D.FA(D.exeType == 0 & D.response1 > 0) = 1;
        
        % create summary table
        T = tapply(D, {'SN', 'isRep'}, ...
            {D.FA, 'nanmedian', 'name', 'FAp'}, ...
            {D.RT, 'nanmedian', 'name', 'mRT', 'subset', D.isError==0 & D.timingError==0 & D.exeType==1});
        %T.FA_bin = round(T.FAp, 4);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'mRT', 'FAp'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.BN, T.SN, T.normFAp, 'length');
        
        subplot(2,2,[1,3]); title('All subjects'); hold on;
        plt.scatter(T.normFAp, T.normmRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xlabel('FA probability'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square; xlim([0 max(T.FAp)]); ylim([200 700]);
        
        % create summary table
        T = tapply(D, {'SN', 'BN', 'isRep'}, ...
            {D.FA, 'nanmedian', 'name', 'FAp'}, ...
            {D.RT, 'nanmedian', 'name', 'mRT', 'subset', D.isError==0 & D.timingError==0 & D.exeType==1});
        %T.FA_bin = round(T.FAp, 4);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'mRT', 'FAp'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.BN, T.SN, T.normFAp, 'length');
        
        subplot(2,2,[2,4]); title('All subjects and blocks'); hold on;
        plt.scatter(T.normFAp, T.normmRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xlabel('FA probability'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square; xlim([0 max(T.FAp)]); ylim([200 700]);
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'ET_RT' % analysis of relationship between ETs and RTs (correlation)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('ET vs RT relationship - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('ET vs RT relationship - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        T = tapply(D, {'SN', 'BN', 'TN', 'isRep'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET', 'RT'}, 'sub');
        
        subplot(2,2,1:4); title('ET vs RT'); hold on;
        plt.scatter(T.normET, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xlabel('ET (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis image;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'first2IPIs_RT' % analysis of relationship between first 2 IPIs and RTs (correlation)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('IPI vs RT relationship - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('IPI vs RT relationship - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % add IPI info
        D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4],1,2);
        D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % create summary table for IPI profile
        T=tapply(D, {'SN', 'isRep'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        for i = 1:size(D.IPI, 2)
            T.IPI(:,i) = eval( sprintf('T.IPI_%d', i));
            T = rmfield(T,sprintf('IPI_%d', i));
            T.IPInum(:,i) = ones(size(T.SN, 1), 1) * i;
            T.SN(:,i) = T.SN(:,1);
            T.isRep(:,i) = T.isRep(:,1);
        end
        T.IPI = reshape(T.IPI, size(T.IPI, 1) * size(T.IPI, 2), 1);
        T.IPInum = reshape(T.IPInum, size(T.IPInum, 1) * size(T.IPInum, 2), 1);
        T.SN = reshape(T.SN, size(T.SN, 1) * size(T.SN, 2), 1);
        T.isRep = reshape(T.isRep, size(T.isRep, 1) * size(T.isRep, 2), 1);
        
        % normalize IPI data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI'}, 'sub');
        
        subplot(2,3,1); title('IPI profile'); hold on;
        plt.line(T.IPInum, T.normIPI, 'split',T.isRep, 'errorbars','shade','style',isrepsty, 'leg',isrepleg);
        xlabel('IPI number'); ylabel('IPI (ms)'); set(gca,'fontsize',fs);
        xlim([0.5, 3.5]); ylim([150, 270]); axis square;
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Relationship with RTs
        T=tapply(D, {'SN', 'BN', 'isRep'},...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            {(D.IPI_1+D.IPI_2), 'nanmedian','name','first2ipi'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'first2ipi', 'RT'}, 'sub');
        
        subplot(2,3,[2,3]); title('Relationship with RTs'); hold on;
        plt.scatter(T.normfirst2ipi, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xlabel('first 2 IPIs (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs);
        xlim([150, 850]); ylim([100, 900]); %axis square;
        
        T=tapply(D, {'SN', 'BN', 'isRep'},...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI_1', 'IPI_2', 'IPI_3', 'RT'}, 'sub');
        
        subplot(2,3,4); title('1st IPI'); hold on;
        plt.scatter(T.normIPI_1, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xlabel('IPI (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs);
        xlim([0, 500]); ylim([150, 850]); axis square;
        
        subplot(2,3,5); title('2nd IPI'); hold on;
        plt.scatter(T.normIPI_1, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xlabel('IPI (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs);
        xlim([0, 500]); ylim([150, 850]); axis square;
        
        subplot(2,3,6); title('3rd IPI'); hold on;
        plt.scatter(T.normIPI_1, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg);
        xlabel('IPI (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs);
        xlim([0, 500]); ylim([150, 850]); axis square;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'N-1_IPI' % analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % add IPI info
        D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4],1,2);
        D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % create summary table for IPI profile
        T=tapply(D, {'SN', 'isRep'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        
        for i = 1:size(D.IPI, 2)
            T.IPI(:,i) = eval( sprintf('T.IPI_%d', i));
            T = rmfield(T,sprintf('IPI_%d', i));
            T.IPInum(:,i) = ones(size(T.SN, 1), 1) * i;
            T.SN(:,i) = T.SN(:,1);
            T.isRep(:,i) = T.isRep(:,1);
        end
        T.IPI = reshape(T.IPI, size(T.IPI, 1) * size(T.IPI, 2), 1);
        T.IPInum = reshape(T.IPInum, size(T.IPInum, 1) * size(T.IPInum, 2), 1);
        T.SN = reshape(T.SN, size(T.SN, 1) * size(T.SN, 2), 1);
        T.isRep = reshape(T.isRep, size(T.isRep, 1) * size(T.isRep, 2), 1);
        
        % normalize IPI data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI'}, 'sub');
        
        % open figure
        if nargin>1; figure('Name',sprintf('N-1 vs IPI analysis - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('N-1 vs IPI analysis - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,1); title(''); hold on;
        plt.line([T.IPInum], T.normIPI, 'split',[T.isRep], 'errorbars','shade','style',isrepsty, 'leg',isrepleg);
        xlabel('Transition number'); ylabel('Inter press interval (ms)'); set(gca,'fontsize',fs);
        
        xlim([0.5, 3.5]);
        ylim([130, 270]);
        axis square;
        
        % stats
        ttest(T.IPI(T.IPInum==1 & T.isRep==0), T.IPI(T.IPInum==1 & T.isRep==1), 2, 'paired');
        ttest(T.IPI(T.IPInum==2 & T.isRep==0), T.IPI(T.IPInum==2 & T.isRep==1), 2, 'paired');
        ttest(T.IPI(T.IPInum==3 & T.isRep==0), T.IPI(T.IPInum==3 & T.isRep==1), 2, 'paired');
        
        %         %-------------------------------------------------------------------------------------------------------------------------------------
        %         % create summary table for IPI profile
        %         T=tapply(D, {'SN', 'isRep', 'BN'},...
        %             {D.IPI_1, 'nanmedian','name','IPI_1'},...
        %             {D.IPI_2, 'nanmedian','name','IPI_2'},...
        %             {D.IPI_3, 'nanmedian','name','IPI_3'},...
        %             'subset',D.isError==0 & D.timingError==0 & D.exeType==1);
        %
        %         for i = 1:size(D.IPI, 2)
        %             T.IPI(:,i) = eval( sprintf('T.IPI_%d', i));
        %             T = rmfield(T,sprintf('IPI_%d', i));
        %             T.IPInum(:,i) = ones(size(T.SN, 1), 1) * i;
        %             T.SN(:,i) = T.SN(:,1);
        %             T.isRep(:,i) = T.isRep(:,1);
        %             T.BN(:,i) = T.BN(:,1);
        %         end
        %         T.IPI = reshape(T.IPI, size(T.IPI, 1) * size(T.IPI, 2), 1);
        %         T.IPInum = reshape(T.IPInum, size(T.IPInum, 1) * size(T.IPInum, 2), 1);
        %         T.SN = reshape(T.SN, size(T.SN, 1) * size(T.SN, 2), 1);
        %         T.isRep = reshape(T.isRep, size(T.isRep, 1) * size(T.isRep, 2), 1);
        %         T.BN = reshape(T.BN, size(T.BN, 1) * size(T.BN, 2), 1);
        %
        %         % normalize IPI data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'IPI'}, 'sub');
        %
        %         subplot(2,2,2); title(''); hold on;
        %         plt.line([T.BN T.IPInum], T.normIPI, 'split',[T.isRep], 'errorbars','shade','style',isrepsty, 'leg',isrepleg);
        %         xlabel('Transition number'); ylabel('Inter press interval (ms)'); set(gca,'fontsize',fs);
        %         ylim([120, 300]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        T = tapply(T, {'SN', 'IPInum'}, ...
            {T.IPI, 'nanmedian', 'name', 'IPI'}, ...
            {T.IPI, 'nanmedian', 'name', 'IPIs', 'subset',T.isRep==0}, ...
            {T.IPI, 'nanmedian', 'name', 'IPIr', 'subset',T.isRep==1});
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPIs', 'IPIr', 'IPI'}, 'sub');
        
        %norms = (T.IPIs ./ T.IPI) * 100;
        %normr = (T.IPIr ./ T.IPI) * 100;
        
        subplot(2,2,2); title(''); hold on;
        %plt.box(T.IPInum, (norms - normr), 'style',boxplotsty);
        %hold on; plt.line(T.IPInum, (norms - normr), 'style',lightgraysty);
        plt.box(T.IPInum, ( (T.IPIs-T.IPIr) ./ T.IPI ) * 100, 'style',boxplotsty);
        hold on; plt.line(T.IPInum, ( (T.IPIs-T.IPIr) ./ T.IPI ) * 100, 'style',lightgraysty);
        xlabel('Transition number'); ylabel('Repetition difference (% of IPI)'); set(gca,'fontsize',fs); axis square;
        xlim([0.5, 3.5]); ylim([-14 31]);
        drawline(0, 'dir','horz', 'linestyle',':');
        
        % stats
        diff2 = T.IPIs(T.IPInum==2) - T.IPIr(T.IPInum==2);
        diff3 = T.IPIs(T.IPInum==3) - T.IPIr(T.IPInum==3);
        ttest(diff2, diff3, 2, 'paired');
        
        diff23 = mean( [( (T.IPIs(T.IPInum==2)-T.IPIr(T.IPInum==2)) ./ T.IPI(T.IPInum==2) ) * 100, ( (T.IPIs(T.IPInum==3)-T.IPIr(T.IPInum==3)) ./ T.IPI(T.IPInum==3) ) * 100], 2);
        diff1  = ( (T.IPIs(T.IPInum==1)-T.IPIr(T.IPInum==1)) ./ T.IPI(T.IPInum==1) ) * 100;
        ttest(diff23, diff1, 2, 'paired');
        
        %         %-------------------------------------------------------------------------------------------------------------------------------------
        %         % N-1 IPI
        %         T=tapply(D, {'SN', 'isRep', 'nm1'},...
        %             {D.IPI_1, 'nanmedian','name','IPI_1'},...
        %             {D.IPI_2, 'nanmedian','name','IPI_2'},...
        %             {D.IPI_3, 'nanmedian','name','IPI_3'},...
        %             'subset', D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        %
        %         % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'IPI_1', 'IPI_2', 'IPI_3'}, 'sub');
        %
        %         subplot(2,2,2); hold on;
        %         plt.box([[ones(numel(T.nm1),1)*1; ones(numel(T.nm1),1)*2; ones(numel(T.nm1),1)*3], [T.nm1;T.nm1;T.nm1]], [T.IPI_1;T.IPI_2;T.IPI_3], 'split',[T.isRep; T.isRep; T.isRep], 'style',isrepstybox, 'leg',isrepleg, 'leglocation','northeast');
        %         axis square;
        %         %ylim([130, 270]);
        %         xt = xticks; xticks([sum(xt(1:2))/2 sum(xt(3:4))/2 sum(xt(5:6))/2 sum(xt(7:8))/2 sum(xt(9:10))/2 sum(xt(11:12))/2]); xticklabels(repmat({'No-go', 'Go'}, 1,3));
        %         xlabel('Previous trial'); ylabel('IPI (ms)'); set(gca,'fontsize',fs);
        %
        %         % stats
        %         T.ANOVA_IPI1 = anovaMixed(T.IPI_1, T.SN,'within', [T.nm1, T.isRep], {'n-1','isRep'});
        %         T.ANOVA_IPI2 = anovaMixed(T.IPI_2, T.SN,'within', [T.nm1, T.isRep], {'n-1','isRep'});
        %         T.ANOVA_IPI3 = anovaMixed(T.IPI_3, T.SN,'within', [T.nm1, T.isRep], {'n-1','isRep'});
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1 IPI
        % create summary table for IPI profile
        T=tapply(D, {'SN', 'isRep', 'nm1'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset',D.isError==0 & D.timingError==0 & D.exeType==1 & D.repNum<2);
        
        for i = 1:size(D.IPI, 2)
            T.IPI(:,i) = eval( sprintf('T.IPI_%d', i));
            T = rmfield(T,sprintf('IPI_%d', i));
            T.IPInum(:,i) = ones(size(T.SN, 1), 1) * i;
            T.SN(:,i) = T.SN(:,1);
            T.isRep(:,i) = T.isRep(:,1);
            T.nm1(:,i) = T.nm1(:,1);
        end
        T.IPI       = reshape(T.IPI, size(T.IPI, 1) * size(T.IPI, 2), 1);
        T.IPInum    = reshape(T.IPInum, size(T.IPInum, 1) * size(T.IPInum, 2), 1);
        T.SN        = reshape(T.SN, size(T.SN, 1) * size(T.SN, 2), 1);
        T.isRep     = reshape(T.isRep, size(T.isRep, 1) * size(T.isRep, 2), 1);
        T.nm1       = reshape(T.nm1, size(T.nm1, 1) * size(T.nm1, 2), 1);
        
        % normalize IPI data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI'}, 'sub');
        
        subplot(2,2,3); title('No-Go'); hold on;
        plt.line(T.IPInum, T.normIPI, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg, 'subset',T.nm1==1);
        xlabel('Transition number'); ylabel('Inter press interval (ms)'); set(gca,'fontsize',fs);
        xlim([0.5, 3.5]);
        ylim([145, 245]);
        axis square;
        
        subplot(2,2,4); title('Go'); hold on;
        plt.line(T.IPInum, T.normIPI, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg, 'subset',T.nm1==2);
        xlabel('Transition number'); ylabel('Inter press interval (ms)'); set(gca,'fontsize',fs);
        xlim([0.5, 3.5]);
        ylim([145, 245]);
        axis square;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'repEff_N-1_smt' % analysis of N-1 trial go/no-go on repetition effect, only sequences with same middle transition (smt)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('smt - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('smt - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % rep eff
        T = tapply(D, {'SN', 'smt'}, ...
            {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        subplot(2,2,1:2); title('Repetition effect'); hold on;
        [~, ~, ~] = plt.bar(T.smt, T.normET, 'style',blacksty);
        xlabel(''); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT', 'Repetition'}); ylim([600, 900]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1
        T = tapply(D, {'SN', 'smt', 'nm1'}, ...
            {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0 & D.nm1>0 );%& D.repNum<2);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        subplot(2,2,3:4); title('N-1'); hold on;
        [x, ~, ~] = plt.bar(T.smt, T.normET, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        xlabel(''); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT', 'Repetition'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2]); ylim([600, 900]);
        
        % stats
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'N-1_IPI_smt' % analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, only sequences with same middle transition (smt)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('N-1 IPI analysis smt - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('N-1 vs IPI analysis smt - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1 IPI smt
        T=tapply(D, {'SN', 'smt'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI_1', 'IPI_2', 'IPI_3'}, 'sub');
        
        subplot(2,3,1); title('IPI 1'); hold on;
        [~, ~, ~] = plt.bar(T.smt, T.normIPI_1, 'style',blacksty);
        xlabel(''); ylabel('First IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT', 'Repetition'}); ylim([150, 280]);
        
        subplot(2,3,2); title('IPI 2'); hold on;
        [~, ~, ~] = plt.bar(T.smt, T.normIPI_2, 'style',blacksty);
        xlabel(''); ylabel('Second IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT', 'Repetition'}); ylim([150, 280]);
        
        subplot(2,3,3); title('IPI 3'); hold on;
        [~, ~, ~] = plt.bar(T.smt, T.normIPI_3, 'style',blacksty);
        xlabel(''); ylabel('Third IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT', 'Repetition'}); ylim([150, 280]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1 IPI smt
        T=tapply(D, {'SN', 'smt', 'nm1'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1 & D.nm1>0 );%& D.repNum<2);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI_1', 'IPI_2', 'IPI_3'}, 'sub');
        
        subplot(2,3,4); title('IPI 1'); hold on;
        [x, ~, ~] = plt.bar(T.smt, T.normIPI_1, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        xlabel(''); ylabel('First IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT', 'Repetition'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2]); ylim([150, 280]);
        
        subplot(2,3,5); title('IPI 2'); hold on;
        [x, ~, ~] = plt.bar(T.smt, T.normIPI_2, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        xlabel(''); ylabel('Second IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT', 'Repetition'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2]); ylim([150, 280]);
        
        subplot(2,3,6); title('IPI 3'); hold on;
        [x, ~, ~] = plt.bar(T.smt, T.normIPI_3, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        xlabel(''); ylabel('Third IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT', 'Repetition'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2]); ylim([150, 280]);
        
        % stats
        %         T.ANOVA_IPI1 = anovaMixed(T.IPI_1, T.SN,'within', [T.nm1, T.smt], {'n-1','smt'});
        %         T.ANOVA_IPI2 = anovaMixed(T.IPI_2, T.SN,'within', [T.nm1, T.smt], {'n-1','smt'});
        %         T.ANOVA_IPI3 = anovaMixed(T.IPI_3, T.SN,'within', [T.nm1, T.smt], {'n-1','smt'});
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'sft_ET' % analysis of same finger transition different sequences (sft)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % combine switch types
        D.sft(D.sft==-1 | D.sft==0) = 0;
        
        % open figure
        if nargin>1; figure('Name',sprintf('sft - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('sft - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        %
        T = tapply(D, {'SN', 'sft'}, ...
            {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        subplot(2,2,1);
        %[x1, y1, e1] = plt.bar(T.sft, T.normET, 'style',blacksty);
        [~, y, ~] = plt.bar(T.sft, T.ET, 'style',boxplotsty, 'dir','horz');
        xlabel('Execution time (ms)'); ylabel(''); set(gca,'fontsize',fs); axis square;
        %yticklabels({'Swc', '0ST', 'S1T', 'S2T', 'S3T', 'S12T', 'S23T', 'Rep'}); xlim([660, 840]);
        yticklabels({'Swc', 'S1T', 'S2T', 'S3T', 'S12T', 'S23T', 'Rep'}); xlim([660, 840]);
        %yticklabels({'Swc', 'S1T', 'S3T', 'Rep'}); xlim([660, 840]);
        drawline(y(1), 'dir','vert', 'linestyle',':', 'color','k');
        %drawline(y(2), 'dir','vert', 'linestyle',':', 'color','k');
        drawline(y(end), 'dir','vert', 'linestyle',':', 'color','k');
        
        % stats
        
        [tab] = pivottable(T.SN, T.sft, T.ET, 'sum');
        swc = tab(:,1);
        s1t = tab(:,2);
        s2t = tab(:,3);
        s3t = tab(:,4);
        %s12 = tab(:,5);
        %s23 = tab(:,6);
        %rep = tab(:,7);
        
        ttest(s1t(~isnan(s1t)), swc(~isnan(s1t)), 2, 'paired');
        ttest(s2t(~isnan(s2t)), swc(~isnan(s2t)), 2, 'paired');
        ttest(s3t(~isnan(s3t)), swc(~isnan(s3t)), 2, 'paired');
        
        % combine switch types
        D.sft(D.sft==1 | D.sft==12) = 1;
        D.sft(D.sft==23 | D.sft==3) = 3;
        
        T = tapply(D, {'SN', 'sft'}, ...
            {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0);
        
        [tab] = pivottable(T.SN, T.sft, T.ET, 'sum');
        swc = tab(:,1);
        s1t = tab(:,2);
        s3t = tab(:,4);
        
        ttest(s1t(~isnan(s1t)), swc(~isnan(s1t)), 2, 'paired');
        ttest(s3t(~isnan(s3t)), swc(~isnan(s3t)), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        %
        %         T = tapply(D, {'SN', 'sft', 'nm1'}, ...
        %             {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1}, ... & D.sft>=0}, ...
        %             'subset', D.isError==0 & D.timingError==0 & D.repNum<2);
        %
        %         % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'ET'}, 'sub');
        %
        %         subplot(2,2,2);
        %         plt.bar(T.sft, T.normET, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','southeast', 'dir','horz');
        %         xlabel('ET (ms)'); ylabel(''); set(gca,'fontsize',fs); axis square;
        %         yt = yticks; yticks([sum(yt(1:2))/2, sum(yt(3:4))/2, sum(yt(5:6))/2, sum(yt(7:8))/2, sum(yt(9:10))/2, sum(yt(11:12))/2, sum(yt(13:14))/2, sum(yt(15:16))/2]);
        %         yticklabels({'Swc', '0ST', 'S1T', 'S2T', 'S3T', 'S12T', 'S23T', 'Rep'}); xlim([660, 840]);
        
        % stats
        %[t,p] = ttest(T.ET(T.sft==0), T.ET(T.sft==3), 2, 'independent');
        %text(x1(1)-.3, max(y1)+sum(e1)/2, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',8);
        
        %[t,p] = ttest(T.ET(T.sft==0), T.ET(T.sft==12), 2, 'independent');
        %text(x1(3)-.3, max(y1)+sum(e1)/2, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',8);
        
        %[t,p] = ttest(T.ET(T.sft==0), T.ET(T.sft==23), 2, 'independent');
        %text(x1(5)-.3, max(y1)+sum(e1)/2, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',8);
        
        %[t,p] = ttest(T.ET(T.sft==0), T.ET(T.sft==123), 2, 'independent');
        %text(x1(7)-.3, max(y1)+sum(e1)/2, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',8);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        %         % N-1
        %         T = tapply(D, {'SN', 'sft', 'nm1'}, ...
        %             {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1 & D.sft>=0}, ...
        %             'subset', D.isError==0 & D.timingError==0 & D.nm1>0);
        %
        %         % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'ET'}, 'sub');
        %
        %         subplot(2,2,3:4); title('N-1'); hold on;
        %         [x, ~, ~] = plt.bar(T.sft, T.normET, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        %         xlabel(''); ylabel('ET (ms)'); set(gca,'fontsize',fs); %axis square;
        %         xticklabels({'0ST', 'S1T', 'S2T', 'S3T', 'S12T', 'S23T', 'Rep'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2, (x(7)+x(8))/2, (x(9)+x(10))/2, (x(11)+x(12))/2, (x(13)+x(14))/2]); ylim([600, 900]);
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'sft_RT' % analysis of same finger transition different sequences (sft)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % define same finger (sf) in different sequences code
        D.sf = zeros(numel(D.SN), 1);
        D.sf(D.isRep == 1 & D.sft == 123) = 4;
        D.sf(D.isRep == 0 & D.sft == 12) = 3;
        D.sf(D.isRep == 0 & D.sft == 1) = 2;
        D.sf(D.isRep == 0 & D.sft ~= 12 & D.sft ~= 1 & D.sff == 1) = 1;
        
        %         % open figure
        %         if nargin>1; figure('Name',sprintf('sft - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('sft - group (N=%d)',ns)); end
        %         set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        %
        T = tapply(D, {'SN', 'sf'}, ...
            {D.RT, 'nanmedian', 'name', 'RT', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        subplot(2,2,3);
        [~,y,~] = plt.bar(T.sf, T.normRT, 'style',boxplotsty, 'dir','horz');
        xlabel('Reaction time (ms)'); ylabel(''); set(gca,'fontsize',fs); axis square;
        yticklabels({'Swc', 'S1F', 'S12F', 'S123F', 'Rep'}); xlim([390, 510]);
        drawline(y(1), 'dir','vert', 'linestyle',':', 'color','k');
        %drawline(y(2), 'dir','vert', 'linestyle',':', 'color','k');
        drawline(y(end), 'dir','vert', 'linestyle',':', 'color','k');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        %
        T = tapply(D, {'SN', 'sf', 'nm1'}, ...
            {D.RT, 'nanmedian', 'name', 'RT', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0 );%& D.repNum<2);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        subplot(2,2,4);
        plt.bar(T.sf, T.normRT, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','southeast', 'dir','horz');
        xlabel('RT (ms)'); ylabel(''); set(gca,'fontsize',fs); axis square;
        yt = yticks; yticks([sum(yt(1:2))/2, sum(yt(3:4))/2, sum(yt(5:6))/2, sum(yt(7:8))/2, sum(yt(9:10))/2]);%, sum(yt(11:12))/2, sum(yt(13:14))/2, sum(yt(15:16))/2]);
        yticklabels({'Swc', 'S1F', 'S12F', 'S123F', 'Rep'}); xlim([390, 510]);
        
        % stats
        %[t,p] = ttest(T.RT(T.sft==0), T.RT(T.sft==3), 2, 'independent');
        %text(x1(1)-.3, max(y1)+sum(e1)/2, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',8);
        
        %[t,p] = ttest(T.RT(T.sft==0), T.RT(T.sft==12), 2, 'independent');
        %text(x1(3)-.3, max(y1)+sum(e1)/2, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',8);
        
        %[t,p] = ttest(T.RT(T.sft==0), T.RT(T.sft==23), 2, 'independent');
        %text(x1(5)-.3, max(y1)+sum(e1)/2, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',8);
        
        %[t,p] = ttest(T.RT(T.sft==0), T.RT(T.sft==123), 2, 'independent');
        %text(x1(7)-.3, max(y1)+sum(e1)/2, sprintf('t = %2.2f, p < %2.2f', t, p), 'fontsize',8);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        %         % N-1
        %         T = tapply(D, {'SN', 'sft', 'nm1'}, ...
        %             {D.RT, 'nanmedian', 'name', 'RT', 'subset', D.exeType == 1 & D.sft>=0}, ...
        %             'subset', D.isError==0 & D.timingError==0 & D.nm1>0);
        %
        %         % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'RT'}, 'sub');
        %
        %         subplot(2,2,3:4); title('N-1'); hold on;
        %         [x, ~, ~] = plt.bar(T.sft, T.normRT, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        %         xlabel(''); ylabel('RT (ms)'); set(gca,'fontsize',fs); %axis square;
        %         xticklabels({'0ST', 'S1T', 'S2T', 'S3T', 'S12T', 'S23T', 'Rep'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2, (x(7)+x(8))/2, (x(9)+x(10))/2, (x(11)+x(12))/2, (x(13)+x(14))/2]); ylim([600, 900]);
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'N-1_IPI_sft' % analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, same transition(s) in different sequences (sft)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('N-1 IPI analysis sft - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('N-1 vs IPI analysis sft - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1 IPI sft
        T=tapply(D, {'SN', 'sft'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI_1', 'IPI_2', 'IPI_3'}, 'sub');
        
        subplot(2,3,1); title('IPI 1'); hold on;
        [~, ~, ~] = plt.bar(T.sft, T.normIPI_1, 'style',blacksty);%, 'subset',ismember(T.sft, [0,1,12,123]));
        xlabel(''); ylabel('First IPI (ms)'); set(gca,'fontsize',fs); axis square;
        %xticklabels({'Swc', 'S1T', 'S12T', 'Rep'}); ylim([135, 295]);
        
        subplot(2,3,2); title('IPI 2'); hold on;
        [~, ~, ~] = plt.bar(T.sft, T.normIPI_2, 'style',blacksty);%, 'subset',ismember(T.sft, [0,2,12,23,123]));
        xlabel(''); ylabel('Second IPI (ms)'); set(gca,'fontsize',fs); axis square;
        %xticklabels({'Swc', 'S2T', 'S12T', 'S23T', 'Rep'}); ylim([135, 295]);
        
        subplot(2,3,3); title('IPI 3'); hold on;
        [~, ~, ~] = plt.bar(T.sft, T.normIPI_3, 'style',blacksty);%, 'subset',ismember(T.sft, [0,3,23,123]));
        xlabel(''); ylabel('Third IPI (ms)'); set(gca,'fontsize',fs); axis square;
        %xticklabels({'Swc', 'S3T', 'S23T', 'Rep'}); ylim([135, 295]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1 IPI sft
        T=tapply(D, {'SN', 'sft', 'nm1'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1 & D.nm1>0 );%& D.repNum<2);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI_1', 'IPI_2', 'IPI_3'}, 'sub');
        
        subplot(2,3,4); title('IPI 1'); hold on;
        [x, ~, ~] = plt.bar(T.sft, T.normIPI_1, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast', 'subset',ismember(T.sft, [0,1,12,123]));
        xlabel(''); ylabel('First IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'S1T', 'S12T', 'Rep'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2, (x(7)+x(8))/2]); ylim([135, 295]);
        
        subplot(2,3,5); title('IPI 2'); hold on;
        [x, ~, ~] = plt.bar(T.sft, T.normIPI_2, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast', 'subset',ismember(T.sft, [0,2,12,23,123]));
        xlabel(''); ylabel('Second IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'S2T', 'S12T', 'S23T', 'Rep'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2, (x(7)+x(8))/2, (x(9)+x(10))/2]); ylim([135, 295]);
        
        subplot(2,3,6); title('IPI 3'); hold on;
        [x, ~, ~] = plt.bar(T.sft, T.normIPI_3, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast', 'subset',ismember(T.sft, [0,3,23,123]));
        xlabel(''); ylabel('Third IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'S3T', 'S23T', 'Rep'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2, (x(7)+x(8))/2]); ylim([135, 295]);
        
        % stats
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'repEff_N-1_smt_fls' % analysis of N-1 trial go/no-go on repetition effect, same middle transition in sequences with first and last press swapped (smt_fls)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('smt_fls - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('smt_fls - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % rep eff
        T = tapply(D, {'SN', 'smt_fls'}, ...
            {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        subplot(2,2,1:2); title('Repetition effect'); hold on;
        [~, ~, ~] = plt.bar(T.smt_fls, T.normET, 'style',blacksty);
        xlabel(''); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT FLS', 'Repetition'}); ylim([600, 950]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1
        T = tapply(D, {'SN', 'smt_fls', 'nm1'}, ...
            {D.ET, 'nanmedian', 'name', 'ET', 'subset', D.exeType == 1}, ...
            'subset', D.isError==0 & D.timingError==0 & D.nm1>0 );%& D.repNum<2);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        subplot(2,2,3:4); title('N-1'); hold on;
        [x, ~, ~] = plt.bar(T.smt_fls, T.normET, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        xlabel(''); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Switch', 'SMT FLS', 'Repetition'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2]); ylim([600, 950]);
        
        % stats
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'N-1_IPI_smt_fls' % analysis of N-1 trial go/no-go on repetition effect, separately for different IPIs, same middle transition in sequences with first and last press swapped (smt_fls)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('N-1 IPI analysis smt_fls - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('N-1 vs IPI analysis smt_fls - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1 IPI smt_fls
        T=tapply(D, {'SN', 'smt_fls'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI_1', 'IPI_2', 'IPI_3'}, 'sub');
        
        subplot(2,3,1); title('IPI 1'); hold on;
        [~, ~, ~] = plt.bar(T.smt_fls, T.normIPI_1, 'style',blacksty);
        xlabel(''); ylabel('First IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'SMT FLS', 'Rep'}); ylim([135, 315]);
        
        subplot(2,3,2); title('IPI 2'); hold on;
        [~, ~, ~] = plt.bar(T.smt_fls, T.normIPI_2, 'style',blacksty);
        xlabel(''); ylabel('Second IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'SMT FLS', 'Rep'}); ylim([135, 315]);
        
        subplot(2,3,3); title('IPI 3'); hold on;
        [~, ~, ~] = plt.bar(T.smt_fls, T.normIPI_3, 'style',blacksty);
        xlabel(''); ylabel('Third IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'SMT FLS', 'Rep'}); ylim([135, 315]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % N-1 IPI smt_fls
        T=tapply(D, {'SN', 'smt_fls', 'nm1'},...
            {D.IPI_1, 'nanmedian','name','IPI_1'},...
            {D.IPI_2, 'nanmedian','name','IPI_2'},...
            {D.IPI_3, 'nanmedian','name','IPI_3'},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1 & D.nm1>0 );%& D.repNum<2);
        
        % normalize data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI_1', 'IPI_2', 'IPI_3'}, 'sub');
        
        subplot(2,3,4); title('IPI 1'); hold on;
        [x, ~, ~] = plt.bar(T.smt_fls, T.normIPI_1, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        xlabel(''); ylabel('First IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'SMT FLS', 'Rep'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2]); ylim([135, 315]);
        
        subplot(2,3,5); title('IPI 2'); hold on;
        [x, ~, ~] = plt.bar(T.smt_fls, T.normIPI_2, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        xlabel(''); ylabel('Second IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'SMT FLS', 'Rep'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2]); ylim([135, 315]);
        
        subplot(2,3,6); title('IPI 3'); hold on;
        [x, ~, ~] = plt.bar(T.smt_fls, T.normIPI_3, 'split',T.nm1, 'style',nm1sty, 'leg',nm1leg, 'leglocation','northeast');
        xlabel(''); ylabel('Third IPI (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'Swc', 'SMT FLS', 'Rep'}); xticks([(x(1)+x(2))/2, (x(3)+x(4))/2, (x(5)+x(6))/2]); ylim([135, 315]);
        
        % stats
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'repEff_scatter' % analysis of relationship between ET, RT and repetition effect
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr3_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr3_all_data.mat'));
        end
        
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('repEff scatter - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('repEff scatter - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        T=tapply(D, {'SN'},...
            {D.RT,'nanmedian', 'name','RT'},...
            {D.RT,'nanmedian', 'name','RTrep', 'subset',D.isRep==1},...
            {D.RT,'nanmedian', 'name','RTswc', 'subset',D.isRep==0},...
            {D.ET,'nanmedian', 'name','ET'},...
            {D.ET,'nanmedian', 'name','ETrep', 'subset',D.isRep==1},...
            {D.ET,'nanmedian', 'name','ETswc', 'subset',D.isRep==0},...
            'subset', D.isError==0 & D.timingError==0 & D.exeType==1);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % RT-repEff RT
        subplot(2,2,1);
        plt.scatter(T.RT, T.RTswc-T.RTrep, 'style',darkgraysty, 'label',subj, 'regression','linear');%, 'subset',~ismember(T.SN,[1,7,14,15]));
        xlabel('Reaction time (ms)'); ylabel('Repetition difference on RT (ms)'); set(gca,'fontsize',fs); axis square;
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % ET-repEff ET
        subplot(2,2,2);
        [T.r2, T.b, T.t, T.p] = plt.scatter(T.ET, ((T.ETswc-T.ETrep)./T.ET)*100, 'style',darkgraysty, 'label',subj, 'regression','linear');%, 'subset',~ismember(T.SN,[1,7,14,15]));
        xlabel('Execution time (ms)'); ylabel('Repetition difference on ET (ms)'); set(gca,'fontsize',fs); axis square;
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % RT-repEff ET
        subplot(2,2,3);
        plt.scatter(T.RT, T.ETswc-T.ETrep, 'style',darkgraysty, 'label',subj);
        xlabel('Reaction time (ms)'); ylabel('Repetition difference on ET (ms)'); set(gca,'fontsize',fs); axis square;
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % ET-repEff RT
        subplot(2,2,4);
        plt.scatter(T.ET, T.RTswc-T.RTrep, 'style',darkgraysty, 'label',subj);
        xlabel('Execution time (ms)'); ylabel('Repetition difference on RT (ms)'); set(gca,'fontsize',fs); axis square;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    otherwise
        error('no such case!')
end
end
