function [D] = sr3_trial(MOV, D, fig, fig_name, varargin)
%% function [D] = sr3_trial(MOV, D, fig, fig_name, varargin)
% Trial routine for SequenceRepetition experiment 3
% called by the Subject routine (sr3_subj.m)

%% extract data
if (isempty(MOV))
    return;
end;
state=MOV(:,1);
%absTime=MOV(:,2);
time=MOV(:,3);
%time=2:2:(size(MOV(:,1),1))*2;
force=smooth_kernel(MOV(:,4:end),4);
RH=[6 7 8 9 10]; %right hand column indices
LH=[1 2 3 4 5]; %left hand column indices

%%
fNames=fieldnames(D);
resp=[]; %which presses?

%%
for i=1:length(fNames)
    if length(fNames{i})>8
        if strcmp(fNames{i}(1:8),'response')
            eval(['resp=[resp,D.',fNames{i},'];']);
        end
    end
end
D.numPress=sum(resp>0); %number of presses
pressTime=zeros(1,D.numPress);
for press=1:D.numPress
    pressNum=['D.pressTime',num2str(press)];
    pressTime(press)=eval(pressNum); %time of presses
end

%% Display trial
if (fig>0)
    figure('Name',fig_name)
    subplot(2,1,1)
    if D.hand==1
        plot(time,force(:,LH),'LineWidth',2);
        title('Force traces for LEFT hand presses','FontSize',20);
    elseif D.hand==2
        plot(time,force(:,RH),'LineWidth',2);
        title('Force traces for RIGHT hand presses','FontSize',20);
    end
    xlabel('Time (ms)'); ylabel('Force (N)'); set(gca,'FontSize',20);
    hold on;
    drawline(pressTime,'dir','vert');
    legend({'Thumb','Index','Middle','Ring','Little'},'FontSize',20)
    hold off;
    subplot(2,1,2)
    plot(time,state,'LineWidth',2);
    ylim([1 7]);
    xlabel('Time (ms)'); ylabel('State'); set(gca,'FontSize',20);
    hold on;
    drawline(pressTime,'dir','vert'); legend state
    KbWait;
end