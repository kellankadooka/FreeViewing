function [betweencorr, p, meanvar, dist] = TestBetweenSingleRange(stim, g1, g2, starti,endi, iterations, graph)
%Simplified version of test between to allow comparison of a single subject
%to one or many others
%%
clear all
stim = 'stim';
g1 = 'stim/7_006.csv';
g2 = '7';
starti = 1;
endi = 59992;
iterations = 1000;
graph = 1;


subjectfiles = dir(strcat(stim,'/*.csv'));

for i = 1:length(subjectfiles)
    subjectfiles(i).name = strcat(stim,'/',subjectfiles(i).name);
end

samples = length(starti:endi);

g1files = []; %Single subject
g2files = []; %Group of subjects
%Identify groups
for i=1:length(subjectfiles)
    if strcmp(subjectfiles(i).name,g1)
        g1files = [g1files i];
    end
    if subjectfiles(i).name(6) == g2
        g2files = [g2files i];
    end
end

%Pull gaze coordinates from group 1 into xpos and ypos
g1xpos = zeros(samples,length(g1files));
g1ypos = zeros(samples,length(g1files));
g1alltimes = zeros(samples,length(g1files));

for i=1:length(g1files)
gazedata = csvread(subjectfiles(g1files(i)).name);
g1xpos(:,i) = gazedata(starti:endi,3);   
g1ypos(:,i) = gazedata(starti:endi,4);
g1alltimes(:,i) = gazedata(starti:endi,6);
end

g1xpos(g1xpos == 0) = NaN;
g1ypos(g1ypos == 0) = NaN;

%Pull gaze coordinates from group 2 into xpos and ypos
g2xpos = zeros(samples,length(g2files));
g2ypos = zeros(samples,length(g2files));
g2alltimes = zeros(samples,length(g2files));
%%
for i=1:length(g2files)
gazedata = csvread(subjectfiles(g2files(i)).name);
g2xpos(:,i) = gazedata(starti:endi,3);   
g2ypos(:,i) = gazedata(starti:endi,4);
g2alltimes(:,i) = gazedata(starti:endi,5);
end

g2xpos(g2xpos == 0) = NaN;
g2ypos(g2ypos == 0) = NaN;

xypos=[g1xpos g1ypos g2xpos g2ypos];
index=isnan(mean(xypos,2))==0;
index = find(index == 1);
eyesopentogether = (length(index)/length(xypos)) %Gives time both eyes open (ratio)
%%
x1 = g1xpos(index,:);
y1 = g1ypos(index,:);
times1 = g1alltimes(index, :);
x2 = g2xpos(index,:);
y2 = g2ypos(index,:);
times2 = g2alltimes(index, :);

g1subj = length(g1files);
g2subj = length(g2files);
numsubj = g1subj + g2subj;

%BETWEEN GROUP CORRELATION
tempMean=mean(x2,2);
temp=corrcoef(tempMean,x1);
corrx2=temp(2,1);

tempMean=mean(y2,2);
temp=corrcoef(tempMean,y1);
corry2=temp(2,1);

betweencorr = (corrx2 + corry2)./2;
%%
%VARIANCE
varx = var(x1);
vary = var(y1);
meanvar = (varx + vary) ./ 2000;

%DISTANCE FROM CENTER
% for i = 1:length(x1)
%     distcenter(i) = distance(x1(i),y1(i),512, 384);
% end
% dist = nanmean(distcenter);

nulldis = zeros(iterations,1);
for j = 1:iterations
    newx = x2;
    newy = y2;
    %Randomly flip each subject
    for i = 1:length(g2files)
        if rand > .5
            newx(:,i) = flipud(newx(:,i));
            newy(:,i) = flipud(newy(:,i));
        end
    end
    
    %Start each scanpath from a random index
    for i = 1:length(g2files)
        startindex = ceil(rand*length(newx));
        xorig = newx(:,i);
        yorig = newy(:,i);
        newx(startindex:end,i) = xorig(1:length(xorig)-startindex+1);
        newx(1:startindex-1,i) = xorig(length(xorig)-startindex+2:end);
        newy(startindex:end,i) = yorig(1:length(yorig)-startindex+1);
        newy(1:startindex-1,i) = yorig(length(yorig)-startindex+2:end);
    end
        
    tempMean=mean(newx,2);
    temp=corrcoef(tempMean,x1);
    nullcorrx2=temp(2,1);

    tempMean=mean(newy,2);
    temp=corrcoef(tempMean,y1);
    nullcorry2=temp(2,1);
    
    nulldis(j) = (nullcorrx2 + nullcorry2)./2;
end

p = sum(betweencorr < nulldis)/iterations;
if p > .5
    p = 1 - p;
end

if graph == 1
    N = hist(nulldis, 50);
    hist(nulldis, 50)
    hold on
    plot(zeros(max(N),1) + betweencorr,1:max(N),'r')
    hold off
end