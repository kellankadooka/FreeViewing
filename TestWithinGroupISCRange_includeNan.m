function [testcorr pval eyesopen alldata] = TestWithinGroupISCRange(stim, starti,endi, iterations, graph)

%%
clear all
%%
stim = 'stim';
starti = 1;
endi = 59992;
iterations = 1000;
graph = 1;
%%
subjectfiles = dir(strcat(stim,'/*.csv'));

for i = 1:length(subjectfiles)
    subjectfiles(i).name = strcat(stim,'/',subjectfiles(i).name);
end
%%
numsubj = length(subjectfiles);
samples = length(starti:endi);

xpos = zeros(samples,numsubj);
ypos = zeros(samples,numsubj);
alltimes = zeros(samples,numsubj);
%%
%Pull gaze coordinates from all subjects into just xpos and ypos
for i=1:length(subjectfiles)
gazedata = csvread(subjectfiles(i).name);
xpos(:,i) = gazedata(starti:endi,3);   
ypos(:,i) = gazedata(starti:endi,4);
alltimes(:,i) = gazedata(1:length(alltimes),6);
end
%%
xpos(xpos <= 0) = NaN;
xpos(xpos > 1280) = NaN;
ypos(ypos <= 0) = NaN;
ypos(ypos > 720) = NaN;

xypos=[xpos ypos];
%%
index=isnan(mean(xypos,2))==0;
%%
index = find(index == 1);
eyesopentogether = (length(index)/length(xypos)) %Gives time both eyes open (ratio)
%%
x = xpos(index,:);
%%
y = ypos(index,:);
times = alltimes(index, :);
%%
for i=1:numsubj
    tempMean=((nansum(xpos,2))-xpos(:,i))/((sum(xpos==xpos,2))-1);
    temp=corrcoef(tempMean,xpos(:,i));
    corrx(i)=temp(2,1); %Correlation within x
    varx(i) = nanvar(xpos(:,i));
end
%%
for i=1:numsubj
    tempMean=((nansum(y,2))-y(:,i))/(numsubj-1);
    temp=corrcoef(tempMean,y(:,i));
    corry(i)=temp(2,1); %Correlation within y
    vary(i) = var(y(:,i));
end


corr=mean([mean(corrx) mean(corry)]);
%%
nulldis = zeros(iterations,1);

for j = 1:iterations
    newx = x;
    newy = y;
    %Randomly flip each subject
    for i = 1:numsubj
        if rand > .5
            newx(:,i) = flipud(newx(:,i));
            newy(:,i) = flipud(newy(:,i));
        end
    end
    
    %Start each scanpath from a random index
    for i = 1:numsubj
        startindex = ceil(rand*length(newx));
        xorig = newx(:,i);
        yorig = newy(:,i);
        newx(startindex:end,i) = xorig(1:length(xorig)-startindex+1);
        newx(1:startindex-1,i) = xorig(length(xorig)-startindex+2:end);
        newy(startindex:end,i) = yorig(1:length(yorig)-startindex+1);
        newy(1:startindex-1,i) = yorig(length(yorig)-startindex+2:end);
    end
        
    for i=1:numsubj
        tempMean=((sum(newx,2))-newx(:,i))/(numsubj-1);
        temp=corrcoef(tempMean,newx(:,i));
        nullcorrx(i)=temp(2,1); %Correlation within x
    end

    for i=1:numsubj
        tempMean=((sum(newy,2))-newy(:,i))/(numsubj-1);
        temp=corrcoef(tempMean,newy(:,i));
        nullcorry(i)=temp(2,1); %Correlation within y
    end
    
    nulldis(j) = mean([mean(nullcorrx) mean(nullcorry)]);
end
    

%%
if graph == 1
    N = hist(nulldis, 50);
    hist(nulldis, 50)
    hold on
    plot(zeros(max(N),1) + corr,1:max(N),'r')
    hold off
end
p = sum(corr < nulldis)/iterations;

if p > .5
    p = 1 - p;
end

if graph == 2 %TWO EXEMPLARS
    time = alltimes(:,2)./1000;
    xdata1 = xpos(:,7);
    xdata2 = xpos(:,9);
    plot(time,xdata1,'b--','LineWidth',1);
    hold on
    plot(time,xdata2,'r--','LineWidth',1);
    hold off
end
    
if graph == 4 %VARIANCE
    time = times(:,2)./1000;
    stdx = smooth(std(y,0,2),120);
    plot(time,stdx,'LineWidth',1);
    stdxsort = sort(stdx);
    upperlim = stdxsort(round(length(stdx) * .95));
    lowerlim = stdxsort(round(length(stdx) * .05));
    hold on
    plot(time(stdx < lowerlim),stdx(stdx < lowerlim),'g*')
    plot(time(stdx > upperlim),stdx(stdx > upperlim),'r*')
end

corrx;
corry;
meancorr = (corrx + corry)./2;
meancorr';
meanvar = (varx + vary)./2000;
meanvar';

if graph == 5
    alldata(:,:,1) = xpos';
    alldata(:,:,2) = ypos';
else
    alldata(:,:,1) = x';
    alldata(:,:,2) = y';
end
pval = p;
testcorr = meancorr';
eyesopen = eyesopentogether;