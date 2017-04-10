function [xvar, yvar, meanvar, top, resulttop, indtop, resultattop] = Variability_Complete(stim, group)
%%
stim = 7;
group = 3;

ppts = csvread(strcat(num2str(stim),'_inclusion.csv'));

id = ppts(:,1);
%include = ppts(:,3);
%inc_id = id(include == 1);
group_num = ppts(:,4);      %select column from inclusion table and value of group 
%child_id = inc_id(group_num==group);
ppt_id = id(group_num == group);

xdata = zeros(3600,numel(ppt_id));
ydata = xdata;

for i = 1:numel(ppt_id)
    filename = sprintf('Variability/%d_%03d.csv',stim,ppt_id(i));
    temp_data = ReadVariability(filename);
    xdata(:,i) = temp_data(:,2);
    ydata(:,i) = temp_data(:,3);
end

xdata(xdata == -999) = NaN;
ydata(ydata == -999) = NaN;

xvar = nanvar(xdata,[],2);
yvar = nanvar(ydata,[],2);

meanvar = mean([xvar yvar],2);

plot(meanvar)

top=flipud(unique(sort(meanvar)));
resulttop=top(1:90);         %top 2.5%
indtop=find(meanvar>=top(90))  ;    %their indices
resultattop=flipud(sortrows([meanvar(indtop) indtop],1));

resultattop(:,[1,2])= resultattop(:,[2,1]);
hold on

plot((resultattop(:,1)), (resultattop(:,2)), 'g*')


bottom=flipud(unique(sort(meanvar)));
resultbottom=bottom(3510:3600);         %top 2.5%
indbottom=find(meanvar<=bottom(3510))  ;    %their indices
resultatbottom=flipud(sortrows([meanvar(indbottom) indbottom],1));

resultatbottom(:,[1,2])= resultatbottom(:,[2,1]);
hold on
plot((resultatbottom(:,1)), (resultatbottom(:,2)), 'g*')


%butter
SF = 30; % sampling frequency
NF = SF/2; % Nyquist frequency
CF = 3; % Cut-off frequency
% initalize normalized cut-off frequency Wn with a value between 0 and 1
Wn = CF/NF; % == 9Hz/30Hz = 0.3
% run butter
[b,a] = butter(2, Wn, 'low'); % 2nd order, low-pass filter
%filter data using filtfilt (zero-phase digital filtering)
var_filtered = filtfilt(b,a,meanvar);SF = 30; % sampling frequency

hold on
plot(var_filtered)

%top & bottom on var_filtered
topf=flipud(unique(sort(var_filtered)));
resulttopf=topf(1:90);         %top 2.5%
indtopf=find(var_filtered>=topf(90))  ;    %their indices
resultattopf=flipud(sortrows([var_filtered(indtopf) indtopf],1));

resultattopf(:,[1,2])= resultattopf(:,[2,1]);
hold on

plot((resultattopf(:,1)), (resultattopf(:,2)), 'r*')


bottomf=flipud(unique(sort(var_filtered)));
resultbottomf=bottomf(3510:3600);         %top 2.5%
indbottomf=find(var_filtered<=bottomf(3510))  ;    %their indices
resultatbottomf=flipud(sortrows([var_filtered(indbottomf) indbottomf],1));

resultatbottomf(:,[1,2])= resultatbottomf(:,[2,1]);
hold on
plot((resultatbottomf(:,1)), (resultatbottomf(:,2)), 'r*')

sorttop = resultattopf(:,1);
sorttop = sort(sorttop, 'descend')

sortbottom = resultatbottomf(:,1);
sortbottom = sort(sortbottom,'descend')

