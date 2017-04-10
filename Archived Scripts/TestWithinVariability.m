function [xvar, yvar, meanvar] = TestWithinVariability(stim, group)


ppts = csvread(strcat(num2str(stim),'_inclusion.csv'));

id = ppts(:,1);
%include = ppts(:,3);
%inc_id = id(include == 1);
group_num = ppts(:,2);      %select column from inclusion table and value of group 
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
