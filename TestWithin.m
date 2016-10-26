function [corr, var, p, eyesopen, available] = TestWithin(stim, group, starti, endi)

ppts = csvread(strcat(num2str(stim),'_inclusion.csv'));

id = ppts(:,1);
group_num = ppts(:,2);

child_id = id(group_num == group);

for i = 1:numel(child_id)
    g1 = sprintf('data/%d_%03d.csv',stim,child_id(i));
    for j = 1:numel(child_id)
        g2 = sprintf('data/%d_%03d.csv',stim,child_id(j));
        [betweencorr(j), p(j), meanvar(j), eyesopentogether(j), available(j)] = TestBetweenSingleRange(g1, g2, starti,endi, 0, 0);
        if i == j
            betweencorr(j) = NaN;
            p(j) = NaN;
            meanvar(j) = NaN;
            eyesopentogether(j) = NaN;
            available(j) = NaN;
        end
    end
    corr(i) = nanmean(betweencorr);
    var(i) = nanmean(meanvar);
    eyesopen(i) = nanmean(eyesopentogether);
    available(i) = nanmean(available);    
end