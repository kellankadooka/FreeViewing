clear all 
stim = 7;
ppts = csvread(strcat(num2str(stim),'_inclusion.csv'));

id = ppts(:,1);
group = ppts(:,2);

adult_id = id(group == 1);
child_id = id(group == 2);

child_between = zeros(length(child_id),1);
%%
for i = 1:numel(child_id)
    g1 = sprintf('data/%d_%03d.csv',stim,child_id(i))
    for j = 1:numel(adult_id)
        g2 = sprintf('data/%d_%03d.csv',stim,adult_id(j))
        [betweencorr(j), p(j), meanvar(j), eyesopentogether(j), available(j)] = TestBetweenSingleRange(g1, g2, 1,59992, 0, 0);
    end
    child_between(i) = mean(betweencorr);
    child_p(i) = mean(p);
    child_meanvar(i) = mean(meanvar);
    child_eyesopen(i) = mean(eyesopentogether);
    child_available(i) = mean(available);
end
%%
for i = 1:numel(child_id)
    g1 = sprintf('data/%d_%03d.csv',stim,child_id(i));
    for j = 1:numel(child_id)
        g2 = sprintf('data/%d_%03d.csv',stim,child_id(j))
        [betweencorr(j), p(j), meanvar(j), eyesopentogether(j)] = TestBetweenSingleRange(g1, g2, 1,59992, 0, 0);
        if i == j
            betweencorr(j) = NaN;
        end
    end
    child_within(i) = nanmean(betweencorr);
end


%%
for i = 1:numel(adult_id)
    g1 = sprintf('data/%d_%03d.csv',stim,adult_id(i));
    for j = 1:numel(adult_id)
        g2 = sprintf('data/%d_%03d.csv',stim,adult_id(j))
        [betweencorr(j), p(j), meanvar(j), eyesopentogether(j)] = TestBetweenSingleRange(g1, g2, 1,59992, 0, 0);
        if i == j
            betweencorr(j) = NaN;
        end
    end
    adult_within(i) = nanmean(betweencorr);
end

