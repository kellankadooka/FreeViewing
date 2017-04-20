function [data] = TestGazeROISaliency(stim, duration, graphs, channel)
%
%clear all
addpath('util')
addpath('ISC_saliency')
% graphs = 0;
% stim = 'sesame';
%duration = 59;
%channel = 'orientation';

subjectfiles = dir(strcat(stim,'/*.csv'));

for i = 1:length(subjectfiles)
    subjectfiles(i).name = strcat(stim,'/',subjectfiles(i).name);
end

numsubj = length(subjectfiles);
samples = duration * 120;

xpos = zeros(samples,numsubj);
ypos = zeros(samples,numsubj);
alltimes = zeros(samples,numsubj);

%Pull gaze coordinates from all subjects into just xpos and ypos
for i=1:length(subjectfiles)
    gazedata = csvread(subjectfiles(i).name);
    xpos(:,i) = gazedata(1:length(xpos),1);   
    ypos(:,i) = gazedata(1:length(ypos),2);
    alltimes(:,i) = gazedata(1:length(alltimes),3);
end

xpos(xpos <= 1) = NaN;
xpos(xpos > 1024) = NaN;
ypos(ypos <= 1) = NaN;
ypos(ypos > 768) = NaN;

xpos = xpos .* .7031;
ypos = ypos .* .6328;

vheight = 486;
vwidth = 720;
npixels = vheight * vwidth;

xrandpool = xpos;
yrandpool = ypos;
xrandpool(isnan(xpos .* ypos)) = [];
yrandpool(isnan(xpos .* ypos)) = [];

files = {'ROI/HumanFace.csv','ROI/YellowMuppet.csv','ROI/PurpleMuppet.csv','ROI/GreenMuppet.csv','ROI/RedMuppet.csv'};
 
%Initialize arrays
   
h1 = csvread(files{1});
m1 = csvread(files{2});
m2 = csvread(files{3});
m3 = csvread(files{4});
m4 = csvread(files{5});
    
h1gaze = NaN(samples, numsubj);
m1gaze = h1gaze; m2gaze = h1gaze; m3gaze = h1gaze; m4gaze = h1gaze;
facegaze = h1gaze;

for j = 1:1:samples
        home
        disp(strcat('Sample:', num2str(j), '/', num2str(samples)))

        f = ceil(j/4); %have to use frames to index masks
        x=xpos(j,:);
        y=ypos(j,:);

        %Load saliency map
        load(strcat('saliency_itti_',channel,'/',num2str(f)))
        % get x and y coordinates of peak saliency value over the entire frame
        [py, px]=find(map==1);
        peaky(j)=py(1);
        peakx(j)=px(1);
    
        % rank saliency values
        salarray = reshape(map,npixels,1);
        ranks = tiedrank(salarray);
        salrankmat = reshape(ranks,vheight,vwidth);
         
        %define masks based on roi data
        h1mask = get_elliptical_mask(h1(f,1:4),[720 486]) .* h1(f,5);
        m1mask = get_elliptical_mask(m1(f,1:4),[720 486]) .* m1(f,5);
        m2mask = get_elliptical_mask(m2(f,1:4),[720 486]) .* m2(f,5);
        m3mask = get_elliptical_mask(m3(f,1:4),[720 486]) .* m3(f,5);
        m4mask = get_elliptical_mask(m4(f,1:4),[720 486]) .* m4(f,5);       
        facemask = h1mask | m1mask | m2mask | m3mask | m4mask;
        h1_mask_size = size(find(h1mask == 1),1);
        masks = {h1mask; m1mask; m2mask; m3mask; m4mask; facemask};
        for k = 1:length(masks)        
          mask = masks{k};
          mask_size = size(find(mask == 1),1);
          if mask_size == 0 
                sal_mean(j,k) = NaN; 
                salrank_mean(j,k) = NaN; 
                haspeak(j,k) = NaN;
            else
                sal_mean(j,k) = (sum(sum(map .* mask)))/mask_size;
                salrank_mean(j,k) = (sum(sum((salrankmat./npixels) .* mask)))/mask_size;
                haspeak(j,k) = mask(peaky(j),peakx(j));
          end;
        end
        
        %Choose random location for baseline roi comparison
        roi_orig = h1(f,1:4);
        half_x = ceil(abs(roi_orig(1)-roi_orig(3))/2);
        half_y = ceil(abs(roi_orig(2)-roi_orig(4))/2);
        new_xc = randi([0+half_x,720-half_x]);
        new_yc = randi([0+half_y,486-half_y]);
        roi_rand = [new_xc-half_x new_yc-half_y new_xc+half_x new_yc+half_y];
        
        randmask = get_elliptical_mask(roi_rand,[720 486]) .* h1(f,5);
        rand_mask_size = size(find(randmask == 1),1);  
        rand_mask_larger(j) = rand_mask_size-h1_mask_size;
        
        if rand_mask_size == 0 
                sal_mean_r(j) = NaN; 
                salrank_mean_r(j) = NaN; 
            else
                sal_mean_r(j) = (sum(sum(map .* randmask)))/rand_mask_size;
                salrank_mean_r(j) = (sum(sum((salrankmat./npixels) .* randmask)))/rand_mask_size;
          end;
          
        for i = 1:numsubj
            if isfinite(xpos(j,i)) && isfinite(ypos(j,i)) && xpos(j,i) > 0 && ypos(j,i)
                xgaze = round(x(i));
                ygaze = round(y(i));
                h1gaze(j,i) = h1mask(ygaze,xgaze);
                m1gaze(j,i) = m1mask(ygaze,xgaze);
                m2gaze(j,i) = m2mask(ygaze,xgaze);
                m3gaze(j,i) = m3mask(ygaze,xgaze);
                m4gaze(j,i) = m4mask(ygaze,xgaze);
                facegaze(j,i) = facemask(ygaze,xgaze);
                
                radius = 24;
                %Gaze saliency
                 mask=get_circle_mask([xgaze ygaze],radius,[720 486]);
                 mask_size = size(find(mask == 1),1);

                if mask_size == 0 
                    gazesal(j,i) = NaN; 
                    gazerank(j,i) = NaN; 
                    gazesal_r(j,i) = NaN; 
                    gazerank_r(j,i) = NaN; 
                else
                    gazesal(j,i) = (sum(sum(map .* mask)))/mask_size;
                    gazerank(j,i) = (sum(sum((salrankmat./npixels) .* mask)))/mask_size;

                    %Create gaze-inspired random point
                    randnum=randsample(length(xrandpool(:,i)),1);
                    randx = ceil(xrandpool(randnum, i));
                    randy = ceil(yrandpool(randnum, i));

                    randmask=get_circle_mask([randx randy],radius,[vwidth vheight]);
                    randmask_size = size(find(mask == 1),1);

                    gazesal_r(j,i) = (sum(sum(map .* randmask)))/randmask_size;
                    gazerank_r(j,i) = (sum(sum((salrankmat./npixels) .* randmask)))/randmask_size; 
                end;   
            end   
                          
        end
        
        if isnan(h1(f,5))
            h1gaze(j,:) = NaN;
        end
        if isnan(m1(f,5))
            m1gaze(j,:) = NaN;
        end
        if isnan(m2(f,5))
            m2gaze(j,:) = NaN;
        end
        if isnan(m3(f,5))
            m3gaze(j,:) = NaN;
        end
        if isnan(m4(f,5))
            m4gaze(j,:) = NaN;
        end
        
            
    if graphs == 1
        
        im = imread(strcat('ISC_saliency/frames/',num2str(f)),'JPG');
        imshow(im)
        hold on
        hImg = imshow(facemask); set(hImg, 'AlphaData', 0.5);
        for i = 1:6
            if not(isnan(x(i) .* y(i)))
                plot(x(i),y(i), 'g+','MarkerSize',7,'LineWidth',2);
            end
        end
        
        
        %saveas(gcf, strcat('saloverlay/',stim,'_',num2str(j),'.png'),'png');
        hold off
        pause(.033)
    end
end
%
%Summary stats
%
%Median split for each roi
roi_salrank_med = nanmedian(salrank_mean)
roi_salrank_mean = nanmean(salrank_mean)

%
threshold = 0:.001:1;
for t = 1:length(threshold)
    hits(t) = nanmean(sal_mean(not(isnan(sal_mean(:,1))),1) >= threshold(t));
    fas(t) = nanmean(sal_mean_r(not(isnan(sal_mean_r))) >= threshold(t));
end
plot(fas,hits)
AUC=-trapz(fas, hits);
%
%
roi_gaze = {h1gaze; m1gaze; m2gaze; m3gaze; m4gaze;facegaze};
for i = 1:numsubj
    if subjectfiles(i).name(8) == '6'
        group(i) = 1;
    elseif subjectfiles(i).name(8) == 'b'
        group(i) = 2;
    elseif subjectfiles(i).name(8) == '1'
        group(i) = 3;
    elseif subjectfiles(i).name(8) == '2'
        group(i) = 4;
    elseif subjectfiles(i).name(8) == 'a'
        group(i) = 5;
    elseif subjectfiles(i).name(8) == 'c'
        group(i) = 6;
    end
    %Individual AUC for predicting gaze
    threshold = 0:.001:1;
    for t = 1:length(threshold)
        hits_gaze(t,i) = nanmean(gazesal(not(isnan(gazesal(:,i))),i) >= threshold(t));
        fas_gaze(t,i) = nanmean(gazesal_r(not(isnan(gazesal_r(:,i))),i) >= threshold(t));
    end
    gazesal_mean(i) = nanmean(gazesal(:,i));
    gazerank_mean(i) = nanmean(gazerank(:,i));
    %plot(fas,hits)
    AUC_gaze(i)=-trapz(fas_gaze(:,i), hits_gaze(:,i));

    for k = 1:6
        kgaze = roi_gaze{k};
        salrank_quartiles = quantile(salrank_mean(not(isnan(salrank_mean(:,k))),k),3);
        
        sal_mean_att(i,k) = nanmean(sal_mean(kgaze(:,i) == 1,k));
        sal_mean_notatt(i,k) = nanmean(sal_mean(kgaze(:,i) == 0,k));
        sal_rank_att(i,k) = nanmean(salrank_mean(kgaze(:,i) == 1,k));
        sal_rank_notatt(i,k) = nanmean(salrank_mean(kgaze(:,i) == 0,k));
        fixate_whenpeak(i,k) = nanmean(kgaze(haspeak(:,k) == 1,i));
        fixate_notsal(i,k) = nanmean(kgaze(salrank_mean(:,k) <= salrank_quartiles(1),i));
        fixate_lesssal(i,k) = nanmean(kgaze(salrank_mean(:,k) < salrank_quartiles(2),i));
        fixate_moresal(i,k) = nanmean(kgaze(salrank_mean(:,k) > salrank_quartiles(2),i));
        fixate_sal(i,k) = nanmean(kgaze(salrank_mean(:,k) >= salrank_quartiles(3),i));
    end
end
id = 1:numsubj;
%
data = [id' group' sal_rank_att sal_rank_notatt fixate_whenpeak fixate_lesssal fixate_moresal fixate_notsal fixate_sal gazesal_mean' gazerank_mean' AUC_gaze'];
%mkdir('ROIGaze')
outfile = strcat('ROIGaze/');
csvwrite(strcat(outfile,channel,'.csv'),data);
save(strcat(outfile,channel))
data = m1gaze;