clear all %Read ROIs from file, turn into files for each ROI

xDoc = xmlread('ROI/ROIs.xml');

faces = xDoc.getElementsByTagName('DynamicAOI');

for j = 0:4

    kfs = faces.item(j).getElementsByTagName('KeyFrame');
    data = zeros(kfs.getLength,5);
    for i = 0:kfs.getLength - 1

       kf = kfs.item(i); 
       data(i+1,1) = str2num(kf.getElementsByTagName('Timestamp').item(0).getTextContent)/1000;
       pts = kf.getElementsByTagName('Point');
       data(i+1,2) = str2num(pts.item(0).getElementsByTagName('X').item(0).getTextContent);
       data(i+1,3) = str2num(pts.item(1).getElementsByTagName('X').item(0).getTextContent);
       data(i+1,4) = str2num(pts.item(0).getElementsByTagName('Y').item(0).getTextContent);
       data(i+1,5) = str2num(pts.item(1).getElementsByTagName('Y').item(0).getTextContent);
       data(i+1,6) = str2num(kf.getElementsByTagName('Visible').item(0).getTextContent);
    end
    
    name = faces.item(j).getElementsByTagName('Name');
    csvwrite(strcat('ROI/',char(name.item(0).getTextContent),'.csv'), data);
    
end

%%

clear all %Take each ROI and resample to correct time
data = csvread('ROI/MuppetFace.csv');

allpts = zeros(1800,5);
data(:,1) = (data(:,1) ./ 1000) .* 30;
for i = 0:length(data)
    
    if i == 0
        fstart = 1;
        fend = ceil(data(i+1,1)) - 1;
        for j = fstart:fend
            allpts(j,:) = [300 300 300 300 0];
        end
    elseif i > 0 && i < length(data)
        fstart = fend+1;
        fend = floor(data(i+1,1));
        x1s = linspace(data(i,2),data(i+1,2),numel(fstart:fend))';
        x2s = linspace(data(i,3),data(i+1,3),numel(fstart:fend))';
        y1s = linspace(data(i,4),data(i+1,4),numel(fstart:fend))';
        y2s = linspace(data(i,5),data(i+1,5),numel(fstart:fend))';
        vis = zeros(numel(fstart:fend),1);
        vis = vis + data(i,6);
        allpts(fstart:fend,:) = [x1s y1s x2s y2s vis];
    elseif i == length(data)
        fstart = fend + 1;
        fend = 1800;
        row = allpts(fstart-1,:);
        for j = fstart:fend
            allpts(j,:) = row;
            allpts(j,5) = 0;
        end
    end

end

csvwrite('ROI/PurpleMuppet.csv',allpts);

%INTERPOLATED FILES = Color-Muppet, Human Face
%%
%Visualize Interpolated ROI
clear all

allpts = csvread('ROI/HumanFace.csv');

for i = 1250:1500
    %d = ds(ds.frame == i);
    im = imread(strcat('frames/',num2str(i)),'JPG');
    imshow(im)
    hold on
    %plot(d.eyex, d.eyey, 'r+', 'MarkerSize', 10, 'MarkerWidth', 2);
     mask = get_elliptical_mask(allpts(i,1:4),[720 486]).*allpts(i,5);
     hImg = imshow(mask); set(hImg, 'AlphaData', 0.5);
    hold off
    %pause(.033)
    saveas(hImg,strcat('saloverlay/',num2str(i),'.jpg'),'jpg');
end
    
%% CALCULATE SALIENCE OF PIXELS WITHIN ROIS
clear all
for i = 1:1800
    
    vheight = 486;
    vwidth = 720;
    npixels = vheight * vwidth;
    
    % load saliency map
    load(strcat('saliency_itti/',num2str(i))) %map
    % rank saliency values
    salarray = reshape(map,npixels,1);
    ranks = tiedrank(salarray);
    salrankmat = reshape(ranks,vheight,vwidth);   
    
    files = {'ROI/HumanFace.csv','ROI/YellowMuppet.csv','ROI/PurpleMuppet.csv','ROI/GreenMuppet.csv','ROI/RedMuppet.csv'};
   for j = 1:5
        allpts = csvread(files{j}); 
        mask = get_elliptical_mask(allpts(i,1:4),[720 486]).*allpts(i,5);
        
        mask_size = size(find(mask == 1),1);
        if mask_size == 0 
            salcircle_mean(i,j) = NaN; 
            salcircleranks_mean(i,j) = NaN; 
            salcircle_max(i,j) = NaN; 
            salcircleranks_max(i,j) = NaN; 
        else
            salcircle_mean(i,j) = (sum(sum(map .* mask)))/mask_size;
            salcircleranks_mean(i,j) = (sum(sum((salrankmat./npixels) .* mask)))/mask_size;
            salcircle_max(i,j) = max(max(map .* mask));
            salcircleranks_max(i,j) = max(max((salrankmat./npixels) .* mask));
        end;
        
        if j == 1
            face = mask;
            muppet = face .* 0; %Set muppet mask to zeros
        elseif j > 1
            muppet = muppet | mask;
            face = face | mask;
        end
        background = not(face);
   end
   
    %Saliency of faces, muppet faces, background
    for j = 6:8
        if j == 6
            mask = face;
        elseif j == 7;
            mask = muppet;
        elseif j == 8;
            mask = background;
        end
        mask_size = size(find(mask == 1),1);
        if mask_size == 0 
            salcircle_mean(i,j) = NaN; 
            salcircleranks_mean(i,j) = NaN; 
            salcircle_max(i,j) = NaN; 
            salcircleranks_max(i,j) = NaN; 
        else
            salcircle_mean(i,j) = (sum(sum(map .* mask)))/mask_size;
            salcircleranks_mean(i,j) = (sum(sum((salrankmat./npixels) .* mask)))/mask_size;
            salcircle_max(i,j) = max(max(map .* mask));
            salcircleranks_max(i,j) = max(max((salrankmat./npixels) .* mask));
        end;
    end

   
end

csvwrite('ROI/c_mean.csv',salcircle_mean);
csvwrite('ROI/c_max.csv',salcircle_max);
csvwrite('ROI/cranks_mean.csv',salcircleranks_mean);
csvwrite('ROI/cranks_max.csv',salcircleranks_max);

%%
muppet_range = 1905:4608;
(length(muppet_range)/120)
    