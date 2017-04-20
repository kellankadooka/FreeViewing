function [data] = TestGazeROIs(stim, duration, graphs)
%
% graphs = 0;
% stim = 'sesameAd';
% duration = 59;

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

files = {'ROI/HumanFace.csv','ROI/YellowMuppet.csv','ROI/PurpleMuppet.csv','ROI/GreenMuppet.csv','ROI/RedMuppet.csv'};
 
%Initialize arrays
   
h1 = csvread(files{1});
m1 = csvread(files{2});
m2 = csvread(files{3});
m3 = csvread(files{4});
m4 = csvread(files{5});
    
h1gaze = NaN(samples, numsubj);
m1gaze = h1gaze; m2gaze = h1gaze; m3gaze = h1gaze; m4gaze = h1gaze;

for j = 1:1:samples
        home
        disp(strcat('Sample:', num2str(j), '/', num2str(samples)))

        f = ceil(j/4); %have to use frames to index masks
        x=xpos(j,:);
        y=ypos(j,:);
        
        h1mask = get_elliptical_mask(h1(f,1:4),[720 486]) .* h1(f,5);
        m1mask = get_elliptical_mask(m1(f,1:4),[720 486]) .* m1(f,5);
        m2mask = get_elliptical_mask(m2(f,1:4),[720 486]) .* m2(f,5);
        m3mask = get_elliptical_mask(m3(f,1:4),[720 486]) .* m3(f,5);
        m4mask = get_elliptical_mask(m4(f,1:4),[720 486]) .* m4(f,5);
        for i = 1:6
            if isfinite(xpos(j,i)) && isfinite(ypos(j,i)) && xpos(j,i) > 0 && ypos(j,i)
                xgaze = round(x(i));
                ygaze = round(y(i));
                h1gaze(j,i) = h1mask(ygaze,xgaze);
                m1gaze(j,i) = m1mask(ygaze,xgaze);
                m2gaze(j,i) = m2mask(ygaze,xgaze);
                m3gaze(j,i) = m3mask(ygaze,xgaze);
                m4gaze(j,i) = m4mask(ygaze,xgaze);
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
        hImg = imshow(h1mask); set(hImg, 'AlphaData', 0.5);
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

outfile = strcat('ROI/',stim,'_');
csvwrite(strcat(outfile,'h1.csv'),h1gaze);
csvwrite(strcat(outfile,'m1.csv'),m1gaze);
csvwrite(strcat(outfile,'m2.csv'),m2gaze);
csvwrite(strcat(outfile,'m3.csv'),m3gaze);
csvwrite(strcat(outfile,'m4.csv'),m4gaze);
data = m1gaze;
