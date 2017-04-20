function mask = get_ellipse(coord, image_size)
% Get circle mask in a specific potion with specific radius.
% mask = get_circle_mask(center, R);
%     center: the center position of the circle.
%             center(1) is horizontal coordinate.
%             center(2) is vertical coordinate.
%     R: radius of the circle
%     image_size: optional. Default is [720 480];
%             image_size(1) is horizontal size.
%             image_size(2) is vertical size.
% Example: mask =
%     mask = get_circle_mask([90 75], 20);
%     imwrite(mask, file_name);
%
if ~exist('image_size','var')
    image_size = [720 480];
end

x1 = round(coord(1));
y1 = round(coord(2));
x2 = round(coord(3));
y2 = round(coord(4));

mask = zeros(image_size(2), image_size(1));  % hard coding for now

% if any(center <= R) || any(center > (image_size - R))
%     warning('The specified circle exceeds the border of the image.');
%     return;
% end


%Define X and Y axes
Rx = round(abs(x1 - x2)/2);
Ry = round(abs(y1 - y2)/2);


%Define center
X = x1 + Rx; Y = y1 + Ry;

for i = Y-Ry:Y+Ry
    for j= X-Rx:X+Rx
        if ((i-Y)^2)/Ry^2 + ((j-X)^2)/Rx^2  == 1
            if i > 0 && j > 0 && i < image_size(2) && j < image_size(1)
                mask(i,j) = 1;
            end
        end
    end
end

area = (pi * Rx * Ry);

mask = logical(mask);
end