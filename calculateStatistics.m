function [area,perimeter,circularity,centroid,mean_hist,std_hist,tv_grey_normarea,normgrad_grey_max,meanlocstd,tv_locstd_normarea] = calculateStatistics(img,levelSet,epsilon_grad)
% Calculation of image statistics (cell interior / negative level-set function)
%  1) Area = Number of pixels inside of the cell
%  2) Perimeter = Length of zero-level-set (number of contour pixels)
%  3) Circularity
%  4) Centroid
%  5) Mean of histogram
%  6) Standard deviation of histogram
%  7) Total variation of grey values normalised by area
%  8) Maximum of norm of gradient of grey values
%  9) Mean of local standard deviation
% 10) Total variation of local standard deviation normalised by area

rp = regionprops((levelSet<0),'Area','Perimeter','Centroid');

area = rp.Area; %1)

perimeter = rp.Perimeter; %2)

circularity = 4*pi*area/perimeter^2; %3)

centroid = rp.Centroid; %4)

mean_hist = mean(img((levelSet<0))); %5)

std_hist = std(double(img((levelSet<0)))); %6)

% From here on, we take the image convolved with a Gaussian kernel into account
if verLessThan('matlab','8.5')
    gaussfilt = fspecial('gaussian',5,1);
    convolution = imfilter(img,gaussfilt,'replicate','conv');
else
    convolution = imgaussfilt(img,1,'FilterSize',5);
end
    
[gradx,grady] = gradient(convolution);
normGrad2 = gradx.^2 + grady.^2 + epsilon_grad;
normgrad_grey = sqrt(normGrad2);

locstd = stdfilt(convolution);
[gradx2,grady2] = gradient(locstd);
normgrad_locstd = sqrt(gradx2.^2 + grady2.^2 + epsilon_grad);

tv_grey_normarea = sum(sum(normgrad_grey((levelSet<0))))/area; %7)

normgrad_grey_max = max(max(normgrad_grey((levelSet<0)))); %8)

meanlocstd = mean2(locstd((levelSet<0))); %9)

tv_locstd_normarea = sum(sum(normgrad_locstd((levelSet<0))))/area; %10)

end