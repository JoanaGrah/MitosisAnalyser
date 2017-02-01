function eventData = mitosisDetection(images,numberOfFrames,mitosisThreshold,radiusMin,radiusMax,sensitivity)

wb = waitbar(0, 'MITOSIS DETECTION: Please wait...');
set(findobj(wb,'type','patch'),'edgecolor','k','facecolor','b');

disp('Starting Mitosis Detection...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function automatically detects mitotic cells in a sequence of     %
%  images. It uses the built-in function "imfindcircles" which is based   %
%  on the circular Hough transform in order to detect the typically very  %
%  round cells immediately before undergoing mitosis. The centres and     %
%  boundary coordinates are saved in the "eventData" structure so that    %
%  the circles surrounding the cells can be used as an initial contour    %
%  for tracking algorithms.                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% images:               grey scale image sequence to be analysed                                                                                                     %
% numberOfFrames:       number of frames in given image sequence                                                                                                     %
% mitosisThreshold:     What is the average mitosis length? / Are there already existing mitotic cells in the previous mitosisThreshold frames?                      %
% radiusMin, radiusMax: only find circles with radii within the interval [radiusMin radiusMax]                                                                       %
% sensitivity:          the higher the sensitivity, the more circular objects are detected / (1-sensitivity) corresponds to threshold of values in accumulator array %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eventData: structure array with fields "frames", "centres", "radii", "metrics" and "boundaries" %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','images:imfindcircles:warnForSmallRadius');
warning('off','images:imfindcircles:warnForLargeRadiusRange');

t = 0:pi/4:2*pi; %Default: 9 points for circle coordinates interpolation
s = size(t,2);

%Preallocation
eventData = struct([]);
eventsPerFrame = zeros(numberOfFrames,1);
eventFrames = [];
centres = [];
radii = [];
metrics = [];

%Search for circles in all frames
for k=1:numberOfFrames
    [centres_new,radii_new,metrics_new] = imfindcircles(images(:,:,k),[radiusMin radiusMax],'Sensitivity',sensitivity);
    if ~isempty(radii_new)
        eventFrames = [eventFrames;k];
        centres_new(:,3) = k;
        radii_new(:,2) = k;
        metrics_new(:,2) = k;
        centres = [centres;centres_new];
        radii = [radii;radii_new];
        metrics = [metrics;metrics_new];
    end
    eventsPerFrame(k) = size(radii_new,1);
    waitbar(k/(numberOfFrames+1));
    drawnow;
end
clear centres_new radii_new metrics_new

%Sort results by metrics and delete circles surrounding identical cells

largestRadius = max(radii(:,1));

allInfo = table(metrics(:,1),radii(:,1),centres(:,1:2),centres(:,3),'VariableNames',{'metrics','radii','centres','frames'});
sortedMetrics = sortrows(allInfo,1,'descend');
thr = floor(mitosisThreshold/2);
k = 1;
while k<=size(sortedMetrics,1)
    result(k,:) = sortedMetrics(k,:);
    curFrame = sortedMetrics.frames(k);
    indices = find(allInfo.frames>=max(1,curFrame-thr) & allInfo.frames<=min(numberOfFrames,curFrame+thr));
    considered = allInfo(indices(1):indices(end),:);
    consideredCentres = considered.centres;
    sizeTemp = size(consideredCentres,1);
    distances = zeros(sizeTemp,1);
    for l=1:sizeTemp
        distances(l) = sqrt((sortedMetrics.centres(k,1)-consideredCentres(l,1)).^2 + (sortedMetrics.centres(k,2)-consideredCentres(l,2)).^2);
    end
    deleteEntries = find(distances~=0 & distances<largestRadius);
    sizeTemp = length(deleteEntries);
    for l=sizeTemp:-1:1
        allInfo(indices(1)+deleteEntries(l)-1,:) = [];
    end
    sortedMetrics = sortrows(allInfo,1,'descend');
    k = k+1;
end

nrEvents = size(result,1);
sortedResult = table2array(sortrows(result,4));

for i = 1:nrEvents
    eventData(i).metrics = sortedResult(i,1);
    eventData(i).radii = sortedResult(i,2);
    eventData(i).centres = sortedResult(i,3:4);
    eventData(i).frames = sortedResult(i,5);
    xcoord = zeros(1,s);
    ycoord = zeros(1,s);
    for j=1:s  
        xcoord(j) = eventData(i).radii * cos(t(j)) + eventData(i).centres(1);
        ycoord(j) = eventData(i).radii * sin(t(j)) + eventData(i).centres(2);
    end
    eventData(i).boundaries = spline(0:1/(s-1):1 ,[xcoord; ycoord], 0:1/499:1);
end

waitbar(1);
drawnow;

warning('on','images:imfindcircles:warnForSmallRadius');
warning('on','images:imfindcircles:warnForLargeRadiusRange');

disp('Mitosis detection done.')

close(wb);

end