function [phiFor,phiForSave,contoursFor,contoursForSave,endFrame,outcome,statistics] = forwardsTracking(images,m,n,frameNr,numberOfFrames,phiBack,radius,mitosisThreshold,lambda1,lambda2,mu,nu,g_adjust_low,g_adjust_high,omega,timeStep,maxIterations,phiUpdate,epsilonNormGradReg,epsilonDeltaReg,Dxpk,Dypk,Dxmk,Dymk,Dxck,Dyck,gpu,err)
% Performance of the forwards tracking algorithm;
% initial boundary around mitotic cell taken from backwards tracking result

%% Preset energy terms

normVel_in = 0;    %Normal velocity term (inside)
normVel_out = 0;   %Normal velocity term (outside)
lengthReg = 0;     %Length regularisation
localStd = 0;      %Local standard deviation term
areaReg = 0;       %Area regularisation

%% Determine statistics for starting frame

f = frameNr+1;
if f>numberOfFrames || err
    phiFor = cell(1);
    phiForSave = cell(1);
    contoursFor = cell(1);
    contoursForSave = cell(1);
    endFrame = 0;
    outcome = 'unknown';
    statistics = struct([]);
    return
end
img = images(:,:,f);
phiFor(:,:,1) = phiBack(:,:,1);
phiFor(:,:,2) = phiFor(:,:,1);
phiForSave(:,:,1) = phiBack(:,:,1);

[area,~,~,~,~,~,~,~,~,~] = calculateStatistics(img,phiFor(:,:,1),epsilonNormGradReg);

crop_size = 100;
area_min = area/4;
timeStepInit = timeStep;

statistics = struct([]);

%% Tracking

c=2; %start counting

while f<=numberOfFrames && f<=frameNr+mitosisThreshold && ~exist('outcome','var') %stop as soon as outcome is known
    
    if exist('phiOld','var')
        clear phiOld
    end
    
    if f~=frameNr+1
        phiFor(:,:,c) = phiFor(:,:,c-1);
    end
    
    phi = phiFor(:,:,c) - 7;
    img = images(:,:,f);
    
    %Create ROI around cell interior (negative level set function) in order
    %to increase speed
    
    try
        massCentre = regionprops(phi<0, 'centroid');
        mC1 = round(massCentre.Centroid(1));
        mC2 = round(massCentre.Centroid(2));
    catch
        outcome = 'unknown';
        break
    end

    phi = phi(max(1,mC2-crop_size):min(m,mC2+crop_size),max(1,mC1-crop_size):min(n,mC1+crop_size)); %Crop phi
    img_roi = img(max(1,mC2-crop_size):min(m,mC2+crop_size),max(1,mC1-crop_size):min(n,mC1+crop_size)); %Crop img

    if lambda1~=0 || lambda2~=0
        normvel = normalVelocity(images,f,numberOfFrames); %Calculate normal velocity image
        normvel = normvel(max(1,mC2-crop_size):min(m,mC2+crop_size),max(1,mC1-crop_size):min(n,mC1+crop_size)); %Crop normvel
    end 
    
    if nu~=0
        if verLessThan('matlab','8.5')
            gaussfilt = fspecial('gaussian',5,1);
            smoothed = imfilter(img_roi,gaussfilt,'replicate','conv');
        else
            smoothed = imgaussfilt(img_roi,1,'FilterSize',5);
        end
        locstd = stdfilt(smoothed);
        scaled = (locstd - min(locstd(:))) ./ (max(locstd(:)) - min(locstd(:)));
        g = imadjust(scaled,[g_adjust_low g_adjust_high],[1 0],3); %Calculate g
    end

    if gpu
        phi = gpuArray(phi); 
        if lambda1~=0 || lambda2~=0
            normvel = gpuArray(normvel);
        end
        if nu~=0
            g = gpuArray(g);
        end
        Dxpk = gpuArray(Dxpk);
        Dypk = gpuArray(Dypk);
        Dxmk = gpuArray(Dxmk);
        Dymk = gpuArray(Dymk);
        Dxck = gpuArray(Dxck);
        Dyck = gpuArray(Dyck);
    end
    
    timeStep = timeStepInit;
    
    for i=1:maxIterations
        
        if lambda1~=0
            c_in = mean2(normvel(phi<0));
            normVel_in = (normvel - c_in).^2;
        end
        if lambda2~=0
            c_out = mean2(normvel(phi>0));
            normVel_out = (normvel - c_out).^2;
        end
        if mu~=0 || nu~=0
            phixp = conv2(phi,Dxpk,'same');
            phiyp = conv2(phi,Dypk,'same');
            phixc = conv2(phi,Dxck,'same');
            phiyc = conv2(phi,Dyck,'same');
            normGradxpyc = sqrt(phixp.^2 + phiyc.^2 + epsilonNormGradReg);
            normGradxcyp = sqrt(phixc.^2 + phiyp.^2 + epsilonNormGradReg);
            normGradc = sqrt(phixc.^2 + phiyc.^2 + epsilonNormGradReg);
        end
        if mu~=0
            lengthReg = conv2(phixp./normGradxpyc,Dxmk,'same') + conv2(phiyp./normGradxcyp,Dymk,'same');
        end
        if nu~=0
            localStd = normGradc .* (conv2(g.*phixc./normGradc,Dxck,'same') + conv2(g.*phiyc./normGradc,Dyck,'same'));
        end
        if omega~=0
            area_phi = bwarea(phi<0);
            if area_phi < area_min
                areaReg = area_phi - area_min;
            else
                areaReg = 0;
            end   
        end
        
        %Regularised Dirac delta function
        deltaPhi = 1/pi * epsilonDeltaReg ./ (epsilonDeltaReg^2 + phi.^2);
        
        %Gradient descent
        phitemp = phi + timeStep .* deltaPhi .* (lambda1 * normVel_in - lambda2 * normVel_out + mu * lengthReg + nu * localStd + omega * areaReg);
        
        %Topology preservation combined with narrow band method
        [xcoord,ycoord,numberofbandpixel]=createbandmapping2D(sign(phitemp.*phi)==-1);
        for bandPixel=1:numberofbandpixel
            x=xcoord(bandPixel);
            y=ycoord(bandPixel);
            phi(x,y)=topologyCheck2D(phi(x-1:x+1,y-1:y+1),phitemp(x,y),gpu);
        end

        %Reinitialise phi every phiUpdate iterations
        if (mod(i,phiUpdate)==0)
            phi = bwdist(phi<0) - bwdist(phi>=0);%signed distance function
        %Stop if level set function does not change significantly anymore
        if exist('phiOld','var')
            changePhi = abs(phi-phiOld);
            changePhi = sum(changePhi(:)) / numel(changePhi);
            if (changePhi < 10e-6)
                break;
            end
        end
        phiOld = phi;
        end
        
    end%loop through iterations
    
    if gpu
        phi = gather(phi);
    end
        
    %Resize phi  
    phiNew = ones(m,n);
    phiNew(max(1,mC2-crop_size):min(m,mC2+crop_size),max(1,mC1-crop_size):min(n,mC1+crop_size)) = phi;
    phi = bwdist(phiNew<0) - bwdist(phiNew>=0);
    phisave = phi;
    phiForSave(:,:,c) = phisave;
    
    %Calculate statistics for cell's interior (negative level-set function)
    [~,~,circularity,~,~,~,~,~,~,~] = calculateStatistics(img,phi,epsilonNormGradReg);
        
    %Create region of interest around the cell boundary
    bw = (phi<0); %black outside of contour and white inside
    stats = regionprops(bw,'BoundingBox'); %find box around white object
    x = round(stats.BoundingBox(1)); %x-coord upper-left corner of ROI
    y = round(stats.BoundingBox(2)); %y-coord upper-left corner of ROI
    x_width = stats.BoundingBox(3); %width of ROI along dimension x
    y_width = stats.BoundingBox(4); %width of ROI along dimension y
    %Make sure that ROI is not out of image bounds
    roi_extend = 20;
    if x>roi_extend
        x_begin = x-roi_extend;
    else
        x_begin = 1;
    end
    if y>roi_extend
        y_begin = y-roi_extend;
    else
        y_begin = 1;
    end
    if x+x_width>n-roi_extend
        x_end = n;
    else
        x_end = x+x_width+roi_extend;
    end
    if y+y_width>m-roi_extend
        y_end = m;
    else
        y_end = y+y_width+roi_extend;
    end 
    roi = img(y_begin:y_end,x_begin:x_end);
    
    %Search for circular objects within ROI
    [centres,radii,metrics] = imfindcircles(roi,[round(radius/2) round(radius)],'Sensitivity',0.95);

    %Make sure circles don't overlap too much or are too far from each other
    averageRadius = sum(radii)/length(radii);
    centreDistances = triu(squareform(pdist(centres)));
    [row,col] = find(centreDistances~=0&centreDistances<averageRadius&centreDistances>3*averageRadius);
    if ~isempty(row)
        deleteEntries = zeros(size(row));
        for i=1:length(row)
            if metrics(row(i)) < metrics(col(i))
                deleteEntries(i) = row(i);
            else
                deleteEntries(i) = col(i);
            end
        end
        deleteEntries = unique(deleteEntries);
        for i=length(deleteEntries):-1:1
            centres(deleteEntries(i),:)=[];
            radii(deleteEntries(i))=[];
            metrics(deleteEntries(i))=[];
        end
    end
                                                                       
    nrCells=size(radii,1);
    if nrCells==2
        outcome='regular';
    elseif nrCells==3
        outcome='3 daughter cells';
    elseif nrCells==0
        outcome='no division';
    elseif circularity<0.75
        outcome='apoptosis';
    end
    
    %Take circle as initial contour for level set function in next frame
    if nrCells==1
        t = 0:pi/4:2*pi; %Default: 9 points for circle coordinates interpolation
        s = size(t,2);
        xcoord = zeros(1,s);
        ycoord = zeros(1,s);
        for j=1:s  
            xcoord(j) = radii * cos(t(j)) + centres(1);
            ycoord(j) = radii * sin(t(j)) + centres(2);
        end
        pp = spline(0:1/(s-1):1 ,[xcoord; ycoord], 0:1/499:1);
        if ~(any(sum(find(round(pp(1,:))<1))) || any(sum(find(round(pp(1,:))>size(roi,2)))) || any(sum(find(round(pp(2,:))<1))) || any(sum(find(round(pp(2,:))>size(roi,1)))))
            phi = getLevelSetFromCoords(pp,size(roi,1),size(roi,2));
            %Resize phi  
            phiNew = ones(m,n);
            phiNew(y_begin:y_end,x_begin:x_end) = phi;
            phi = bwdist(phiNew<0) - bwdist(phiNew>=0);
        end
    elseif nrCells>1
        for k=1:nrCells
            t = 0:pi/4:2*pi; %Default: 9 points for circle coordinates interpolation
            s = size(t,2);
            xcoord = zeros(k,s);
            ycoord = zeros(k,s);
            for j=1:s  
                xcoord(k,j) = radii(k) * cos(t(j)) + centres(k,1);
                ycoord(k,j) = radii(k) * sin(t(j)) + centres(k,2);
            end
            pp(:,:,k) = spline(0:1/(s-1):1 ,[xcoord(k,:); ycoord(k,:)], 0:1/499:1);
            for i=1:500
                pp(1,i,k) = pp(1,i,k) + y_begin;
                pp(2,i,k) = pp(2,i,k) + x_begin;
            end
        end
    end
    
    phiFor(:,:,c) = phi;
    
    [area,perimeter,circularity,centroid,mean_hist,std_hist,tv_grey_normarea,normgrad_grey_max,meanlocstd,tv_locstd_normarea] = calculateStatistics(img,phisave,epsilonNormGradReg);
    statistics(c-1).area = area;
    statistics(c-1).perimeter = perimeter;
    statistics(c-1).circularity = circularity;
    statistics(c-1).centroid = centroid;
    statistics(c-1).mean_hist = mean_hist;
    statistics(c-1).std_hist = std_hist;
    statistics(c-1).tv_grey_normarea = tv_grey_normarea;
    statistics(c-1).normgrad_grey_max = normgrad_grey_max;
    statistics(c-1).meanlocstd = meanlocstd;
    statistics(c-1).tv_locstd_normarea = tv_locstd_normarea;    
        
    f = f+1; %frame number
    c = c+1; %counter
    
end%loop through frames

endFrame = f-1;

try
    contoursFor = cell(1,c-1);
    for i=1:c-1
        contoursFor(i) = bwboundaries((phiFor(:,:,i)<0),'noholes');
    end

    contoursForSave = cell(1,c-1);
    for i=1:c-1
        contoursForSave(i) = bwboundaries((phiForSave(:,:,i)<0),'noholes');
    end
catch
    contoursFor = cell(1,c-1);
    contoursForSave = cell(1,c-1);
    outcome = 'unknown';
end

if ~exist('outcome','var')
    outcome = 'unknown';
end

end
