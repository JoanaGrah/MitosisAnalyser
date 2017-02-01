function [phiBack,contoursBack,beginFrame,statistics,err] = backwardsTracking(images,m,n,frameNr,numberOfFrames,boundaries,mitosisThreshold,lambda1,lambda2,mu,nu,g_adjust_low,g_adjust_high,omega,timeStep,maxIterations,phiUpdate,epsilonNormGradReg,epsilonDeltaReg,Dxpk,Dypk,Dxmk,Dymk,Dxck,Dyck,gpu)
% Performance of the backwards tracking algorithm starting from the circle
% around the detected mitotic cell

err = 0;

%% Initialisation of Level Set Function

initPhi = getLevelSetFromCoords(boundaries,m,n);

%% Preset energy terms

normVel_in = 0;    %Normal velocity term (inside)
normVel_out = 0;   %Normal velocity term (outside)
lengthReg = 0;     %Length regularisation
localStd = 0;      %Local standard deviation term
areaReg = 0;       %Area regularisation

%% Determine statistics for starting frame

f = frameNr;
img = images(:,:,f);
phiBack(:,:,1) = initPhi - 5;

[area,~,circularity,~,~,~,~,~,~,~] = calculateStatistics(img,phiBack(:,:,1),epsilonNormGradReg);

area_min = area/2;
timeStepInit = timeStep;
crop_size = 100;
thr_circularity = 0.9;
contourExtension = 10;

statistics = struct([]);

%% Tracking

c = 1;%start counting

while c==1 || (f>0 && f>frameNr-mitosisThreshold && circularity>thr_circularity) %stop as soon as cell reaches regular flat state
    
    if exist('phiOld','var')
        clear phiOld
    end
    
    %Reinitialise level set function by slightly extending the previous one
    if f~=frameNr     
        phiBack(:,:,c) = phiBack(:,:,c-1) - contourExtension;
    end
    
    phi = phiBack(:,:,c);
    img = images(:,:,f);
    
    %Create ROI around cell interior (negative level set function) in order
    %to increase speed
    
    massCentre = regionprops(phi<0, 'centroid');
    mC1 = round(massCentre.Centroid(1));
    mC2 = round(massCentre.Centroid(2));

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
            if gpu
                normGradxpyc = bsxfun(@hypot,phixp,phiyc)+epsilonNormGradReg;
                normGradxcyp = bsxfun(@hypot,phixc,phiyp)+epsilonNormGradReg;
                normGradc = bsxfun(@hypot,phixc,phiyc)+epsilonNormGradReg;
            else
                normGradxpyc = sqrt(phixp.^2 + phiyc.^2 + epsilonNormGradReg);
                normGradxcyp = sqrt(phixc.^2 + phiyp.^2 + epsilonNormGradReg);
                normGradc = sqrt(phixc.^2 + phiyc.^2 + epsilonNormGradReg);
            end
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
    
    %Don't get too close to image boundary
    tmp = phi<contourExtension;
    [m1,n1] = size(phi);
    boundPix = vertcat(tmp(1:m1,1),tmp(1,1:n1)',tmp(m1,1:n1)',tmp(1:m1,n1));
    if any(boundPix)
        err = 1;
        break
    end
    
    if gpu
        phi = gather(phi);
    end
    
    %Resize phi  
    phiNew = ones(m,n);
    phiNew(max(1,mC2-crop_size):min(m,mC2+crop_size),max(1,mC1-crop_size):min(n,mC1+crop_size)) = phi;
    phi = bwdist(phiNew<0) - bwdist(phiNew>=0);
    
    %Calculate statistics for cell's interior (negative level-set function)
    [area,perimeter,circularity,centroid,mean_hist,std_hist,tv_grey_normarea,normgrad_grey_max,meanlocstd,tv_locstd_normarea] = calculateStatistics(img,phi,epsilonNormGradReg);
    
    phiBack(:,:,c) = phi;
    statistics(c).area = area;
    statistics(c).perimeter = perimeter;
    statistics(c).circularity = circularity;
    statistics(c).centroid = centroid;
    statistics(c).mean_hist = mean_hist;
    statistics(c).std_hist = std_hist;
    statistics(c).tv_grey_normarea = tv_grey_normarea;
    statistics(c).normgrad_grey_max = normgrad_grey_max;
    statistics(c).meanlocstd = meanlocstd;
    statistics(c).tv_locstd_normarea = tv_locstd_normarea;    
   
    f = f-1;%frame number
    c = c+1;%counter
    
end%loop through frames

beginFrame = f+1;

try
    contoursBack = cell(1,c-1);
    for i=1:c-1
        result = bwboundaries((phiBack(:,:,i)<0),'noholes');
        contoursBack(i) = result(1);
    end
catch
    contoursBack = cell(1,c-1);
    phiBack(:,:,c-1) = zeros(m,n);
    err = 1;
end

end
