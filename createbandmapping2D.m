function [xcoord,ycoord,numberofbandpixel] = createbandmapping2D(band)
%-> input "band" in tracking function is sign(phitemp.*phi)==-1
%-> that is, only potential simple points are taken into account

%Code adopted from "Sabine" framework by Michael Möller et al.
%Reference: M. Möller, M. Burger, P. Dieterich, A. Schwab. "A framework for
%automated cell tracking in phase contrast microscopic videos based on
%normal velocities." Journal of Visual Communication and Image
%Representation 25.2 (2014): 396-409.

[b1,b2] = size(band);    
%Keep the coordinates away from the boundary of the image
band(1,:) = 0;
band(b1,:) = 0;
band(:,1) = 0;
band(:,b2) = 0;
xcoord = [];
ycoord = [];
[xcoordTmp,ycoordTmp] = find(band>0);
xcoord = [xcoord;xcoordTmp];
ycoord = [ycoord;ycoordTmp];
numberofbandpixel = numel(xcoord);

end
