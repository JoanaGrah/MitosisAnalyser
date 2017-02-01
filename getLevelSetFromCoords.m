function levelSet = getLevelSetFromCoords(coordinates,m,n)
%Get Level Set Function (Signed Distance Function) from Coordinates

%Code adopted from "Sabine" framework by Michael Möller et al.
%Reference: M. Möller, M. Burger, P. Dieterich, A. Schwab. "A framework for
%automated cell tracking in phase contrast microscopic videos based on
%normal velocities." Journal of Visual Communication and Image
%Representation 25.2 (2014): 396-409.

roundpp = round(coordinates);
temp = ones(m,n);
for k=1:size(roundpp,2)
    if roundpp(2,k)>0 && roundpp(2,k)<=n && roundpp(1,k)>0 && roundpp(1,k)<=m
        temp(roundpp(2,k),roundpp(1,k)) = 0;
    end
end
se = strel('disk',2); %circular structuring element for erosion
temp = imerode((temp==1),se); %erosion
classitemp = bwlabel(temp);
if sum(sum((classitemp==1)))>sum(sum((classitemp==2))) %smaller region is supposed to be the cell
    temp = (classitemp==2);
else
    temp = (classitemp==1);
end
levelSet = -bwdist(1-temp) + bwdist(temp); %signed distance function

end