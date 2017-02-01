function value = topologyCheck2D(windowPhi,tempValue,gpu)
%-> windowPhi = phi(x-1:x+1,y-1:y+1)
%-> tempValue = phitemp(x,y)
%-> x = xcoord(bandPixel), y = ycoord(bandPixel)

%Code adopted from "Sabine" framework by Michael Möller et al.
%Reference: M. Möller, M. Burger, P. Dieterich, A. Schwab. "A framework for
%automated cell tracking in phase contrast microscopic videos based on
%normal velocities." Journal of Visual Communication and Image
%Representation 25.2 (2014): 396-409.

%create windowTempPhi
%that is, create a 3x3 window with the same values as in windowPhi
%except the value in the centre which is set to tempValue
%(where we want to check whether it is a simple point or not)
windowTempPhi = windowPhi;
windowTempPhi(2,2) = tempValue;

%convert the windows to logicals
windowPhi = windowPhi < 0;
windowTempPhi = windowTempPhi < 0;

if gpu
    %check whether the change of sign also changed the number of object
    %connected components (8-connectivity)
    [~,T1] = bwlabel(windowPhi,8);
    %T1 (topological number): connected object components in windowPhi
    [~,T2] = bwlabel(windowTempPhi,8);
    %T2 (topological number): connected object components in windowTempPhi

    %check whether the change of sign also changed the number of background
    %connected components (4-connectivity)
    [~,T3] = bwlabel(~windowPhi,4);
    %T3 (topological number): connected background components in windowPhi
    [~,T4] = bwlabel(~windowTempPhi,4);
    %T4 (topological number): connected background components in windowTempPhi
else
    %check whether the change of sign also changed the number of object
    %connected components (8-connectivity)
    [~,T1] = bwlabeln(windowPhi,8);
    %T1 (topological number): connected object components in windowPhi
    [~,T2] = bwlabeln(windowTempPhi,8);
    %T2 (topological number): connected object components in windowTempPhi

    %check whether the change of sign also changed the number of background
    %connected components (4-connectivity)
    [~,T3] = bwlabeln(~windowPhi,4);
    %T3 (topological number): connected background components in windowPhi
    [~,T4] = bwlabeln(~windowTempPhi,4);
    %T4 (topological number): connected background components in windowTempPhi
end

if (T1~=T2 || T3~=T4 ) % NO SIMPLE POINT -> Topology would have changed
    value = -10e-2*sign(tempValue);
else % SIMPLE POINT
    value = tempValue;
end

end
