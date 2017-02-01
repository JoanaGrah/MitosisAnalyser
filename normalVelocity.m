function normvel = normalVelocity(images,currentFrame,numberOfFrames)
%This function calculates the normal velocity for the current frame in a
%given image sequence based on the previous and subsequent image, if
%existent.

%normal velocity = |d/dt f(x,t)|/|grad(f(x,t))|,
%where the time derivative is discretised with central differences and the
%spatial derivative is discretised with backwards differences using an
%epsilon regularisation for the gradient

%Code adopted from "Sabine" framework by Michael Möller et al.
%Reference: M. Möller, M. Burger, P. Dieterich, A. Schwab. "A framework for
%automated cell tracking in phase contrast microscopic videos based on
%normal velocities." Journal of Visual Communication and Image
%Representation 25.2 (2014): 396-409.

epsilon = 0.001;

if currentFrame==1
    timederivative = abs(images(:,:,currentFrame) - images(:,:,currentFrame+1));
elseif currentFrame==numberOfFrames
    timederivative = abs(images(:,:,currentFrame-1) - images(:,:,currentFrame));
else
    timederivative = abs(images(:,:,currentFrame-1) - images(:,:,currentFrame+1));
end

gradval = sqrt(imfilter(images(:,:,currentFrame),[-1 1 0],'same').^2 + imfilter(images(:,:,currentFrame),[-1 1 0]','same').^2 + epsilon);
normalvelocity = sum(timederivative./gradval,3);
normvel = normalvelocity .* (normalvelocity<2) + (normalvelocity>2);

end