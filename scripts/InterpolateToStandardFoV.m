function [ StdImage ] = InterpolateToStandardFoV( Image,SAG,COR,TRA,FOVST,NST,MaxSeparation )
% INTERPOLATETOSTANDARDFOV Summary of this function goes here
% Detailed explanation goes here

Dim = size(squeeze(Image));
if (length(Dim)<3)
    error('3D Data expected')
end

% Sum of squares of higher dimensions
if(length(Dim)>4)
    for cDim=5:length(Dim)
        Image = sum(abs(Image).^2,cDim);
    end
    Image = sqrt(Image);
    warning('Summing dimensions that are higher than 4')
end


StdImage = InterpolateToStandardFoVMex( double(Image),...
                                        double(SAG),double(COR),double(TRA),...
                                        double(FOVST),...
                                        double(NST),...
                                        double(MaxSeparation)); 
end

