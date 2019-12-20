function [filtered_k_space] = hannFilt(k_space, mode,factor)
% function [filtered_k_space] = hannFilt(k_space)
%  
% Filters the k-space with a Hanning filter to reduce ringing effects.
% Supports 3D datasets
%
%   input:
%       - k_space: the unfiltered k-space
%       -    mode: 'rad' or 'sep'
%
%   output:
%       - filtered_k_space: the filtered k-space
%

    if ndims(k_space) > 3
        error('Sorry only datasets up to 3D supported!');
    end
    
    if nargin < 2
        mode = 'rad';
    end
    
    if nargin < 3
        factor = 1;
    end
    
    sz = [0; 0; 0];
    sz(1:ndims(k_space)) = size(k_space);
    
    kx = linspace(-pi,pi,sz(2))*factor;
    ky = linspace(-pi,pi,sz(1))*factor;
    if sz(3) ~= 0
        kz = linspace(-pi,pi,sz(3))*factor;
        [KX,KY,KZ] = meshgrid(kx,ky,kz);
    else
         kz = 0;
        [KX,KY] = meshgrid(kx,ky);
        KZ=0;
    end
    
    if strcmpi(mode,'rad')
        R = sqrt(KX.^2+KY.^2+KZ.^2);
        X   = R>pi;
        R(X) = pi;
        filtered_k_space = k_space.*my_hann(R);
    elseif strcmpi(mode,'sep')
        
        filtered_k_space = k_space ...
                           .* my_hann(KX) .* my_hann(KY) .* my_hann(KZ);
    else
        error('Sorry filter mode not supported!');
    end
        
end

function omega = my_hann(n)
    n(abs(n)>pi)=pi;
    omega=0.5.*(cos(n)+1);
    %+0.5warning('Hannfilter changed+0.5 ')
end