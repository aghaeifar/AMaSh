function [recon,wfull,cmap]=openadapt_3D(im,st,donorm,rn)
%   Adaptive recon based on Walsh et al.
%   Walsh DO, Gmitro AF, Marcellin MW.
%   Adaptive reconstruction of phased array MR imagery. 
%   Magn Reson Med. 2000 May;43(5):682-90.
%
%    and
%
%   Mark Griswold, David Walsh, Robin Heidemann, Axel Haase, Peter Jakob. 
%   The Use of an Adaptive Reconstruction for Array Coil Sensitivity Mapping and Intensity Normalization, 
%   Proceedings of the Tenth  Scientific Meeting of the International Society for Magnetic Resonance in Medicine pg 2410 (2002)
%
%   IN:         im                 image to be reconstructed          (#coils, Ny, Nx)
%               donorm             flag determining whether to normalize image intensity    
%               rn                 input noise covariance matrix.     [#coils, #coils)  
%   
%   OUT:        recon              Reconstructed image                (Ny, Nx)    
%               cmap               "Coil maps"                        (# coils,  Ny, Nx) 
%
%   This non-optimized function will calculate adaptively estimated coil maps
%   based on the Walsh algorithm for use in either optimal combination of 
%   array images, or for parallel imaging applications. The donorm flag can be
%   used to include the normalization described in the abstract above. This is
%   only optimal for birdcage type arrays, but will work reasonably for many other geometries.
%   This normalization is only applied to the recon'd image, not the coil maps.
%
%   The rn matrix should be the noise covariances as decribed in Walsh et al.
%   This could also be a covariance of a region of distrubing artifact as 
%   described by Walsh et al.
%
%   The default block size is 4x4. One can also use interpolation to speed
%   up the calculation (see code), but this could cause some phase errors in practice.
%   Just pay attention to what you are doing.
%
%   Please read the license text at the bottom of this program. By using this program, you 
%   implicity agree with the license. 
%
%   The main points of the license:
%
%   1) This code is strictly for non-commercial applications. The code is protected by
%      multiple patents.
%   2) This code is strictly for research purposes, and should not be used in any
%      diagnostic setting.
%

%   10/1/2001  Mark Griswold
%   22/10/2004  MG - Updated help and changed interpolation to be more
%   stable.

[nc,ny,nx,nz] = size(im);

[~,maxcoil]   = max(sum(abs(im(:,:)),2));  %find coil with maximum intensity
                                           %for correcting the phase of all 
                                           %of the other coils.    
if ~exist('donorm','var')
    donorm=0;
end
if ~exist('rn','var')
    rn     = eye(nc);
end

% in 3D können wir blocksize reduzieren, da generell mehr pixel fuer 
% statistik zur verfuegung stehen (b_y*b_x*bs_z statt bs_y*bs_x)
if nz > 1
    bs_y = min(8,2*floor(ny/2));  %y-block size
    bs_x = min(8,2*floor(nx/2));  %x-block size
   % bs_z = min(8,2*floor(nz/2));  %x-block size
    bs_z = min(min(4,2*ceil(nz/4)),2*floor(nz/2)+1);  %z-block size
%     bs_y = min(16,2*floor(ny/2));  %y-block size
else
    bs_y = min(8,2*floor(ny/2));  %y-block size
    bs_x = min(8,2*floor(nx/2));  %x-block size
    bs_z =         1;             %z-block size
end

%increase to set interpolation step size
st_y = min(2,ny);
st_x = min(2,nx);
st_z = min(1,nz);
if exist('st','var')
    st_y = st(1);
    if numel(st)==1
        st_x = st(1);
        st_z = st(1);
    end
    if numel(st)>=2
        st_y = st(2);
    end
    if numel(st)>=3
        st_z = st(3);
    end
end

nysmall = round(ny/st_y);
nxsmall = round(nx/st_x);
nzsmall = round(nz/st_z);

small_sz  = [nc,nysmall,nxsmall,nzsmall];

% vectorize cmapsmall for now
cmapsmall = zeros(small_sz(1),prod(small_sz(2:end)), class(im));

%Collect block for calculation of blockwise values. 
%Edges are cropped so the results
%near the edges of the image could 
%be in error. Not normally a problem.
%But watch out for aliased regions.
%precalculate index ranges for faster calculation
ymin = zeros(1,nysmall); ymax = ymin; ly = ymin;
xmin = zeros(1,nxsmall); xmax = xmin; lx = xmin;
for y=1:nysmall
    ymin(y) = max(st_y*y-floor(bs_y/2), 1);
    ymax(y) = min(st_y*y+ceil(bs_y/2), ny);
    ly(y)   = ymax(y)-ymin(y)+1;
end
for x=1:nxsmall
    xmin(x) = max(st_x*x-floor(bs_x/2), 1);
    xmax(x) = min(st_x*x+ceil(bs_x/2), nx);
    lx(x)   = xmax(x)-xmin(x)+1;
end
    
inv_rn = rn\eye(nc); % =inv(rn)

cnt = 0;
for z=1:nzsmall
    zmin = max(st_z*z-floor(bs_z/2), 1);
    zmax = min(st_z*z+ceil(bs_z/2), nz);
    lz   = zmax-zmin+1;
    for x=1:nxsmall
        for y=1:nysmall
            m1 = reshape(im(:,ymin(y):ymax(y),xmin(x):xmax(x),zmin:zmax),nc,ly(y)*lx(x)*lz);

            m = m1*m1';                                %Calculate signal covariance
            
            [v,d]   = eig(inv_rn*m);                   %Eigenvector with max eigenval gives
            [~,ind] = max(diag(d));                    %the correct combination coeffs.
            
            cnt = cnt+1;
            cmapsmall(:,cnt) = v(:,ind);
%             cmapsmall(:,cnt) = conj(v(:,ind));
        end
    end
end

%Correct phase based on coil with max intensity
cmapsmall = bsxfun(@times,cmapsmall,exp(-1j*angle(cmapsmall(maxcoil,:))));

% calc weights for coil combination
wsmall = cmapsmall;
if ~(donorm || nargout > 2)
    % we don't need cmapsmall anymore...
    clear cmapsmall;
end
for k=1:size(wsmall,2)
    tmp = wsmall(:,k);
    wsmall(:,k) = tmp/(tmp'*inv_rn*tmp);
end


if 1%(st_x*st_y*st_z ~= 1)
    %Now have to interpolate these weights up to the full resolution. This is done separately for
    %magnitude and phase in order to avoid 0 magnitude pixels between +1 and -1 pixels.
    [y,x,z] = ndgrid(1:(nysmall-1)/(ny-1):nysmall,1:(nxsmall-1)/(nx-1):nxsmall,1:(nzsmall-1)/(nz-1):nzsmall);
    wfull   = zeros(nc,ny,nx,nz,class(im));
    wsmall  = reshape(wsmall,small_sz);
    % permute wsmall/cmapsmall for faster access
    wsmall  = permute(wsmall,[2 3 4 1]);
    if (donorm || nargout > 2)
        cmap      = zeros(nc,ny,nx,nz,class(im));
        cmapsmall = reshape(cmapsmall,small_sz);
        cmapsmall = permute(cmapsmall,[2 3 4 1]);
    end
    for c=1:nc
        if nz > 1
            wfull(c,:,:,:)    = interpn(abs(   wsmall(:,:,:,c)),y,x,z,'linear') .* exp(-1j.*interpn(angle(   wsmall(:,:,:,c)),y,x,z,'nearest'));
        else % for 2D scans
            wfull(c,:,:,:)    = interpn(abs(   wsmall(:,:,:,c)),y,x,'linear') .* exp(-1j.*interpn(angle(   wsmall(:,:,:,c)),y,x,'nearest'));
        end
        if (donorm || nargout > 2)
            if nz > 1
                cmap(c,:,:,:) = interpn(abs(cmapsmall(:,:,:,c)),y,x,z,'linear') .* exp( 1j.*interpn(angle(cmapsmall(:,:,:,c)),y,x,z,'nearest'));
            else % for 2D scans
                cmap(c,:,:,:) = interpn(abs(cmapsmall(:,:,:,c)),y,x,'linear') .* exp( 1j.*interpn(angle(cmapsmall(:,:,:,c)),y,x,'nearest'));
            end
        end
    end
    clear y x z;

else
    wsmall = reshape(wsmall,small_sz);
    wfull  = wsmall;
    if (donorm || nargout > 2)
        cmapsmall = reshape(cmapsmall,small_sz);
        cmap      = cmapsmall;
    end
end
clear wsmall cmapsmall;

recon = squeeze(sum(wfull.*im));   %Combine coil signals. 

if donorm
    recon = recon.*squeeze(sum(abs(cmap))).^2;   %This is the normalization proposed in the abstract 
                                                 %referenced in the header.
end

% You should carefully read the following terms and conditions before installing or using the 
% software. Unless you have entered into a separate written license agreement with 
% Universit�t W�rzburg providing otherwise, installation or use of the software indicates your 
% agreement to be bound by these terms and conditions. 
% 
% Use of the software provided with this agreement constitutes your acceptance of these terms. 
% If you do NOT agree to the terms of this agreement, promptly remove the software together 
% with all copies from your computer. User's use of this software is conditioned upon compliance 
% by user with the terms of this agreement. 
% 
% Upon ordering, downloading, copying, installing or unencrypting any version of the software, you
% are reaffirming that you agree to be bound by the terms of this agreement. 
% 
% License to use 
% 
% Universit�t W�rzburg grants to you a limited, non-exclusive, non-transferable and non-assignable 
% license to install and use this software for research purposes. Use of this software for any 
% diagnostic imaging procedure is strictly forbidden.
% 
% License to distribute 
% 
% Please feel free to offer the non-commercial version of this software on any website, CD, or 
% bulletin board, demonstrate the non-commercial version of the software and its capabilities, or 
% give copies of the non-commercial version of the software to other potential users, so that others 
% may have the opportunity to obtain a copy for use in accordance with the license terms contained
% here. 
% 
% You agree you will only copy the non-commercial version of the software in whole with this 
% license and all delivered files, but not in part. 
% 
% Termination 
% 
% This license is effective until terminated. You may terminate it at any point by destroying 
% the software together with all copies of the software. 
% 
% If you have acquired a non-commercial version, the license granted herein shall automatically 
% terminate if you fail to comply with any term or condition of this Agreement. 
% 
% Also, Universit�t W�rzburg has the option to terminate any license granted herein if you fail 
% to comply with any term or condition of this Agreement. 
% 
% You agree upon such termination to destroy the software together with all copies of the software.
% 
% 
% Copyright 
% 
% The software is protected by copyright law. You acknowledge that no title to the intellectual 
% property in the software is transferred to you. You further acknowledge that title and full 
% ownership rights to the software will remain the exclusive property of Universit�t W�rzburg, 
% and you will not acquire any rights to the software except as expressly set forth in this 
% license. You agree that any copies of the software will contain the same proprietary notices 
% which appear on and in the software. 
% 
% Rent, lease, loan 
% 
% You may NOT rent, lease or loan the software without first negotiating a specific license
% for that purpose with Universit�t W�rzburg. 
%     
% No warranties 
% 
% Universit�t W�rzburg does NOT warrant that the software is error free. Universit�t W�rzburg 
% disclaims all warranties with respect to the software, either express or implied, including 
% but not limited to implied warranties of merchantability, fitness for a particular purpose and 
% noninfringement of third party rights. The software is provided "AS IS." 
% 
% No liability for consequential damages 
% 
% In no event will Universit�t W�rzburg be liable for any loss of profits, business, use, or data 
% or for any consequential, special, incidental or indirect damages of any kind arising out of 
% the delivery or performance or as a result of using or modifying the software, even if 
% Universit�t W�rzburg has been advised of the possibility of such damages. In no event will 
% Universit�t W�rzburg's liability for any claim, whether in contract, negligence, tort or any 
% other theory of liability, exceed the license fee paid by you, if any. 
% The licensed software is not designed for use in high-risk activities requiring fail-safe 
% performance. Universit�t W�rzburg disclaims any express or implied warranty of fitness for 
% high-risk activities. 
% 
% Severability 
% 
% In the event of invalidity of any provision of this license, the parties agree that such 
% invalidity shall not affect the validity of the remaining portions of this license.
% 
