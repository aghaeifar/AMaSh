function [SHIM, data] = shimMapping_v2
% author: Ali Aghaeifar (ali.aghaeifar <at> tuebingen.mpg.de)
%
% This program needs uipickfiles scripts: https://mathworks.com/matlabcentral/fileexchange/10867-uipickfiles--uigetfile-on-steroids
%
% Discription: 
%   This script recounstructs basismaps by reading acquired rawdata. All files must be selected in a correct order (first one is the reference) 
% output:
%   SHIM: contains basis maps
%   data: contains additional median parameters
%

%%
SHIM.STD_FOV      = 300;
SHIM.STD_NPixel   = 151;
    
if exist('\\10.41.60.40\Upload\9T_Ali\', 'dir') == 7
    startFolder = '\\10.41.60.40\Upload\9T_Ali\';
else
    startFolder = pwd;
end
data.FileName = uipickfiles('FilterSpec', startFolder, 'REFilter', '.dat');
if isempty(data.FileName) || (iscell(data.FileName) == 0)
    disp('No file is selected')
    return;
end

%% Load background B0 map as base
SHIM.STD_B0 = single(ones(SHIM.STD_NPixel, SHIM.STD_NPixel, SHIM.STD_NPixel));
SHIM.Raw_B0 = [];
SHIM.Raw_B0_unwrapped = [];
for i=1:numel(data.FileName)
    [Image, SHIM.mrprot, SHIM.SODA] = RecoImage(data.FileName{i});
    if SHIM.mrprot.Config.NEco ~= 2
       error('We only accept two echoes'); 
    end
    if SHIM.SODA.Dim == 2% measurment is 2D, swap dimension of slices with partitions
        Image = permute(Image, [1 2 5 4 3 6:ndims(Image)]); % Col, Lin, Par, Cha, Slc/Slb, ... -> Col, Lin, Slc/Slb, Cha, Par, ...
        for k=1:SHIM.SODA.NSlices
            SAG(:,:,k) = SHIM.SODA.Coords{k}.SAG(:,:,2); % 2 = center of the slice
            COR(:,:,k) = SHIM.SODA.Coords{k}.COR(:,:,2);
            TRA(:,:,k) = SHIM.SODA.Coords{k}.TRA(:,:,2);
        end
    else
        SAG = SHIM.SODA.Coords{1}.SAG;
        COR = SHIM.SODA.Coords{1}.COR;
        TRA = SHIM.SODA.Coords{1}.TRA;
    end
    data.MaxSep = max([SHIM.SODA.PixelSizeReadout, SHIM.SODA.PixelSizePhase, SHIM.SODA.PixelSizeSlice])/2; %Compute maximal sepration of data points in the original image   

%      base = applyHannFilt(squeeze(Image)); % Col, Lin, Slc/Slb/Par, Cha, Echo, Rep    
     base = squeeze(Image); % Col, Lin, Slc/Slb/Par, Cha, Echo, Rep
    
    % Use zero SHIM to compute coil combination weights
    [~,data.wfull,~] = openadapt_3D(permute(base(:,:,:,:,1,1),[4 1 2 3]),[3 3 3],true); % adaptive combine
    data.wfull = permute(data.wfull,[2 3 4 1]);
    baseComb   = squeeze(sum(bsxfun(@times, base, data.wfull), 4));
    % Initialize variables
    for cMeas = 2:SHIM.mrprot.Config.NRepMeas        
        temp = squeeze(angle(baseComb(:,:,:,2,cMeas) .* conj(baseComb(:,:,:,1,cMeas)) .* baseComb(:,:,:,1,1) .* conj(baseComb(:,:,:,2,1))));
        SHIM.Raw_B0 = cat(4, SHIM.Raw_B0, single(temp));
        temp = UnWrap_mex(single(temp));
        SHIM.Raw_B0_unwrapped = cat(4, SHIM.Raw_B0_unwrapped, single(temp));
        temp = temp /((SHIM.mrprot.MeasYaps.alTE{2}-SHIM.mrprot.MeasYaps.alTE{1})*1e-6)/(2*pi) / 0.5;
        temp = InterpolateToStandardFoV(temp, SAG, COR, TRA, SHIM.STD_FOV, SHIM.STD_NPixel, data.MaxSep);
        SHIM.STD_B0 = cat(4, SHIM.STD_B0, single(temp));
    end
    % clear temp baseComb base
end
SHIM.STD_B0 = single(SHIM.STD_B0);
fclose all;
    
SHIM.NShimCha = numel(data.FileName);
%% Mask
% data.RAW_Mask = imerode(bet2(single(sqrt(sum(sum(abs(data.seriesComb).^2,4),5)))), ones(1,1,1));
data.RAW_Mask = sqrt(sum(abs(baseComb(:,:,:,:,1)).^2,4));
data.RAW_Mask = data.RAW_Mask/max(data.RAW_Mask(:));
data.T = 10*mean(col(data.RAW_Mask(1:4,1:4,1:4)));
data.RAW_Mask(data.RAW_Mask<data.T) = 0;
data.RAW_Mask(data.RAW_Mask>=data.T) = 1; 
data.RAW_Mask = padarray(data.RAW_Mask,[1 1], 1, 'both');
data.RAW_Mask = imfill(data.RAW_Mask, 'holes');
data.RAW_Mask = imerode(data.RAW_Mask, ones(3,3,3));
data.RAW_Mask = double(bwareaopen(data.RAW_Mask, round(sum(col(data.RAW_Mask))/2)));
data.RAW_Mask = imdilate(data.RAW_Mask, ones(3,3,3));
data.RAW_Mask = imfill(data.RAW_Mask, 'holes');
data.RAW_Mask = imerode(data.RAW_Mask, ones(2,2,2));
data.RAW_Mask = data.RAW_Mask(2:end-1, 2:end-1, :);
vin(data.RAW_Mask)                                
SHIM.STD_Mask = InterpolateToStandardFoV(data.RAW_Mask, SAG, COR, TRA, SHIM.STD_FOV, SHIM.STD_NPixel, data.MaxSep); 
SHIM.STD_Mask(SHIM.STD_Mask>0.5) = 1;                                        
SHIM.STD_Mask(SHIM.STD_Mask<=0.5) = 0;  

%%
data.RAW_Mask =sqrt(sum(abs(baseComb(:,:,:,:,1)).^2,4));
data.RAW_Mask(data.RAW_Mask < .81) = 0;
data.RAW_Mask(data.RAW_Mask >= .81) = 1;
data.RAW_Mask = padarray(data.RAW_Mask,[1 1], 1, 'both');
data.RAW_Mask = imfill(data.RAW_Mask, 'holes');
data.RAW_Mask = data.RAW_Mask(2:end-1, 2:end-1,:);
data.RAW_Mask = imerode(data.RAW_Mask, ones(3,3,3));
% data.RAW_Mask = imdilate(data.RAW_Mask, strel('disk',3)); % imdilate(data.RAW_Mask, ones(5,5,5));
vin(data.RAW_Mask);
% vin(data.RAW_Mask .* baseComb(:,:,:,1));
SHIM.STD_Mask = InterpolateToStandardFoV(data.RAW_Mask, SAG, COR, TRA, SHIM.STD_FOV, SHIM.STD_NPixel, data.MaxSep); 
SHIM.STD_Mask(SHIM.STD_Mask>0.5) = 1;                                        
SHIM.STD_Mask(SHIM.STD_Mask<=0.5) = nan;  
SHIM.STD_Mask = single(SHIM.STD_Mask);

%%
data.RAW_Mask =sqrt(sum(abs(baseComb(:,:,:,:,1)).^2,4));
if isunix
    data.RAW_Mask = double(betFSL(abs(data.RAW_Mask)));
else
    data.RAW_Mask = bet2(single(abs(data.RAW_Mask)));
end
% data.RAW_Mask = imdilate(data.RAW_Mask, strel('disk',3)); % imdilate(data.RAW_Mask, ones(5,5,5));
vin(data.RAW_Mask);
% vin(data.RAW_Mask .* baseComb(:,:,:,1));
SHIM.STD_Mask = InterpolateToStandardFoV(data.RAW_Mask, SAG, COR, TRA, SHIM.STD_FOV, SHIM.STD_NPixel, data.MaxSep); 
SHIM.STD_Mask(SHIM.STD_Mask>0.5) = 1;                                        
SHIM.STD_Mask(SHIM.STD_Mask<=0.5) = nan;  
SHIM.STD_Mask = single(SHIM.STD_Mask);
%%
% temp = SHIM.STD_Mask; temp(temp==0) = NaN;
% SHIM.STD_B0 = bsxfun(@times, SHIM.STD_B0, temp);

SHIM.HardLimitsMax = [1500 5000*ones(1,3) 1.5*ones(1,16)];
SHIM.HardLimitsMin = -SHIM.HardLimitsMax;  
% temp = data.RAW_Mask; temp(temp==0) = NaN;
% SHIM.RAW_B0 = bsxfun(@times, data.RAW_B0, temp);

