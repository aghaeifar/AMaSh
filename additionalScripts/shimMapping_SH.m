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
clear SHIM
SHIM.STD_FOV      = 300;
SHIM.STD_NPixel   = 151;
    
root = 'D:\Projects\Shimming\Basis_Maps\Rawdata\2ndOrderSH_basismaps_25.06.2019';
SHIM.files.Ref = 'meas_MID732_Ref_Sag_2_0mm_FID17820.dat';
SHIM.files.X   = 'meas_MID734_X_50_Sag_2_0mm_FID17822.dat';
SHIM.files.Y   = 'meas_MID737_Y50_Sag_2_0mm_FID17825.dat';
SHIM.files.Z   = 'meas_MID740_Z50_Sag_2_0mm_FID17828.dat';
SHIM.files.A20 = 'meas_MID743_A20_50_Sag_2_0mm_FID17831.dat';
SHIM.files.A21 = 'meas_MID745_A21_50_Sag_2_0mm_FID17833.dat';
SHIM.files.B21 = 'meas_MID747_B21__50_Sag_2_0mm_FID17835.dat';
SHIM.files.A22 = 'meas_MID749_A22__50_Sag_2_0mm_FID17837.dat';
SHIM.files.B22 = 'meas_MID751_B22_50_Sag_2_0mm_FID17839.dat';
scales_applied = [-50, 50, 50, 50, 50, -50, -50, 50];
field_name = {'Ref', 'X', 'Y', 'Z', 'A20', 'A21', 'B21', 'A22', 'B22'};

SHIM.HardLimitsMin  = [9999 -5000 -5000 -5000 -9359 -4679 -4620 -4620 -4559];
SHIM.HardLimitsMax  = [9999 5000 5000 5000 9359 4679 4620 4620 4559];
SHIM.indSH = 1:9;
SHIM.indMC = [];
SHIM.Name = field_name;

%% Load background B0 map as base
[base, SHIM.mrprot, SHIM.SODA, SHIM.SAG, SHIM.COR, SHIM.TRA] = simple_reco(fullfile(root, SHIM.files.Ref));
ShimSetting_ref = getShimSetting(SHIM.mrprot);

SHIM.MaxSep = max([SHIM.SODA.PixelSizeReadout, SHIM.SODA.PixelSizePhase, SHIM.SODA.PixelSizeSlice])/2; %Compute maximal sepration of data points in the original image   
[~,wfull,~] = openadapt_3D(permute(base(:,:,:,:,1,1),[4 1 2 3]),[3 3 3],true); % adaptive combine
wfull = permute(wfull,[2 3 4 1]);
baseComb   = squeeze(sum(bsxfun(@times, base, wfull), 4));

SHIM.STD_B0 = ones(SHIM.STD_NPixel, SHIM.STD_NPixel, SHIM.STD_NPixel);
SHIM.Raw_B0 = ones(size(baseComb,1), size(baseComb,2), size(baseComb,3));
SHIM.Raw_B0_unwrapped = SHIM.Raw_B0 ;
%% Reconstruct basis-maps
ShimSettings = zeros(numel(field_name), 9);
for i=2:numel(field_name)
    if isempty(field_name{i})
       continue; 
    end
    [raw_data, mrprot] = simple_reco(fullfile(root, SHIM.files.(field_name{i})));
    temp = getShimSetting(mrprot);
    ShimSettings(i,:) = temp.SH;
    scales = ShimSettings(i,i) - ShimSetting_ref.SH(i); % i=1 is frequencey
    
    raw_data = squeeze(sum(bsxfun(@times, raw_data, wfull), 4));
    temp = squeeze(angle(raw_data(:,:,:,2) .* conj(raw_data(:,:,:,1)) .* baseComb(:,:,:,1) .* conj(baseComb(:,:,:,2))));
    SHIM.Raw_B0 = cat(4, SHIM.Raw_B0, single(temp));
    temp = UnWrap_mex(single(temp));
    SHIM.Raw_B0_unwrapped = cat(4, SHIM.Raw_B0_unwrapped, single(temp));
    temp = temp /((SHIM.mrprot.MeasYaps.alTE{2}-SHIM.mrprot.MeasYaps.alTE{1})*1e-6)/(2*pi) / scales;
    temp = InterpolateToStandardFoV(temp, SHIM.SAG, SHIM.COR, SHIM.TRA, SHIM.STD_FOV, SHIM.STD_NPixel, SHIM.MaxSep);
    SHIM.STD_B0 = cat(4, SHIM.STD_B0, temp);
end
SHIM.STD_B0 = single(SHIM.STD_B0);
SHIM.Raw_B0 = single(SHIM.Raw_B0);
SHIM.Raw_B0_unwrapped = single(SHIM.Raw_B0_unwrapped);
SHIM.NShimCha = size(SHIM.STD_B0, 4);
fclose all;
    
%% Masking Method 1
data.RAW_Mask = sqrt(sum(abs(baseComb(:,:,:,:,1)).^2,4));
data.RAW_Mask = data.RAW_Mask/max(data.RAW_Mask(:));
data.T = 30*mean(col(data.RAW_Mask(1:4,1:4,1:4)));
data.RAW_Mask(data.RAW_Mask<data.T) = 0;
data.RAW_Mask(data.RAW_Mask>=data.T) = 1; 
data.RAW_Mask = padarray(data.RAW_Mask,[1 1], 1, 'both');
data.RAW_Mask = imfill(data.RAW_Mask, 'holes');
data.RAW_Mask = imerode(data.RAW_Mask, ones(3,3,3));
data.RAW_Mask = double(bwareaopen(data.RAW_Mask, double(round(sum(col(data.RAW_Mask))/2))));
data.RAW_Mask = imdilate(data.RAW_Mask, ones(3,3,3));
data.RAW_Mask = imfill(data.RAW_Mask, 'holes');
data.RAW_Mask = imerode(data.RAW_Mask, ones(2,2,2));
data.RAW_Mask = data.RAW_Mask(2:end-1, 2:end-1, :);
vin(data.RAW_Mask)                                
SHIM.STD_Mask = InterpolateToStandardFoV(data.RAW_Mask, SHIM.SAG, SHIM.COR, SHIM.TRA, SHIM.STD_FOV, SHIM.STD_NPixel, SHIM.MaxSep); 
SHIM.STD_Mask(SHIM.STD_Mask>0.5) = 1;                                        
SHIM.STD_Mask(SHIM.STD_Mask<=0.5) = nan;  
SHIM.STD_Mask = single(SHIM.STD_Mask);

%% Masking Method 2
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
SHIM.STD_Mask = InterpolateToStandardFoV(data.RAW_Mask, SHIM.SAG, SHIM.COR, SHIM.TRA, SHIM.STD_FOV, SHIM.STD_NPixel, SHIM.MaxSep); 
SHIM.STD_Mask(SHIM.STD_Mask>0.5) = 1;                                        
SHIM.STD_Mask(SHIM.STD_Mask<=0.5) = nan;  
SHIM.STD_Mask = single(SHIM.STD_Mask);

%% Masking Method 3
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
function [base, mrprot, SODA, SAG, COR, TRA] = simple_reco(addr)
    [Image, mrprot, SODA] = RecoImage(addr);
    if mrprot.Config.NEco ~= 2
       error('We only accept two echoes'); 
    end
    if SODA.Dim == 2% measurment is 2D, swap dimension of slices with partitions
        Image = permute(Image, [1 2 5 4 3 6:ndims(Image)]); % Col, Lin, Par, Cha, Slc/Slb, ... -> Col, Lin, Slc/Slb, Cha, Par, ...
        for k=1:SODA.NSlices
            SAG(:,:,k) = SODA.Coords{k}.SAG(:,:,2); % 2 = center of the slice
            COR(:,:,k) = SODA.Coords{k}.COR(:,:,2);
            TRA(:,:,k) = SODA.Coords{k}.TRA(:,:,2);
        end
    else
        SAG = SODA.Coords{1}.SAG;
        COR = SODA.Coords{1}.COR;
        TRA = SODA.Coords{1}.TRA;
    end
    base = squeeze(Image); % Col, Lin, Slc/Slb/Par, Cha, Echo, Rep
end
