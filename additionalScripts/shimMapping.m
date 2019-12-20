function [SHIM, data] = shimMapping
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
[Image, SHIM.mrprot, SHIM.SODA] = RecoImage(data.FileName{1});
if strcmp(SHIM.mrprot.MeasYaps.sKSpace.ucDimension, '0x2') || isequal(SHIM.mrprot.MeasYaps.sKSpace.ucDimension, 2)% measurment is 2D, swap dimension of slices with partitions
    Image = permute(Image, [1 2 5 4 3 6:ndims(Image)]); % Col, Lin, Par, Cha, Slc/Slb, ... -> Col, Lin, Slc/Slb, Cha, Par, ...
end
base = applyHannFilt(squeeze(Image));
shim_setting = getShimSetting(SHIM.mrprot, 36);

% Use zero SHIM to compute coil combination weights
[~,data.wfull,~] = openadapt_3D(permute(base(:,:,:,:,1),[4 1 2 3]),[3 3 3],true); % adaptive combine
data.wfull = permute(data.wfull,[2 3 4 1]);

baseComb = zeros(size(base,1), size(base,2), size(base,3), SHIM.mrprot.Config.NEco);
for cEco = 1:SHIM.mrprot.Config.NEco % Number of echos
    baseComb(:,:,:,cEco) = sum(bsxfun(@times, base(:,:,:,:,cEco), data.wfull), 4);
end

%
% Initialize variables
baseCombSz    = size(baseComb);
shimChaComb   = zeros(baseCombSz);
data.Shims    = ones([baseCombSz(1:3) SHIM.mrprot.Config.NEco-1 numel(data.FileName)]);
data.RAW_B0   = zeros(size(data.Shims));
data.MaxSep   = max([SHIM.SODA.PixelSizeReadout, SHIM.SODA.PixelSizePhase, SHIM.SODA.PixelSizeSlice])/2; %Compute maximal sepration of data points in the original image
SHIM.NRXCha   = SHIM.mrprot.Config.NChaMeas;
SHIM.NShimCha = numel(data.FileName);
SHIM.scaling  = [1 50 50 50 ones(1,32)];
SHIM.STD_FOV      = 300;
SHIM.STD_NPixel   = 151;
SHIM.STD_B0   = ones(SHIM.STD_NPixel, SHIM.STD_NPixel, SHIM.STD_NPixel, SHIM.mrprot.Config.NEco-1, numel(data.FileName));

if SHIM.SODA.Dim==2  % Stupid 2D data
    for i=1:SHIM.SODA.NSlices
        SAG(:,:,i) = SHIM.SODA.Coords{i}.SAG(:,:,2); % 2 = center of the slice
        COR(:,:,i) = SHIM.SODA.Coords{i}.COR(:,:,2);
        TRA(:,:,i) = SHIM.SODA.Coords{i}.TRA(:,:,2);
    end
else
    SAG = SHIM.SODA.Coords{1}.SAG;
    COR = SHIM.SODA.Coords{1}.COR;
    TRA = SHIM.SODA.Coords{1}.TRA;
end

%
for i=2:numel(data.FileName) % start after base
    SHIM.Name{i} = ['Ch' num2str(i-1)];
    disp(SHIM.Name{i});
    [Image, SHIM.mrprot, SHIM.SODA] = RecoImage(data.FileName{i});
    shim_setting2 = getShimSetting(SHIM.mrprot, 36);
    if shim_setting2.SH ~=  shim_setting.SH
       [shim_setting.SH', shim_setting2.SH']
       error('Shim settings have been changed');        
    end
    if strcmp(SHIM.mrprot.MeasYaps.sKSpace.ucDimension, '0x2') || isequal(SHIM.mrprot.MeasYaps.sKSpace.ucDimension, 2)% measurment is 2D, swap dimension of slices with partitions
        Image = permute(Image, [1 2 5 4 3 6:ndims(Image)]); % Col, Lin, Par, Cha, Slc/Slb, ... -> Col, Lin, Slc/Slb, Cha, Par, ...
    end
    shimCha = applyHannFilt(squeeze(Image));
    for cEco = 1:SHIM.mrprot.Config.NEco % Number of echos
        shimChaComb(:,:,:,cEco)=sum(bsxfun(@times,shimCha(:,:,:,:,cEco),data.wfull),4);
    end
    for cEco = 1:SHIM.mrprot.Config.NEco-1
        data.Shims(:,:,:,cEco,i) = angle(baseComb(:,:,:,cEco).* conj(baseComb(:,:,:,cEco+1))... % No SHIM
                                  .*conj(shimChaComb(:,:,:,cEco)).*shimChaComb(:,:,:,cEco+1));
        data.Shims(:,:,:,cEco,i) = UnWrap_mex(single(data.Shims(:,:,:,cEco,i))); 
        data.RAW_B0(:,:,:,cEco,i)= data.Shims(:,:,:,cEco,i) /((SHIM.mrprot.MeasYaps.alTE{cEco+1}-SHIM.mrprot.MeasYaps.alTE{cEco})*1e-6)/(2*pi);
    end                   
    data.RAW_B0(:,:,:,:,i) = data.RAW_B0(:,:,:,:,i)/50; %/SHIM.scaling(i); % Apply scaling factors
    SHIM.STD_B0(:,:,:,:,i) = InterpolateToStandardFoV(data.RAW_B0(:,:,:,:,i),SAG, COR, TRA, SHIM.STD_FOV, SHIM.STD_NPixel, data.MaxSep); 
end
SHIM.STD_B0 = squeeze(SHIM.STD_B0);
fclose all;
       
%% Mask
% data.RAW_Mask = imerode(bet2(single(sqrt(sum(sum(abs(data.seriesComb).^2,4),5)))), ones(1,1,1));
data.RAW_Mask = sqrt(sum(abs(baseComb).^2,4));
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
data.RAW_Mask =sqrt(sum(abs(baseComb).^2,4));
data.RAW_Mask(data.RAW_Mask < 0.01) = 0;
data.RAW_Mask(data.RAW_Mask >= 0.01) = 1;
data.RAW_Mask = padarray(data.RAW_Mask,[1 1], 1, 'both');
data.RAW_Mask = imfill(data.RAW_Mask, 'holes');
data.RAW_Mask = data.RAW_Mask(2:end-1, 2:end-1,:);
data.RAW_Mask = imerode(data.RAW_Mask, ones(3,3,3));
% data.RAW_Mask = imdilate(data.RAW_Mask, strel('disk',3)); % imdilate(data.RAW_Mask, ones(5,5,5));
vin(data.RAW_Mask);
% vin(data.RAW_Mask .* baseComb(:,:,:,1));
SHIM.STD_Mask = InterpolateToStandardFoV(data.RAW_Mask, SAG, COR, TRA, SHIM.STD_FOV, SHIM.STD_NPixel, data.MaxSep); 
SHIM.STD_Mask(SHIM.STD_Mask>0.5) = 1;                                        
SHIM.STD_Mask(SHIM.STD_Mask<=0.5) = 0;  

%%
% temp = SHIM.STD_Mask; temp(temp==0) = NaN;
% SHIM.STD_B0 = bsxfun(@times, SHIM.STD_B0, temp);
SHIM.STD_Mask(SHIM.STD_Mask == 0) = nan;
SHIM.STD_B0 = single(SHIM.STD_B0);
SHIM.STD_Mask = single(SHIM.STD_Mask);
SHIM.HardLimitsMax = [1500 5000*ones(1,3) 4000*ones(1,5) 1.5*ones(1,16)];
SHIM.HardLimitsMin = -SHIM.HardLimitsMax;  
% temp = data.RAW_Mask; temp(temp==0) = NaN;
% SHIM.RAW_B0 = bsxfun(@times, data.RAW_B0, temp);

%%
load('D:\Dropbox\Matlab\B0GUI\basisMaps\SHIM_MC32opt_SH1_2D_201p_21.09.2018_Siemens_multiRef.mat');

SHIM.mrprot = SHIM_01_04.mrprot;
SHIM.SODA   = SHIM_01_04.SODA;
SHIM.STD_NPixel = 151;

SHIM.STD_B0 = zeros(151,151,151,36);
SHIM.STD_B0(:,:,:,1:4)   = squeeze(SHIM_X_Y_Z.STD_B0);
SHIM.STD_B0(:,:,:,5:8)   = squeeze(SHIM_01_04.STD_B0(:,:,:,:,2:end));
SHIM.STD_B0(:,:,:,9:16)  = squeeze(SHIM_05_12.STD_B0(:,:,:,:,2:end));
SHIM.STD_B0(:,:,:,17:24) = squeeze(SHIM_13_20.STD_B0(:,:,:,:,2:end));
SHIM.STD_B0(:,:,:,25:32) = squeeze(SHIM_21_28.STD_B0(:,:,:,:,2:end));
SHIM.STD_B0(:,:,:,33:36) = squeeze(SHIM_29_32.STD_B0(:,:,:,:,2:end));

SHIM.STD_Mask = SHIM_01_04.STD_Mask;
