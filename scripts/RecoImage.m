function [Image, mrprot, SODA, twix] = RecoImage(filename)

%% Map raw data
%
if(nargin==0)
   twix = mapVBVD();
else
   twix = mapVBVD(filename);
end

if iscell(twix) % for VD
    twix = twix{2};
end

twix.image.flagDoAverage = true;
mrprot  = twix.hdr;
data    = twix.image();
SODA    = SODA_OBJ('mrprot',mrprot);  % Calculate SODA 

%% Partial Fourier
if(~strcmp(mrprot.MeasYaps.sKSpace.ucSlicePartialFourier,'0x10') && sum(mrprot.MeasYaps.sKSpace.ucSlicePartialFourier ~= 16)) 
    error('Slice Partial Fourier')
end
if(~strcmp(mrprot.MeasYaps.sKSpace.ucPhasePartialFourier,'0x10') && sum(mrprot.MeasYaps.sKSpace.ucPhasePartialFourier ~= 16)) 
    error('Phase Partial Fourier')
end
% check if zero filling is required
if twix.image.NCol ~= twix.hdr.Meas.iRoFTLength || twix.image.NLin ~= twix.hdr.Meas.iPEFTLength
    shift = floor([twix.hdr.Meas.iRoFTLength - twix.image.NCol, twix.hdr.Meas.iPEFTLength - twix.image.NLin] /2);
    data  = padarray(data,[shift(1), 0, shift(2)], 0, 'both'); % fill with zero - this function is part of Image Processing Toolbox
end

%% Fourier transform
data = permute(data,[1 3 4 2 5 6 7 8 9]);  % -> Col, Lin, Par, Cha, Slc, Ave, Phs, Eco, Rep, ...
data = flip(flip(data,1),2);

if strcmp(mrprot.MeasYaps.ucReadOutMode, '0x2') % bipolar readout
    warning('Bipolar Readout')
%     data(:,:,:,:,:,:,:,2:2:end,:,:,:,:,:,:) = flip(data(:,:,:,:,:,:,:,2:2:end,:,:,:,:,:,:),1);
end

if(SODA.Dim==2)
    Image = FFT2D(data);
else
    Image = FFT3D(data);
end

Image = ReorderSlab(Image, twix); % 3D acquisition with one slab doesn't require reorder

%% Cut of oversampling
% CutOff = twix.image.NCol/4; % this doesn't work when zero filling is required
CutOff = (size(Image, 1) - twix.hdr.MeasYaps.sKSpace.lBaseResolution)/2;
% CutOff = twix.hdr.MeasYaps.sKSpace.lBaseResolution/2; % twix.hdr.MeasYaps.sKSpace.lBaseResolution == twix.hdr.Meas.iRoFTLength
Image = Image(CutOff+1:end-CutOff,:,:,:,:,:,:,:,:);
Image = flip(flip(Image,1),2);
