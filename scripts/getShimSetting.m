function ShimSetting = getShimSetting(twix, nCha)

if isfield(twix, 'hdr')
    mrprot = twix.hdr;
elseif isfield(twix, 'MeasYaps')
     mrprot = twix;
else
    error('Wrong input');
end
    
ShimSetting.A00 = mrprot.MeasYaps.sTXSPEC.asNucleusInfo{1}.lFrequency;

% X
try
    ShimSetting.A11 = mrprot.Phoenix.sGRADSPEC.lOffsetX.*mrprot.Phoenix.sGRADSPEC.flSensitivityX*1000;
catch
    ShimSetting.A11 = 0; % sometimes 'lOffsetX' doesn't exist in rawdata
end
% Y
try
    ShimSetting.B11 = mrprot.Phoenix.sGRADSPEC.lOffsetY.*mrprot.Phoenix.sGRADSPEC.flSensitivityY*1000;
catch
    ShimSetting.B11 = 0;
end
% Z
try
    ShimSetting.A10 = mrprot.Phoenix.sGRADSPEC.lOffsetZ.*mrprot.Phoenix.sGRADSPEC.flSensitivityZ*1000;
catch
    ShimSetting.A10 = 0;
end

% Run editmeasperm.exe on the scanner and type "d flshim". Record the values of aflSHIMSensitivity[0-4] and flSHIMRefRadius
% sGRADSPEC.alShimCurrent[0-4]*aflSHIMSensitivity[0-4]/1000/(flSHIMRefRadius^2)
% Z2

Z2ConvFactor = 1000/1335;
ShimSetting.A20 = mrprot.Phoenix.sGRADSPEC.alShimCurrent{1}*Z2ConvFactor;
if isempty(ShimSetting.A20)
    ShimSetting.A20 = 0;
end
% ZX
ZXConvFactor = 1010/2697;
ShimSetting.A21 = mrprot.Phoenix.sGRADSPEC.alShimCurrent{2}*ZXConvFactor;
if isempty(ShimSetting.A21)
    ShimSetting.A21 = 0;
end
% ZY
ZYConvFactor = 1020/2759;
ShimSetting.B21 = mrprot.Phoenix.sGRADSPEC.alShimCurrent{3}*ZYConvFactor;
if isempty(ShimSetting.B21)
    ShimSetting.B21 = 0;
end
% X2-Y2
X2Y2ConvFactor = 1030/2786;
ShimSetting.A22 = mrprot.Phoenix.sGRADSPEC.alShimCurrent{4}*X2Y2ConvFactor;
if isempty(ShimSetting.A22)
    ShimSetting.A22 = 0;
end
% XY
XYConvFactor = 1040/2850;
ShimSetting.B22 = mrprot.Phoenix.sGRADSPEC.alShimCurrent{5}*XYConvFactor;
if isempty(ShimSetting.B22)
    ShimSetting.B22 = 0;
end

ShimSetting.SH = [ShimSetting.A00, ShimSetting.A11, ShimSetting.B11, ShimSetting.A10, ShimSetting.A20, ShimSetting.A21, ShimSetting.B21, ShimSetting.A22, ShimSetting.B22];
if nargin > 1 && isfield(twix, 'image')
    % slcNo = twix.hdr.MeasYaps.sSliceArray.lSize; % in the case of 3D, it constains number of the slabs
    try
        iceParam = twix.image.iceParam(:,1:twix.image.dataSize(8):end); % twix.image.dataSize(8) = echo
        iceParam = iceParam(:,1:twix.image.dataSize(4):end); % twix.image.dataSize(4) = slices per slab(partitions)
        iceparam_sz1 = size(twix.image.iceParam, 1); % it is 4 for VB17
        nCha_d4 = ceil(nCha/iceparam_sz1);
        iceParam = reshape(iceParam, [iceparam_sz1 twix.image.sqzSize(3) twix.image.dataSize(5)]); % twix.image.sqzSize(3) = Line  & twix.image.dataSize(5) = Slice
        iceParam = transpose(reshape(iceParam(:,1:nCha_d4,:), [iceparam_sz1*nCha_d4 twix.image.dataSize(5)]));
        iceParam(iceParam>2^15) = iceParam(iceParam>2^15) - 2^16; % twix.image.iceParam is unsigned int 
    catch
        iceParam = zeros(twix.image.dataSize(5),nCha);
    end
    ShimSetting.Dynamic = [iceParam(:,1) iceParam(:,2:4) iceParam(:,5:end)/1000];  
end
