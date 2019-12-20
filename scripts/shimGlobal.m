function [coef, mapB0_shimmed, mean, std] = shimGlobal(mapB0, mapShim, opt_method, LL, UL)

%% Startup checking
if nargin<2
   error('No input') 
end

sz_mapB0   = size(mapB0);
sz_mapShim = size(mapShim);
nd         = ndims(mapShim);
if numel(LL) == 1
    sz_mapShim(end+1) = 1; 
end
% if ~isequal(sz_mapB0, sz_mapShim(1:end-1))
%     error('Mismatch in the size between B0map and shim basis-maps ')
% end

if nargin == 2
   opt_method = 1;
end

mapB0_shimmed = mapB0;
coef = zeros(sz_mapShim(end),1);
if sum(~isnan(mapB0)) == 0
    cprintf('_red', 'Slab/Slice with empty ROI detected.\n') 
    mean = nan;
    std = nan;
    return;
end

%%
% Prepare data in row format for solver

mapShimT = double(permute(mapShim, [nd 1:nd-1])); % change to [Cha Lin Col Slc]
mapShimT = reshape(mapShimT, size(mapShimT, 1), []);
mapB0T = double(reshape(mapB0, 1, []));

% Remove NAN values
mapShimT(:, isnan(mapB0T)) = [];
mapB0T(isnan(mapB0T)) = [];
if sum(isnan(mapShimT(1,:))) > 0
    cprintf('_red', 'basis-maps volume is smaller than ROI (%.1f%%).\n', 100*sum(isnan(mapShimT(1,:)))/numel(mapB0T));
end
mapB0T(isnan(mapShimT(1,:))) = [];
mapShimT(:, isnan(mapShimT(1,:))) = [];

%% Calc shim coef
switch opt_method
    case 1
        coef = calcShim_pseudoInverse(mapB0T, mapShimT, 'inv'); 
    case 2
        coef = calcShim_pseudoInverse(mapB0T, mapShimT, 'svd'); 
    case 3
        coef = calcShim_consTru(mapB0T, mapShimT, LL, UL); 
    case 4
        coef = calcShim_lsqlin(mapB0T, mapShimT, LL, UL, 'off'); 
    case 5 % fmincon (interior-point)
        coef = calcShim_fmincon(mapB0T, mapShimT, LL, UL, 1);
    case 6 % fmincon (sqp)
        coef = calcShim_fmincon(mapB0T, mapShimT, LL, UL, 2, 'off');
    case 7 % fmincon_std
        coef = calcShim_fmincon(mapB0T, mapShimT, LL, UL, 3);
    otherwise
        coef = zeros(sz_mapShim(end),1);
end

%% 
if nargout > 1
    mapShim_sum = repmat(permute(coef, [2:nd 1]), [sz_mapShim(1:end-1) 1]).* mapShim;
    mapShim_sum = squeeze(sum(mapShim_sum, numel(sz_mapShim)));
    mapB0_shimmed = mapB0 + mapShim_sum;
end
%%
if nargout > 2
    mean =  nanmean(mapB0_shimmed(:));
end
if nargout > 3
    std  =  nanstd(mapB0_shimmed(:));
end


