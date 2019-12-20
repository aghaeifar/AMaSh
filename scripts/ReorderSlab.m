function [B, Reorder] = ReorderSlab(A, twix)

B = [];
slcOrder = twix.hdr.Meas.chronSliceIndices(twix.hdr.Meas.chronSliceIndices>-1)+1; % +1 -> Index starts from 1
Reorder = zeros(1, numel(slcOrder));
for i=1:numel(slcOrder)
    Reorder(i) = find(slcOrder==i); 
end
if nargin > 1
    B = A(:,:,:,:,Reorder,:,:,:,:,:,:,:,:,:,:,:);
end