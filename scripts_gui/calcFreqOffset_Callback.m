% ------------------------------------------------------------------
% -----------------     Calculating Freq Offset      ---------------

function calcFreqOffset_Callback(hObject)

handles = guidata(hObject);
if ~handles.flags.isShimMapsLoaded || ~handles.flags.SODALoaded
    disp(['BasisMaps are loaded=' num2str(handles.flags.isShimMapsLoaded) ', SODA is loaded=' num2str(handles.flags.SODALoaded)]);
    h = msgbox(['BasisMaps are loaded=' num2str(handles.flags.isShimMapsLoaded) ', SODA is loaded=' num2str(handles.flags.SODALoaded)], 'Warning','Warn');
    set(h, 'WindowStyle', 'modal')
    return;
end

setStatus(hObject, 1);
% for calculation of offset we don't employ neighbor pixels
myInterpolate(hObject, false, 0);
handles = guidata(hObject);
handles.TARGET.DeltaShimSettings = zeros(handles.TARGET.slcSize, handles.SHIM.NShimCha);
handles.TARGET.B0_SHIMMED = handles.TARGET.MAP_B0.*handles.TARGET.Mask;

for cSlc = 1:handles.TARGET.slcSize
    handles.TARGET.DeltaShimSettings(cSlc,1) = -nanmean(col(handles.TARGET.B0_SHIMMED(:,:,:,cSlc)));
    if isnan(handles.TARGET.DeltaShimSettings(cSlc,1)) % for the slices not covered by the mask
        handles.TARGET.DeltaShimSettings(cSlc,1) = 0;
        cprintf('_red', 'Slab/Slice with empty ROI detected.\n')
    end
    handles.TARGET.B0_SHIMMED(:,:,:,cSlc) = handles.TARGET.B0_SHIMMED(:,:,:,cSlc) + handles.TARGET.DeltaShimSettings(cSlc,1);
end
handles.flags.isCalcShimDone = true;
guidata(hObject, handles);
afterCalc(hObject);
setStatus(hObject, 0);