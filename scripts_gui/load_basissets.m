% ------------------------------------------------------------------
% -------------------      Load Shim BasisMaps        --------------

function load_basissets(hObject)   
% Get correct handles (handles provided by the function arguments are not working)
handles = guidata(hObject);
setStatus(hObject, 1);
% Load file
[FileName,PathName,FilterIndex] = uigetfile('*.mat', 'Select shim basis-maps');
if FilterIndex
    load(fullfile(PathName, FileName), 'SHIM');
    if SHIM.STD_FOV ~= handles.settings.FOV || SHIM.STD_NPixel ~= handles.settings.NPixel
        h = msgbox(['Mismatch between resolution of shim-maps (' num2str(SHIM.STD_NPixel) ') and standard space (' num2str(handles.settings.NPixel) '). Try to update it...'], 'Error','error');  set(h, 'WindowStyle', 'modal')
        handles.settings.FOV = SHIM.STD_FOV;
        handles.settings.NPixel = SHIM.STD_NPixel;
        handles.flags.isDataLoaded = false; % Since the old imported B0 has different resolution
    end
    handles.SHIM = SHIM;
    handles.settings.HardLimitsMax = handles.SHIM.HardLimitsMax;
    handles.settings.HardLimitsMin = handles.SHIM.HardLimitsMin;     
    handles.settings.shim_channels = handles.SHIM.indMC; %[handles.SHIM.indMC(1) handles.SHIM.indMC(end)];
    if isempty(handles.settings.shim_channels)
        handles.settings.shim_channels = handles.SHIM.indSH;
    end
    handles.flags.isShimMapsLoaded = true;
    set(handles.editCrosshairPos, 'String', sprintf('%d %d %d', floor(handles.settings.NPixel/2 * ones(1,3))) );
    clear SHIM;   
end
setStatus(hObject, 0);
% Update handles
guidata(hObject, handles);