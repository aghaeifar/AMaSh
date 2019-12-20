% ------------------------------------------------------------------
% -------------------     Calculating Shim coefficients     --------

function calcShim_Callback(hObject)

handles = guidata(hObject);
% Check requirements 
if ~handles.flags.isShimMapsLoaded || ~handles.flags.SODALoaded || ~handles.flags.isDataLoaded
    h = msgbox(['BasisMaps are loaded = ' num2str(handles.flags.isShimMapsLoaded) ', SODA is loaded = ' num2str(handles.flags.SODALoaded)], 'Warning','Warn'); set(h, 'WindowStyle', 'modal');
    return;
end

setStatus(hObject, 1);
handles.flags.isDuplicateSession = false; % remove duplicate flag if the session is duplicated
guidata(hObject, handles);

myInterpolate(hObject, true, handles.settings.adjacent_slices); % true = interpolate basis-maps too
handles = guidata(hObject);
shimCha = handles.settings.shim_channels; %handles.settings.shim_channels(1):handles.settings.shim_channels(2);   % Select shim channels
handles.TARGET.DeltaShimSettings = zeros(handles.TARGET.ShimBox.NSlices, handles.SHIM.NShimCha);
handles.TARGET.B0_SHIMMED        = cell(handles.TARGET.ShimBox.NSlices, 1);

cprintf('text', 'Calculating shim currents...\n')
for cSlc = 1:handles.TARGET.ShimBox.NSlices
    MAP_B0_MSK    = handles.TARGET.MAP_B0{cSlc}.* handles.TARGET.Mask{cSlc}; % It's very important to multiply Mask here, 0 -> nan
%     try
        SoftLimitsMin = handles.settings.HardLimitsMin(shimCha) - handles.B0_OBJ.SHIM.Dynamic(abs(handles.B0_OBJ.SODA.Dim-3)*cSlc + (handles.B0_OBJ.SODA.Dim-2), shimCha);
        SoftLimitsMax = handles.settings.HardLimitsMax(shimCha) - handles.B0_OBJ.SHIM.Dynamic(abs(handles.B0_OBJ.SODA.Dim-3)*cSlc + (handles.B0_OBJ.SODA.Dim-2), shimCha);
%     catch
%         SoftLimitsMin = handles.settings.HardLimitsMin(shimCha);
%         SoftLimitsMax = handles.settings.HardLimitsMax(shimCha);
%     end
    [handles.TARGET.DeltaShimSettings(cSlc, shimCha), handles.TARGET.B0_SHIMMED{cSlc}] = shimGlobal(MAP_B0_MSK(:),...
                                                                                                    handles.TARGET.MAP_Shims{cSlc}(:,shimCha),...
                                                                                                    handles.settings.shim_algorithm_no,...
                                                                                                    SoftLimitsMin, SoftLimitsMax);                                                                                                
    handles.TARGET.B0_SHIMMED{cSlc} = reshape(handles.TARGET.B0_SHIMMED{cSlc}, size(handles.TARGET.MAP_B0{cSlc}));
    if isfield(handles.SHIM, 'scale')
        handles.TARGET.DeltaShimSettings(cSlc, shimCha) = handles.TARGET.DeltaShimSettings(cSlc, shimCha) .* handles.SHIM.scale(shimCha);
    end
end
cprintf('green', 'Done\n');

handles.flags.isCalcShimDone = true;
guidata(hObject, handles);
afterCalc(hObject); 
setStatus(hObject, 0);
end









