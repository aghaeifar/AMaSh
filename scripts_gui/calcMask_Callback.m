% ------------------------------------------------------------------
% -------------------       Brain Mask           -------------------

function calcMask_Callback(hObject)   

handles =  guidata(hObject);
if handles.flags.isDataLoaded == false
    h = msgbox('Rawdata is not loaded!', 'Warning','Warn');
    set(h, 'WindowStyle', 'modal')
    return;
end
setStatus(hObject, 1);
% Get values from GUI
% Threshold       = str2num(get(handles.editThreshold,'String'));  % [Percent]
% Steepness       = str2num(get(handles.editSteepness,'String'));  % [Percent]
% ErodeKernelSize = str2num(get(handles.editKernelSize,'String')); % [Pixel]

% Calculate mask with class method
mask_methods = {'none', 'bet', 'threshold'}; % the orders must be identical with setting GUI
handles.B0_OBJ.mask_mode = mask_methods(handles.settings.mask_algorithm_no);
handles.B0_OBJ  = CalcMask(handles.B0_OBJ, handles.settings.mask_erosion_sz );
handles.DISPLAY.Mask = handles.B0_OBJ.STD_Mask;

handles.flags.enMask = true;
set(handles.toolbar_use_mask, 'State', 'on');
% Update
guidata(hObject,handles);
UpdateShimmedField(hObject)     % Always must be called before UpdateMainImages
UpdateB0Map(hObject);
UpdateMainImages(hObject);
setStatus(hObject, 0);