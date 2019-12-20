function varargout = mySettings(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mySettings_OpeningFcn, ...
                   'gui_OutputFcn',  @mySettings_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end




function mySettings_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);
handles = guidata(hObject);
if ~isempty(varargin)
    handles.hObject_b = varargin{1};
    handles2 = guidata(varargin{1});  
	set(handles.edit_shim_channel_a,        'String', num2str(handles2.settings.shim_channels));
    % set(handles.edit_shim_channel_a,        'String', num2str(handles2.settings.shim_channels(1)));
    % set(handles.edit_shim_channel_b,        'String', num2str(handles2.settings.shim_channels(2)));
    set(handles.popupmenu_shim_algorithm,   'Value' , handles2.settings.shim_algorithm_no);
    set(handles.edit_adjacent_slices,       'String', num2str(handles2.settings.adjacent_slices));
    set(handles.popupmenu_reco_mode,        'Value' , handles2.settings.reco_mode);
    set(handles.popupmenu_unwrapping_method,'Value' , handles2.settings.unwrapping_method);
    set(handles.edit_save_path,             'String', handles2.settings.path2save_coef );
    set(handles.edit_soda_path,             'String', handles2.settings.path2soda);
    set(handles.popupmenu_mask_algorithm,   'Value' , handles2.settings.mask_algorithm_no);
    set(handles.edit_mask_erosion_size,     'String', num2str(handles2.settings.mask_erosion_sz));  
    set(handles.edit_interpolation_fov,     'String', num2str(handles2.settings.FOV));
    set(handles.edit_interpolation_npixels, 'String', num2str(handles2.settings.NPixel));
    set(handles.edit_ampl_upper_boundry,    'String', num2str(handles2.settings.HardLimitsMax, '%2.1f  '));
    set(handles.edit_ampl_lower_boundry,    'String', num2str(handles2.settings.HardLimitsMin, '%2.1f  '));
end
guidata(hObject, handles);


function varargout = mySettings_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


function pushbutton_ok_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    update_settings(hObject, handles.hObject_b);
    delete(gcbf);

function pushbutton_cancel_Callback(hObject, eventdata, handles)
    delete(gcbf);

function pushbutton_apply_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    update_settings(hObject, handles.hObject_b);
    

function update_settings(hObject, hObject_b)
handles = guidata(hObject);
handles_b = guidata(hObject_b);

handles_b.settings.shim_channels     = str2num(get(handles.edit_shim_channel_a, 'String'));
%handles_b.settings.shim_channels(1)  = str2num(get(handles.edit_shim_channel_a, 'String'));
%handles_b.settings.shim_channels(2)  = str2num(get(handles.edit_shim_channel_b, 'String'));
handles_b.settings.shim_algorithm_no = get(handles.popupmenu_shim_algorithm, 'Value');
handles_b.settings.reco_mode         = get(handles.popupmenu_reco_mode, 'Value');
handles_b.settings.unwrapping_method = get(handles.popupmenu_unwrapping_method, 'Value');
handles_b.settings.adjacent_slices   = str2num(get(handles.edit_adjacent_slices, 'String'));
handles_b.settings.path2save_coef    = get(handles.edit_save_path, 'String');
handles_b.settings.path2soda         = get(handles.edit_soda_path, 'String');
handles_b.settings.mask_algorithm_no = get(handles.popupmenu_mask_algorithm, 'Value');
handles_b.settings.mask_erosion_sz   = str2num(get(handles.edit_mask_erosion_size, 'String'));
handles_b.settings.FOV               = str2num(get(handles.edit_interpolation_fov, 'String'));
handles_b.settings.NPixel            = str2num(get(handles.edit_interpolation_npixels, 'String'));
handles_b.settings.HardLimitsMax     = str2num(get(handles.edit_ampl_upper_boundry, 'String'));
handles_b.settings.HardLimitsMin     = str2num(get(handles.edit_ampl_lower_boundry, 'String'));
% set(handles_b.popupmenu_optimizationMethod, 'Value', 1); % To test

guidata(hObject_b, handles_b);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
handles = guidata(hObject);
currChar = get(handles.figure1,'CurrentCharacter');
if isequal(currChar,char(13)) %char(13) == enter key
    pushbutton_ok_Callback(hObject, eventdata, handles);
end
