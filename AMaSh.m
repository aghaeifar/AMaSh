function varargout = AMaSh(varargin)

gui_Singleton = 0;  % Allow to duplicate session
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AMaSh_OpeningFcn, ...
                   'gui_OutputFcn',  @AMaSh_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
   gui_mainfcn(gui_State, varargin{:});
end

function AMaSh_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for AMaSh
handles.output = hObject;

guidata(hObject, handles);
% Get correct handles (handles provided by the function arguments are not working)
handles = guidata(hObject);


% aa: Update Path
f_path = fileparts(mfilename('fullpath'));
addpath(fullfile(f_path, 'scripts_gui'));
addpath(fullfile(f_path, 'scripts'));
addpath(fullfile(f_path, 'scripts/reco3D'));
addpath(fullfile(f_path, 'scripts/reco3D/mapVBVD'));
addpath(fullfile(f_path, 'scripts/bet'));

% Flags
handles.flags.isShimMapsLoaded      = false;
handles.flags.isCalcShimDone        = false;
handles.flags.isDataLoaded          = false;
handles.flags.SODALoaded            = false;
handles.flags.is2D              	= false;
handles.flags.enMask                = false;
handles.flags.enUnwrap              = false;
handles.flags.isDuplicateSession    = false;
handles.flags.isShimBoxSODA         = false;
% Settings
handles.settings.shim_channels      = 5:36; %[5 36];
handles.settings.shim_algorithm_no  = 4;
handles.settings.unwrapping_method  = 2; % 1:None 2:Spatial 3:Temporal
handles.settings.adjacent_slices    = 1;
handles.settings.mask_algorithm_no  = 2; % 1:None 2:Bet 3:Threshold
handles.settings.mask_erosion_sz    = 2;
handles.settings.FOV                = 300;
handles.settings.NPixel             = 151;
handles.settings.HardLimitsMax      = [0.0 0.0];
handles.settings.HardLimitsMin      = [0.0 0.0];
handles.settings.reco_mode          = 1; % 1:sum of squares  2:adaptive combine, 
if isunix
    handles.settings.path2soda      = '/mnt/mrz9t/SODA_ADJ_DATA';
    handles.settings.path2data      = '/mnt/mrz9t/9T_Ali/Rawdata'; 
    handles.settings.path2save_coef = '/mnt/temp/amoco/shim';     
elseif ispc
    handles.settings.path2soda      = 'M:\Upload9T\USERS\Ali\SODA_ADJ_DATA';
    handles.settings.path2data      = 'M:\Upload9T\USERS\Ali';
    handles.settings.path2save_coef = 'M:\Upload9T\USERS\Ali';    
end

handles.settings.shim_scope.options = {'Global'; 'Slice-Wise'; 'Group-Wise'};
handles.settings.shim_scope.active  = handles.settings.shim_scope.options{1};
guidata(hObject, handles);
set(handles.toolbar_select_volume_shimming, 'State', 'On');

set(handles.editCrosshairPos, 'String', sprintf('%d %d %d', round((handles.settings.NPixel+1)/2)*ones(1,3)));

set(handles.uitable_coef,'Data',cell(0,0));
set(handles.uitable_oldShim,'Data',cell(0,0));

% Loc current slice
handles.cSAGSlice = (handles.settings.NPixel+1)/2;
handles.cCORSlice = (handles.settings.NPixel+1)/2;
handles.cTRASlice = (handles.settings.NPixel+1)/2;

% Image
handles.DISPLAY.Image       = zeros(handles.settings.NPixel, handles.settings.NPixel, handles.settings.NPixel,3,'single');
handles.DISPLAY.ShimmedMap  = zeros(handles.settings.NPixel, handles.settings.NPixel, handles.settings.NPixel,3,'single');
handles.DISPLAY.B0Map       = zeros(handles.settings.NPixel, handles.settings.NPixel, handles.settings.NPixel,3,'single');
handles.DISPLAY.Overlay     = zeros(handles.settings.NPixel, handles.settings.NPixel, handles.settings.NPixel,'single');
handles.DISPLAY.Alpha       = zeros(handles.settings.NPixel, handles.settings.NPixel, handles.settings.NPixel,'single');
handles.DISPLAY.Mask        = ones(handles.settings.NPixel,  handles.settings.NPixel, handles.settings.NPixel,'single');
handles.DISPLAY.Max = 1;
handles.DISPLAY.Min = 0;

% handles.TARGET.STD_B0_SHIMMED = ones(handles.settings.NPixel,handles.settings.NPixel,handles.settings.NPixel,'single');

% Black out images
axes(handles.axesSAG)
imagesc(zeros(handles.settings.NPixel, handles.settings.NPixel),[0 1]); axis off, axis xy;
colormap(gray)

axes(handles.axesCOR)
imagesc(zeros(handles.settings.NPixel, handles.settings.NPixel),[0 1]); axis off, axis xy;

axes(handles.axesTRA)
imagesc(zeros(handles.settings.NPixel, handles.settings.NPixel),[0 1]); axis off, axis xy;

axes(handles.axesSAG2)
imagesc(zeros(handles.settings.NPixel, handles.settings.NPixel),[0 1]); axis off, axis xy;

axes(handles.axesCOR2)
imagesc(zeros(handles.settings.NPixel, handles.settings.NPixel),[0 1]); axis off, axis xy;

axes(handles.axesTRA2)
imagesc(zeros(handles.settings.NPixel,handles.settings.NPixel),[0 1]); axis off, axis xy;

axes(handles.axesSAG4)
imagesc(zeros(handles.settings.NPixel, handles.settings.NPixel),[0 1]); axis off, axis xy;

axes(handles.axesCOR4)
imagesc(zeros(handles.settings.NPixel, handles.settings.NPixel),[0 1]); axis off, axis xy;

axes(handles.axesTRA4)
imagesc(zeros(handles.settings.NPixel, handles.settings.NPixel),[0 1]); axis off, axis xy;


if ~isempty(varargin)
%   assignin('base','b',varargin);   
    strList = fieldnames(varargin{1});
    % idea of recovering is inspired from http://www.mathworks.com/matlabcentral/fileexchange/24669
    def.checkbox = {'value'};
    def.edit = {'string'};
    def.text = {'string'};
    def.slider = {'value','max','min'};
    def.popupmenu = {'value','string'};
    for i=1:numel(strList)
       if ~ishandle(varargin{1}.(strList{i})) % dynamic field references instead of using "getfield"
           handles.(strList{i}) = varargin{1}.(strList{i});
       elseif strcmp(get(varargin{1}.(strList{i}), 'type'), 'uicontrol')
           style = get(varargin{1}.(strList{i}), 'style');
           if isfield(def, style)
               set(handles.(strList{i}), def.(style), get(varargin{1}.(strList{i}), def.(style)))
           end
       elseif strcmp(get(varargin{1}.(strList{i}), 'type'), 'uitoggletool')
           set(handles.(strList{i}), 'State', get(varargin{1}.(strList{i}), 'State'))
       end
    end
    guidata(hObject, handles);
    UpdateShimmedField(hObject)
    UpdateB0Map(hObject);
    UpdateMainImages(hObject);
    updateShimInfo(hObject);
end
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = AMaSh_OutputFcn(hObject, e, h) 
% Get correct handles (handles provided by the function arguments are not working)
handles = guidata(hObject);
% Get default command line output from handles structure
varargout{1} = handles.output;

    


% ------------------------------------------------------------------
% -------------------   Color Range       --------------------------
function editB0Max_Callback(hObject, e, h)
    UpdateB0Map(hObject);
    UpdateShimmedField(hObject)
    UpdateMainImages(hObject);
    
% ------------------------------------------------------------------
% -------------------   Phase Shift       --------------------------    
function edit_phase_shift_2pi_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    if handles.flags.isDataLoaded
        shift = str2double(get(handles.edit_phase_shift_2pi, 'String'));
        if isnan(shift)
            set(handles.edit_phase_shift_2pi, 'String', '0');
            return;
        end
        handles.B0_OBJ = unwrapB0(handles.B0_OBJ, shift);
        if handles.flags.enUnwrap
            handles.B0_OBJ.MAP_B0 = handles.B0_OBJ.MAP_B0_Unwrapped;
            handles.B0_OBJ.STD_MAP_B0 = handles.B0_OBJ.STD_MAP_B0_Unwrapped; 
        end
        guidata(hObject, handles);
        UpdateB0Map(hObject);
        UpdateMainImages(hObject);
    end
    
% --------------------------------------------------------------------
% -------------------             Menu            --------------------
% --------------------------------------------------------------------
function contexMenu1_Plot_Results_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    if handles.flags.isCalcShimDone   
        NSli        = handles.TARGET.ShimBox.NSlices;
        Reorder = 1:NSli;
%         if handles.flags.SODALoaded && handles.TARGET.SLCVOL.IsInterleaved % are slices ordered in a interleaved manner?             
%             Order = [~mod(NSli,2)+1:2:NSli, mod(NSli,2)+1:2:NSli]; % depends on number of slices, it is odd or even!
%             for i=1:NSli
%                 Reorder(i) = find(Order==i); 
%             end
%         end   
        sz = round(get(0,'ScreenSize')/4);
        figure('WindowStyle', 'modal', 'Name', 'Evaluate Shim Performance', 'Position', [sz(3), sz(4), 2*sz(3), 2*sz(4)]); % need modal to avoid "WindowButtonMotionFcn"
        
        subplot(211);  hold on;  ub = max(abs(handles.TARGET.MEAN(:)));     
        h1=plot(handles.TARGET.MEAN(1,Reorder), 'r'); 
        h2=plot(handles.TARGET.MEAN(2,Reorder), 'b'); 
        legend([h1, h2], 'Before Shim', 'Simulation'); ylim([-ub ub]);
        if handles.flags.isDuplicateSession
            h3=plot(handles.TARGET.MEAN(3,Reorder), 'g'); 
            legend([h1, h2, h3], 'Before Shim', 'Simulation', 'Measurment');
        end
        title('Mean')
        
        subplot(212); hold on; 
        h1=plot(handles.TARGET.STD(1,Reorder), 'r');
        h2=plot(handles.TARGET.STD(2,Reorder), 'b');
        legend([h1, h2], 'Before Shim', 'Simulation');
        if handles.flags.isDuplicateSession
            h3 = plot(handles.TARGET.STD(3,Reorder), 'g'); 
            legend([h1, h2, h3],'Before Shim', 'Simulation', 'Measurment');
        end        
        title('STD')
    end
    
% --------------------------------------------------------------------
function contexMenu1_Plot_Currents_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
if handles.flags.isCalcShimDone    
    h = figure;
    set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
    nSlc = size(handles.TARGET.NewShimSettings,1);
    for i=1:size(handles.TARGET.NewShimSettings,2)
        subplot(5,4,i);
        plot(1:nSlc,handles.TARGET.NewShimSettings(:,i), 'b-o', 'MarkerFaceColor',[153 255 51]/255);
        xlim([0 nSlc+1]); ylim([ handles.SHIM.HardLimitsMin(i), handles.SHIM.HardLimitsMax(i)]*1.5);        
        set(gca, 'box','off','XTick',[])
        line([0 nSlc+1], [handles.SHIM.HardLimitsMin(i), handles.SHIM.HardLimitsMin(i)], 'Color','red','LineStyle','--');
        line([0 nSlc+1], [handles.SHIM.HardLimitsMax(i), handles.SHIM.HardLimitsMax(i)], 'Color','red','LineStyle','--');
    end
end

% --------------------------------------------------------------------
function menuFile_OpenRawdata_Callback(hObject, eventdata, handles)
    recoB0Map(hObject)
    
% --------------------------------------------------------------------
function menuFile_OpenDICOM_Callback(hObject, eventdata, handles)
    msgbox('Under development...','Info','Help')
    
% --------------------------------------------------------------------
function menuFile_OpenMatlab_Callback(hObject, eventdata, handles)
    msgbox('Under development...','Info','Help')

% --------------------------------------------------------------------
function menuFile_OpenIce_Callback(hObject, eventdata, handles)
    msgbox('Under development...','Info','Help')
    
% --------------------------------------------------------------------
function menuFile_save_Callback(hObject, eventdata, handles)
    saveSimulation_Callback(hObject);
    
% --------------------------------------------------------------------
function menuFile_exit_Callback(hObject, eventdata, handles)
    delete(gcbf);
    
% --------------------------------------------------------------------
% --------------------------------------------------------------------
function menuTools_AdvancedUserAccess_Callback(hObject, eventdata, handles)
% Master pass
handles = guidata(hObject);
c           = clock;
master_pass = sprintf('%d',c(1:4));
prompt      = {'Enter Password:'};
dlg_title   = 'Advanced Options';
num_lines   = [1 45];
defaultans  = {'******'};

while 1
    answer = passcode;
    if ~isempty(answer) 
        if strcmp(answer, master_pass)
            set(handles.menuTools_MakeBasismaps,        'Enable', 'on');
            set(handles.menuTools_GenerateShimCoef,     'Enable', 'on');
            set(handles.menuTools_GetSODAfromSeparateRawdata, 'Enable', 'on');
            set(handles.menuTools_GetSODAfromCurrentRawdata,  'Enable', 'on');
            set(handles.menuTools_CheckSODAchange,   'Enable', 'on'); 
            set(handles.menuTools_updateShimInfo,       'Enable', 'on');
            set(handles.menuFile_OpenDICOM,             'Enable', 'on');
            set(handles.menuFile_OpenMatlab,            'Enable', 'on');
            set(handles.menuFile_OpenIce,               'Enable', 'on');           
            set(handles.toolbar_select_group_shimming,  'Enable', 'on');
            set(handles.toolbar_select_slice_shimming,  'Enable', 'on');
            set(handles.menuTools_AdvancedUserAccess,   'Enable', 'off');
            set(handles.toolbar_advanced_settings,   'Enable', 'off');

            guidata(hObject, handles);
            uiwait(msgbox('Welcome to advanced user area.'));
            break;
        else
            uiwait(msgbox('Wrong Password', 'Error','error'));
        end 
    else
        break;
    end
end
% --------------------------------------------------------------------
function menuTools_MakeBasismaps_Callback(hObject, eventdata, handles)
    [SHIM, data] = shimMapping;
% --------------------------------------------------------------------
function menuTools_GenerateShimCoef_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menuTools_updateShimInfo_Callback(hObject, eventdata, handles)
    updateShimInfo(hObject);
    
% --------------------------------------------------------------------
function menuTools_GetSODAfromSeparateRawdata_Callback(hObject, eventdata, handles)
    updateSODA_Callback(hObject, 1);
 
% --------------------------------------------------------------------
function menuTools_GetSODAfromCurrentRawdata_Callback(hObject, eventdata, handles)  
    updateSODA_Callback(hObject, 2);
    
% --------------------------------------------------------------------
function menuTools_CheckSODAchange_Callback(hObject, eventdata, handles)
    updateSODA_Callback(hObject, 0, 1);
    
% --------------------------------------------------------------------
% --------------------------------------------------------------------
function menuHelp_About_Callback(hObject, eventdata, handles)
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
msgbox({'Advacned Magnetic field Shimming (AMaSh)',...
        'Authors: Ali Aghaeifar & Christian Mirkes',...
        'https://github.com/Aghaeifar/AMaSh',...
        'For questions contact:',...
        'ali.aghaeifar@tuebingen.mpg.de'},...
        'About', CreateStruct);
% --------------------------------------------------------------------
% --------------------------------------------------------------------
function toolbar_open_ClickedCallback(hObject, eventdata, handles)
    recoB0Map(hObject);
    
% --------------------------------------------------------------------
function toolbar_save_ClickedCallback(hObject, eventdata, handles)
    saveSimulation_Callback(hObject);
    
% --------------------------------------------------------------------
function toolbar_duplicate_session_ClickedCallback(hObject, eventdata, handles)
    handles = guidata(hObject);
    handles.flags.isDuplicateSession = true;
    B0GUI(handles);
    
% --------------------------------------------------------------------
function toolbar_load_basis_Maps_ClickedCallback(hObject, eventdata, handles)
    load_basissets(hObject);
    
% --------------------------------------------------------------------
function toolbar_see_basis_Maps_ClickedCallback(hObject, eventdata, handles)
    handles = guidata(hObject);
    if handles.flags.isShimMapsLoaded == false
       h = msgbox('Basis-maps are not loaded!', 'Warning','Warn'); 
       set(h, 'WindowStyle', 'modal')
       return; 
    end
    BasisMapsDisp(handles.SHIM.STD_B0 .* repmat(handles.SHIM.STD_Mask,[1 1 1 handles.SHIM.NShimCha]), handles.DISPLAY.Mask);
    
% --------------------------------------------------------------------
function toolbar_creat_mask_ClickedCallback(hObject, eventdata, handles)
    calcMask_Callback(hObject);
    
% --------------------------------------------------------------------
function toolbar_update_soda_ClickedCallback(hObject, eventdata, handles)
    updateSODA_Callback(hObject, 0);
 
% --------------------------------------------------------------------
function toolbar_check_soda_ClickedCallback(hObject, eventdata, handles)
    

% --------------------------------------------------------------------
function toolbar_calculate_shim_ClickedCallback(hObject, eventdata, handles)
    calcShim_Callback(hObject);

% --------------------------------------------------------------------    
function toolbar_calculate_freq_ClickedCallback(hObject, eventdata, handles)
    calcFreqOffset_Callback(hObject);
    
% --------------------------------------------------------------------
function toolbar_evaluate_shim_ClickedCallback(hObject, eventdata, handles)  
    handles = guidata(hObject);
    if handles.flags.isShimMapsLoaded && handles.flags.isCalcShimDone && handles.flags.isDataLoaded 
        sz_mapShim = size(handles.SHIM.STD_B0(1:2:end, 1:2:end, 1:2:end, :));
        dataTemp = zeros([sz_mapShim(1:end-1), handles.TARGET.ShimBox.NSlices]);
        for i=1:handles.TARGET.ShimBox.NSlices
            coef = handles.TARGET.DeltaShimSettings(i, :)';
            mapShim_sum = repmat(permute(coef, [2:4 1]), [sz_mapShim(1:end-1) 1]).* handles.SHIM.STD_B0(1:2:end, 1:2:end, 1:2:end, :);
            mapShim_sum = squeeze(sum(mapShim_sum, numel(sz_mapShim)));
            dataTemp(:,:,:,i)  = mapShim_sum + handles.B0_OBJ.STD_MAP_B0(1:2:end, 1:2:end, 1:2:end) .* handles.B0_OBJ.STD_Mask(1:2:end,1:2:end,1:2:end);
        end
     BasisMapsDisp(dataTemp);
    end

% --------------------------------------------------------------------
function toolbar_write_shim_ClickedCallback(hObject, eventdata, handles)
    writeShimCoef_Callback(hObject);
    
% --------------------------------------------------------------------
function toolbar_settings_ClickedCallback(hObject, eventdata, handles)
    mySettings(hObject);
    
% --------------------------------------------------------------------
function toolbar_use_mask_OffCallback(hObject, eventdata, handles)
    handles = guidata(hObject);
    handles.flags.enMask  = false;
    guidata(hObject, handles)
    if handles.flags.isDataLoaded == false
       return; 
    end
    UpdateShimmedField(hObject)     % Always must be called before UpdateMainImages
    UpdateB0Map(hObject);
    UpdateMainImages(hObject);

% --------------------------------------------------------------------
function toolbar_use_mask_OnCallback(hObject, eventdata, handles)
    handles = guidata(hObject);
    handles.flags.enMask = true;
    guidata(hObject, handles)
    if handles.flags.isDataLoaded == false
       return; 
    end
    UpdateShimmedField(hObject)     % Always must be called before UpdateMainImages
    UpdateB0Map(hObject);
    UpdateMainImages(hObject);

% --------------------------------------------------------------------
function toolbar_unwrap_OffCallback(hObject, eventdata, handles)
    handles = guidata(hObject);
    if handles.flags.isDataLoaded == false
        return;
    end
    handles.flags.enUnwrap  = false;
    handles.B0_OBJ.MAP_B0 = handles.B0_OBJ.MAP_B0_Wrapped;
    handles.B0_OBJ.STD_MAP_B0 = handles.B0_OBJ.STD_MAP_B0_Wrapped; 
    guidata(hObject, handles)
    UpdateShimmedField(hObject);
    UpdateB0Map(hObject);
    UpdateMainImages(hObject);
  
% --------------------------------------------------------------------
function toolbar_unwrap_OnCallback(hObject, eventdata, handles)    
    handles = guidata(hObject);
    if handles.flags.isDataLoaded == false
        return;
    end
    handles.flags.enUnwrap  = true;
    handles.B0_OBJ.MAP_B0 = handles.B0_OBJ.MAP_B0_Unwrapped;
    handles.B0_OBJ.STD_MAP_B0 = handles.B0_OBJ.STD_MAP_B0_Unwrapped;
    guidata(hObject, handles)
    UpdateShimmedField(hObject);
    UpdateB0Map(hObject);
    UpdateMainImages(hObject); 

% --------------------------------------------------------------------
function toolbar_select_volume_shimming_OffCallback(hObject, eventdata, handles)
    shim_scope(hObject, handles.settings.shim_scope.options{1}, 0);

% --------------------------------------------------------------------
function toolbar_select_volume_shimming_OnCallback(hObject, eventdata, handles)
    shim_scope(hObject, handles.settings.shim_scope.options{1}, 1);

% --------------------------------------------------------------------
function toolbar_select_slice_shimming_OffCallback(hObject, eventdata, handles)
    shim_scope(hObject, handles.settings.shim_scope.options{2}, 0);
    
% --------------------------------------------------------------------
function toolbar_select_slice_shimming_OnCallback(hObject, eventdata, handles)  
    shim_scope(hObject, handles.settings.shim_scope.options{2}, 1);

% --------------------------------------------------------------------
function toolbar_select_group_shimming_OffCallback(hObject, eventdata, handles)
    shim_scope(hObject, handles.settings.shim_scope.options{3}, 0);
    
% --------------------------------------------------------------------
function toolbar_select_group_shimming_OnCallback(hObject, eventdata, handles)
    shim_scope(hObject, handles.settings.shim_scope.options{3}, 1);
      
% --------------------------------------------------------------------
function toolbar_check_soda_change_OnCallback(hObject, eventdata, handles)    
    updateSODA_Callback(hObject, 0, 1);
    
%%################################################################################################################  
function editCrosshairPos_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    pos = get(handles.editCrosshairPos, 'string');
    pos = round(sscanf(pos,'%f'));
    if numel(pos) ~= 3 
        pos = round(handles.settings.NPixel/2)*ones(1,3);
    end
    pos(pos>handles.settings.NPixel) = handles.settings.NPixel;
    pos(pos<1) = 1;
    set(handles.editCrosshairPos, 'string', sprintf('%d ', pos));
    handles.cTRASlice = pos(3);
    handles.cCORSlice = pos(2);
    handles.cSAGSlice = pos(1);    
    guidata(hObject, handles);   % Update handles  
    if handles.flags.isDataLoaded == false
        return;
    end    
    b0_value = num2str(handles.B0_OBJ.STD_MAP_B0(handles.cSAGSlice, handles.cCORSlice, handles.cTRASlice), '%.0f');
    set(handles.textCrosshairValue, 'String', b0_value);
    guidata(hObject, handles);   % Update handles 
    UpdateShimmedField(hObject)  % Always must be called before UpdateMainImages
    UpdateB0Map(hObject);        % Update B0 map    
    UpdateMainImages(hObject);   % Update main images    


% 

% hObject    handle to menuTools_GetSODAfromCurrentRawdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% hObject    handle to menuTools_CheckSODAchange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% hObject    handle to toolbar_is_soda_changed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% hObject    handle to toolbar_check_soda_change (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
