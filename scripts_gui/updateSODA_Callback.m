function updateSODA_Callback(hObject, isRawdata, checkChange)
    %
    % isRawdata == 0 -> read from .ini file
    % isRawdata == 1 -> read from a matlab file
    %
    if nargin < 3
        checkChange = false;
    end
    
    handles = guidata(hObject);    
    setStatus(hObject, 1);
    if isRawdata == 1
        [FileName,PathName,FilterIndex] = uigetfile(fullfile(handles.settings.path2data,'*.dat'));
        if FilterIndex == 0
            setStatus(hObject, 0);
            return;
        end
        twix = mapVBVD(fullfile(PathName, FileName));
        if iscell(twix) % for VD/VE
            twix = twix{2};
        end
        mrprot = twix.hdr;
    elseif isRawdata == 2
        mrprot = handles.B0_OBJ.mrprot;
    elseif isRawdata == 0
        if exist(handles.settings.path2soda, 'dir') ~= 7
            handles.settings.path2soda = uigetdir;
            if handles.settings.path2soda == 0
                setStatus(hObject, 0);
                return;
            end
        end
    end
    
    %%
    boxNames = {'SLCVOL', 'ADJVOL'}; % Yellow, Green and Blue windows in the scanner
    if strcmp(get(handles.toolbar_select_group_shimming, 'State'), 'on') 
        boxNames{3} = 'NAVVOL';
    end
    % check if there is any change in SODA
    if checkChange && ~islogical(isRawdata) && handles.flags.SODALoaded
        while strcmp(get(handles.toolbar_check_soda_change, 'State'), 'on')
%             msg = cell(numel(boxNames), 1);
%             res = {'Change found', 'No Change'};
%             for i=1:numel(boxNames)
%                 box_new = SHIMBOX_OBJ(fullfile(handles.settings.path2soda, [boxNames{i}, '.ini']));
%                 box_old = handles.TARGET.(boxNames{i});
%                 temp = isequal(box_new.Position, box_old.Position);
%                 msg{i} = [boxNames{i}, ' = ' res{temp + 1}];
%             end
%             msgbox(msg, 'Results');
            temp = 0;
            for i=1:numel(boxNames)
                box_new = SHIMBOX_OBJ(fullfile(handles.settings.path2soda, [boxNames{i}, '.ini']));
                box_old = handles.TARGET.(boxNames{i});
                temp = temp+isequal(box_new.Position, box_old.Position);
            end
            if temp == numel(boxNames)
                cprintf('green', 'No Change Found\n');
            else
                cprintf('red', 'Warning, Change Found\n');
            end
            pause(1);
        end
        setStatus(hObject, 0);         
        return;
    end
    % read SODA
    colorInd = [0.5, 0.25, 0.75];
    handles.DISPLAY.Overlay = zeros(handles.settings.NPixel, handles.settings.NPixel, handles.settings.NPixel); 
    
    for i=1:numel(boxNames)
        if isRawdata ~= 0 % check if SODA source is from ini file
            handles.TARGET.(boxNames{i}) = SHIMBOX_OBJ(mrprot, boxNames{i});
        else
            handles.TARGET.(boxNames{i}) = SHIMBOX_OBJ(fullfile(handles.settings.path2soda, [boxNames{i}, '.ini']));       
        end
        if handles.TARGET.(boxNames{i}).NSlices == 0
            continue;
        end
        
        ShimBox = handles.TARGET.(boxNames{i});        
        coor_x = []; coor_y = []; coor_z = []; 
        for cSlc=1:ShimBox.NSlices % # of slices or groups or adjustment volume
            coor_x = cat(1, coor_x, ShimBox.Coords{cSlc}.SAG(:)); % I have to make a single vector for coordinate, otherwise "interpn" must be called per slice which make the whole calculation too slow
            coor_y = cat(1, coor_y, ShimBox.Coords{cSlc}.COR(:));
            coor_z = cat(1, coor_z, ShimBox.Coords{cSlc}.TRA(:));                      
        end
                
        MaxSizePixel = max(ShimBox.FOV(:)./ShimBox.NPixel(:))/2; 
        tempOverlay = InterpolateToStandardFoVMex(ones(size(coor_x)),coor_x, coor_y, coor_z, handles.settings.FOV, handles.settings.NPixel, MaxSizePixel); 
        tempOverlay(tempOverlay>1) = 1; % in case of overlap between slices
        
        tempOverlay = abs(tempOverlay - circshift(tempOverlay,[1 1 1])); % Fastest method to calculate 3d gradient for a cube while keeps the size
        handles.DISPLAY.Overlay(tempOverlay>0) = colorInd(i);
    end
    
    handles.DISPLAY.Alpha = handles.DISPLAY.Overlay*0; 
    handles.DISPLAY.Alpha(handles.DISPLAY.Overlay>0) = 0.7;
   
    handles.flags.SODALoaded = true;
    handles.flags.isCalcShimDone = false;   
    % Update handles
    guidata(hObject, handles);    
    updateShimInfo(hObject);  
    UpdateShimmedField(hObject); 
    UpdateB0Map(hObject);
    UpdateMainImages(hObject);
    setStatus(hObject, 0); 
end
