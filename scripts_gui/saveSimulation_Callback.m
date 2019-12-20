% ------------------------------------------------------------------
% ---------------         Save data to mat file        -------------    

function saveSimulation_Callback(hObject)

handles = guidata(hObject);
if ~handles.flags.SODALoaded
    disp(['SODA is loaded=' num2str(handles.flags.SODALoaded)]);
    return;
end

setStatus(hObject, 1);
TARGET          = handles.TARGET;
TARGET.MAP_B0   = TARGET.MAP_B0;
TARGET.Mask     = TARGET.Mask;
TARGET.twix     = handles.B0_OBJ.twix;
TARGET.settings = handles.settings;
if handles.flags.isCalcShimDone
    TARGET.B0_SHIMMED = TARGET.B0_SHIMMED;
end
if handles.flags.isShimMapsLoaded
    TARGET = rmfield(TARGET, 'MAP_Shims');
end
fileNames = {'MC_Global_Sim', 'MC_Global_Meas','..........................', ...
    'MC_Dynamic_Sim', 'MC_Dynamic_Meas', 'MC_Dynamic_Meas_FreqAdj',...
    'MC_SH01_Dynamic_Sim', 'MC_SH01_Dynamic_Meas', 'MC_SH01_Dynamic_Meas_FreqAdj', '..........................',...
    'MC_Dynamic_16_2mm_Sim', 'MC_Dynamic_16_2mm_Meas', 'MC_SH01_Dynamic_16_2mm_Sim', 'MC_SH01_Dynamic_16_2mm_Meas','..........................',...
    'MC_Dynamic_16_15mm_Sim', 'MC_Dynamic_16_15mm_Meas', 'MC_SH01_Dynamic_16_15mm_Sim', 'MC_SH01_Dynamic_16_15mm_Meas','..........................',...
    'MC_CSI_Sim', 'MC_CSI_Meas', 'MC_SVS_Sim', 'MC_SVS_Meas'};
% [Selection,ok] = listdlg('PromptString','Select a file Name:', 'SelectionMode','single','ListString', fileNames, 'ListSize', [210 370]);
[pathstr,name,~] = fileparts(siemensRawRenamer(handles.B0_OBJ.filename{1}));
% if ok && ~ismember(Selection, [3 10 15 20])
%     name =  fileNames{Selection};
% end
[FileName, PathName, FilterIndex] = uiputfile(fullfile(pathstr, ['TARGET_' name '.mat']),'Save file name');
if FilterIndex
    save(fullfile(PathName,FileName), 'TARGET');
    %         [~,name,ext] = fileparts(FileName);
    %         SODA = TARGET.SODA;
    %         ADJVOL = TARGET.ADJVOL;
    %         save(fullfile(PathName,[name, 'SODA', ext]), '_SODA', 'ADJVOL');
end
setStatus(hObject, 0);