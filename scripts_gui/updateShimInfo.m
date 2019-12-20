% ------------------------------------------------------------------
% -------------------    Display Shimming information     ----------
function updateShimInfo(hObject)
    handles = guidata(hObject);

    protocol = [];
%     oldCurr = [];
    Scanner_Shim_Setup = [];
    if handles.flags.isDataLoaded
        type1 = {'0x1', '0x2', '0x4'};
        type2 = {'Ascending', 'Desceinding', 'Interleaved'};
        dataSeries = type2{ismember(type1, handles.B0_OBJ.twix.hdr.MeasYaps.sSliceArray.ucMode)};
        type1 = {'0x2', '0x4'};
        type2 = {'2D', '3D'};
        d2d3d = type2{ismember(type1, handles.B0_OBJ.twix.hdr.MeasYaps.sKSpace.ucDimension)};
        dGroup = handles.B0_OBJ.twix.hdr.MeasYaps.sGroupArray.lSize;
        dSlabSlc = handles.B0_OBJ.twix.image.dataSize(5);
        dPar = handles.B0_OBJ.twix.image.dataSize(4);
        dTE  = cell2mat(handles.B0_OBJ.twix.hdr.MeasYaps.alTE(1:handles.B0_OBJ.twix.hdr.MeasYaps.lContrasts))/1000;
        dTR  = handles.B0_OBJ.twix.hdr.MeasYaps.alTR{1}/1000;
        protocol = ['Series:' dataSeries ', ' d2d3d ', Group:' num2str(dGroup) ', Slab/Slice:' num2str(dSlabSlc) ', Partition:' num2str(dPar) ', TE/TR: [' num2str(dTE, '%.1f ') ']/' num2str(dTR)];

        SH = handles.B0_OBJ.SHIM.SH;
        Scanner_Shim_Setup = sprintf('%d | %s | %s', SH(1), num2str(SH(2:4),'  %.2f'), num2str(SH(5:9),' %.2f'));
        
%         for i=1:size(handles.B0_OBJ.SHIM.Dynamic,1)
%             shim_sh0 = num2str(round(handles.B0_OBJ.SHIM.Dynamic(i,1)), '%d');
%             shim_sh1 = num2str(round(handles.B0_OBJ.SHIM.Dynamic(i,2:4)), ' %d');
%             shim_mc  = num2str(handles.B0_OBJ.SHIM.Dynamic(i,5:end), ' %.1f');
%             oldCurr = sprintf('%s%s %s %s\n', oldCurr, shim_sh0, shim_sh1, shim_mc);
%         end
        data_sh0    = num2cell(round(handles.B0_OBJ.SHIM.Dynamic(:,1)));     % Freq
        data_sh1    = num2cell(round(handles.B0_OBJ.SHIM.Dynamic(:,2:4)));   % Gradients
        data_mc     = num2cell(handles.B0_OBJ.SHIM.Dynamic(:,5:end));        % Multi-Coil
        
        data_sh0    = cellfun(@(x) num2str(x,'%d')  , data_sh0  , 'UniformOutput', false); % convert all numbers to string to display in the table
        data_sh1    = cellfun(@(x) num2str(x,'%d')  , data_sh1  , 'UniformOutput', false);
        data_mc     = cellfun(@(x) num2str(x,'%.2f'), data_mc   , 'UniformOutput', false);
        data        = [data_sh0, data_sh1, data_mc]; 
        set(handles.uitable_oldShim,'Data', data);
        tLen        = num2cell(35*ones(1,4 + size(data_mc,2)));
        set(handles.uitable_oldShim,'ColumnWidth', tLen);
        tColName    = {'Freq'; 'X'; 'Y'; 'Z'; strcat('Ch', num2str((1:size(data_mc,2))'))};
        set(handles.uitable_oldShim,'ColumnName', tColName); 
        
    end
    
    sodaInfo = [];
    if handles.flags.SODALoaded
        sodaInfo = sprintf('Series:?, %dD, Slab/Slice:%d, Partition:?, Tickness:%.1fmm', handles.TARGET.SLCVOL.Dim, handles.TARGET.SLCVOL.NSlices, handles.TARGET.SLCVOL.Thickness(1));
    end
    
    temp = sprintf('[Rawdata Info]\n%s\n[Next Protocol]\n%s\n[Scanner Shim info]\n%s', protocol, sodaInfo, Scanner_Shim_Setup);
    set(handles.edit_shimInfo, 'String', temp);
 
    %% 
    if handles.flags.isCalcShimDone
        handles.TARGET.NewShimSettings = zeros(size(handles.TARGET.DeltaShimSettings));
        for cSlc=1:size(handles.TARGET.DeltaShimSettings, 1)
%             try
%                 handles.TARGET.NewShimSettings(cSlc,2:end) = handles.TARGET.DeltaShimSettings(cSlc,2:end) + handles.B0_OBJ.SHIM.Dynamic(abs(handles.B0_OBJ.SODA.Dim-3)*cSlc + (handles.B0_OBJ.SODA.Dim-2),2:size(handles.SHIM.STD_B0, 4));
%                 handles.TARGET.NewShimSettings(cSlc,1) = round(-handles.TARGET.DeltaShimSettings(cSlc,1) + handles.B0_OBJ.SHIM.Dynamic(abs(handles.B0_OBJ.SODA.Dim-3)*cSlc + (handles.B0_OBJ.SODA.Dim-2),1));     
%             catch
                handles.TARGET.NewShimSettings(cSlc,:) = handles.TARGET.DeltaShimSettings(cSlc,:);
%                 handles.TARGET.NewShimSettings(cSlc,1) = round(-handles.TARGET.DeltaShimSettings(cSlc,1));
%             end
        end
        if max(handles.settings.shim_channels) <= numel(handles.SHIM.Name) % selected channels have a name
            selected_ch_name = handles.SHIM.Name(handles.settings.shim_channels);
        end 
        
        data_all  = num2cell(handles.TARGET.NewShimSettings(:,handles.settings.shim_channels)); % all coef, SH and MC

        ind =  2*handles.flags.isDuplicateSession + 1;% 1 when is not duplicated session and 3 when it is
        data_mean = num2cell(transpose(round(handles.TARGET.MEAN([ind,2],:))));   % Offset
        data_std  = num2cell(transpose(handles.TARGET.STD([ind,2],:)));           % Deviation 
        
        data_all  = cellfun(@custom_num2str, data_all , 'UniformOutput', false);
        data_mean = cellfun(@custom_num2str, data_mean, 'UniformOutput', false);
        data_std  = cellfun(@custom_num2str, data_std , 'UniformOutput', false);
        
        data      = [data_mean, data_std, data_all]; 
        set(handles.uitable_coef,'Data', data);      
        tLen      = num2cell(35*ones(1,4 + numel(selected_ch_name))); % 4 is for 'Mean'; ''; 'STD'; ''; see below
        set(handles.uitable_coef,'ColumnWidth', tLen);
        tColName  = cat(1, {'Mean'; ''; 'STD'; ''}, selected_ch_name(:));
        set(handles.uitable_coef,'ColumnName', tColName);                      
    end

    guidata(hObject, handles);
end    

function x = custom_num2str(x)
    if abs(x)>100
        x= sprintf('%d',round(x));
    elseif abs(x)>10
        x= sprintf('%.1f',x);
    else
        x= sprintf('%.2f',x);
    end
end