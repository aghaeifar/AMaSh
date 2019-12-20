function shim_scope(hObject, option, active)

handles = guidata(hObject);

switch option
    case handles.settings.shim_scope.options{1} % global shimming
        if active
            handles.settings.shim_scope.active = handles.settings.shim_scope.options{1};
            guidata(hObject, handles);  % have to save 'handles.settings.shim_scope.active' since callbacks are invoked when 'State' changes
            set(handles.toolbar_select_slice_shimming, 'State', 'Off'); % this cause to call 'toolbar_select_slice_shimming_OffCallback'
            set(handles.toolbar_select_group_shimming, 'State', 'Off');            
        elseif isequal(handles.settings.shim_scope.active, handles.settings.shim_scope.options{1})
            set(handles.toolbar_select_volume_shimming, 'State', 'On');
        end
    case handles.settings.shim_scope.options{2} % slice-wise shimming
        if active
            handles.settings.shim_scope.active = handles.settings.shim_scope.options{2};
            guidata(hObject, handles);
            set(handles.toolbar_select_volume_shimming, 'State', 'Off');
            set(handles.toolbar_select_group_shimming, 'State', 'Off');
            
        elseif isequal(handles.settings.shim_scope.active, handles.settings.shim_scope.options{2})
            set(handles.toolbar_select_slice_shimming, 'State', 'On');
        end   
    case handles.settings.shim_scope.options{3} % group-wise shimming
        if active
            handles.settings.shim_scope.active = handles.settings.shim_scope.options{3}; 
            guidata(hObject, handles);
            set(handles.toolbar_select_volume_shimming, 'State', 'Off');
            set(handles.toolbar_select_slice_shimming, 'State', 'Off');
            
        elseif isequal(handles.settings.shim_scope.active, handles.settings.shim_scope.options{3})
            set(handles.toolbar_select_group_shimming, 'State', 'On');
        end
end

guidata(hObject, handles);