function mainFigure_WindowButtonDownFcn(hObject, eventdata, handles2)
    
    handles = guidata(hObject); 
    if ~strcmp(get(handles.mainFigure,'SelectionType'), 'normal') % only accept mouse left click
        return
    end

    handles.updateSliceShow = true;
    guidata(hObject, handles);
    mainFigure_WindowButtonMotionFcn(hObject, eventdata, handles2);
    
    