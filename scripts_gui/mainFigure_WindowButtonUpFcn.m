% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function mainFigure_WindowButtonUpFcn(hObject, eventdata, handles)

    handles = guidata(hObject);
    handles.updateSliceShow = false;
    guidata(hObject, handles);