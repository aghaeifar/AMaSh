% ------------------------------------------------------------------
% -------------------    Busy or Ready status     ------------------

function setStatus(hObject, stat)
    handles = guidata(hObject);
    if stat == 0
        set(handles.text_readyBusy, 'String', 'Ready', 'ForegroundColor', [0 0 0]);
    elseif stat == 1
        set(handles.text_readyBusy, 'String', 'Busy', 'ForegroundColor', 'r');
    end
    drawnow;
    guidata(hObject, handles);