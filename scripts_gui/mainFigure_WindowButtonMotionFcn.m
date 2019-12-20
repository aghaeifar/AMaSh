function mainFigure_WindowButtonMotionFcn(hObject, e, h)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);
    if ~isfield(handles, 'B0_OBJ')
        return;
    end

    cTRASlice = handles.cTRASlice;
    cCORSlice = handles.cCORSlice;
    cSAGSlice = handles.cSAGSlice;
    xrange = handles.settings.NPixel;
    yrange = handles.settings.NPixel;
%     xrange = range(get(handles.axesSAG, 'Xlim'));
%     yrange = range(get(handles.axesSAG, 'Ylim'));
     %% SAG
    mousePoint = get(handles.axesSAG, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cTRASlice = yrange-mouseY+1;
        cCORSlice = mouseX;                
    end
    
    mousePoint = get(handles.axesSAG2, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cTRASlice = yrange-mouseY+1;
        cCORSlice = mouseX;
    end
    
    mousePoint = round(get(handles.axesSAG4, 'CurrentPoint'));
    mouseX = mousePoint(1,1);
    mouseY = mousePoint(1,2);
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cTRASlice = yrange-mouseY+1;
        cCORSlice = mouseX;   
    end
    
    %% COR
    mousePoint = get(handles.axesCOR, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cTRASlice = yrange-mouseY+1;
        cSAGSlice = mouseX;
    end
    
    mousePoint = get(handles.axesCOR2, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cTRASlice = yrange-mouseY+1;
        cSAGSlice = mouseX;
    end
    
    mousePoint = get(handles.axesCOR4, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cTRASlice = yrange-mouseY+1;
        cSAGSlice = mouseX;
    end
    
    %% TRA
    mousePoint = get(handles.axesTRA, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cCORSlice = mouseY;
        cSAGSlice = mouseX;
    end  
    
    mousePoint = get(handles.axesTRA2, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cCORSlice = mouseY;
        cSAGSlice = mouseX;
    end 
    
    mousePoint = get(handles.axesTRA4, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        cCORSlice = mouseY;
        cSAGSlice = mouseX;
    end 
    
    %     if ~isequal(cSlice, [handles.cTRASlice, handles.cCORSlice, handles.cTRASlice]) % A click in valid areas occured
    set(handles.editCrosshairPos, 'string', sprintf('%d %d %d', cSAGSlice, cCORSlice, cTRASlice));
    if  handles.flags.isDataLoaded
        b0_value = num2str(handles.B0_OBJ.STD_MAP_B0(cSAGSlice, cCORSlice, cTRASlice), '%.0f');
        b0shimmed_value = [];
        if handles.flags.isCalcShimDone
            b0shimmed_value = [' / ' num2str(handles.TARGET.STD_B0_SHIMMED(cSAGSlice, cCORSlice, cTRASlice),'%.0f')];
        end
        set(handles.textCrosshairValue, 'String', ['B0 = ' b0_value b0shimmed_value]);
    end
    guidata(hObject, handles);
    handles = guidata(hObject);
    
    if isfield(handles, 'updateSliceShow') && handles.updateSliceShow
        handles.cTRASlice = cTRASlice;
        handles.cCORSlice = cCORSlice;
        handles.cSAGSlice = cSAGSlice;        
        guidata(hObject, handles);        
        UpdateShimmedField(hObject)     % Always must be called before UpdateMainImages
        UpdateB0Map(hObject);
        UpdateMainImages(hObject);
    end
    guidata(hObject, handles);
           

%     end