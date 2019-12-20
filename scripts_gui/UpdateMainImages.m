function UpdateMainImages(hObject)
    % This function displays images in patient coordinate system   
    handles = guidata(hObject); 

    %% Colormap [Black, Green, Yellow, Blue White]
    Colormap = [0,0,0; 0,0.5,0; 0.8 0.8 0; 0.2,1,1; 1 1 1];
    % We need these lines, to update sliders and crosshair pos when UpdateMainImages called individually
    set(handles.editCrosshairPos, 'String', sprintf('%d %d %d', handles.cSAGSlice, handles.cCORSlice, handles.cTRASlice));
    
    %% Update standard image and shim basis sets if needed
    disp_img = handles.DISPLAY.Image;
    if handles.flags.enMask
        disp_img = handles.DISPLAY.Image .* repmat(handles.DISPLAY.Mask, [1 1 1 3]); 
    end
    
    % Sagital View
    cla(handles.axesSAG,'reset')   
    imagesc(permute(squeeze(disp_img(handles.cSAGSlice,:,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesSAG);
    hold(handles.axesSAG,'on')
    I=imagesc(permute(squeeze(handles.DISPLAY.Overlay(handles.cSAGSlice,:,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesSAG);
    colormap(handles.axesSAG,Colormap)
    set(handles.axesSAG,'CLim',[0 1], 'xtick',[],'ytick',[])
    set(I, 'AlphaData', permute(squeeze(handles.DISPLAY.Alpha(handles.cSAGSlice,:,end:-1:1,:)),[2 1 3]));
    
    axes(handles.axesSAG);
    line([0 handles.settings.NPixel ],[handles.settings.NPixel-handles.cTRASlice handles.settings.NPixel-handles.cTRASlice],'Color',[1 1 0],'LineWidth',1)
    line([handles.cCORSlice handles.cCORSlice],[0 handles.settings.NPixel ],'Color',[1 1 0],'LineWidth',1)
    
    % Coronal View
    cla(handles.axesCOR,'reset')
    imagesc(permute(squeeze(disp_img(:,handles.cCORSlice,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesCOR);
    hold(handles.axesCOR,'on')
    I2=imagesc(permute(squeeze(handles.DISPLAY.Overlay(:,handles.cCORSlice,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesCOR);
    colormap(handles.axesCOR,Colormap)
    set(handles.axesCOR,'CLim',[0 1],'xtick',[],'ytick',[])
    set(I2, 'AlphaData', permute(squeeze(handles.DISPLAY.Alpha(:,handles.cCORSlice,end:-1:1,:)),[2 1 3]));
    
    axes(handles.axesCOR);
    line([0 handles.settings.NPixel ],[handles.settings.NPixel-handles.cTRASlice handles.settings.NPixel-handles.cTRASlice],'Color',[1 1 0],'LineWidth',1)
    line([handles.cSAGSlice handles.cSAGSlice],[0 handles.settings.NPixel ],'Color',[1 1 0],'LineWidth',1)
    
    % Axial View
    cla(handles.axesTRA,'reset')
    imagesc(permute(squeeze(disp_img(:,:,handles.cTRASlice,:)),[2 1 3]), 'Parent', handles.axesTRA);
    hold(handles.axesTRA,'on')
    I3=imagesc(permute(squeeze(handles.DISPLAY.Overlay(:,:,handles.cTRASlice,:)),[2 1 3]), 'Parent', handles.axesTRA);
    colormap(handles.axesTRA,Colormap)
    set(handles.axesTRA,'CLim',[0 1],'xtick',[],'ytick',[])
    set(I3, 'AlphaData', permute(squeeze(handles.DISPLAY.Alpha(:,:,handles.cTRASlice,:)),[2 1 3]));
    
    axes(handles.axesTRA);
    line([0 handles.settings.NPixel ],[handles.cCORSlice handles.cCORSlice],'Color',[1 1 0],'LineWidth',1)
    line([handles.cSAGSlice handles.cSAGSlice],[0 handles.settings.NPixel ],'Color',[1 1 0],'LineWidth',1)        
    
    %% Always update images for B0 Map 
    
%     cla(handles.axesSAG2,'reset')
    imagesc(abs(permute(squeeze(handles.DISPLAY.B0Map(handles.cSAGSlice,:,end:-1:1,:)),[2 1 3])), 'Parent', handles.axesSAG2);
    set(handles.axesSAG2,'Visible','off');
    
%     cla(handles.axesCOR2,'reset')
    imagesc(permute(squeeze(handles.DISPLAY.B0Map(:,handles.cCORSlice,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesCOR2);
    set(handles.axesCOR2,'Visible','off');
    
%     cla(handles.axesTRA2,'reset')
    imagesc(permute(squeeze(handles.DISPLAY.B0Map(:,:,handles.cTRASlice,:)),[2 1 3]), 'Parent', handles.axesTRA2);
    set(handles.axesTRA2,'Visible','off');
    
    %% Update images for the shimmed brain   
    
%     cla(handles.axesSAG4,'reset')
    imagesc(abs(permute(squeeze(handles.DISPLAY.ShimmedMap(handles.cSAGSlice,:,end:-1:1,:)),[2 1 3])), 'Parent', handles.axesSAG4);
    set(handles.axesSAG4,'Visible','off');
    
%     cla(handles.axesCOR4,'reset')
    imagesc(permute(squeeze(handles.DISPLAY.ShimmedMap(:,handles.cCORSlice,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesCOR4);
    set(handles.axesCOR4,'Visible','off');
    
%     cla(handles.axesTRA4,'reset')
    imagesc(permute(squeeze(handles.DISPLAY.ShimmedMap(:,:,handles.cTRASlice,:)),[2 1 3]), 'Parent', handles.axesTRA4);
    set(handles.axesTRA4,'Visible','off');
    
    %% Update handles structure
    guidata(hObject, handles);
  
