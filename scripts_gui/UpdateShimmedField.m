function UpdateShimmedField(hObject)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);   
    try
        %% Convert to truecolor
        cmap_size = 1024;
        ConversionMap = bipolar(cmap_size);
        ConversionMap = [0 0 0;ConversionMap];
        SAGMag  = handles.TARGET.STD_B0_SHIMMED(handles.cSAGSlice,:,:);
        CORMag  = handles.TARGET.STD_B0_SHIMMED(:,handles.cCORSlice,:);
        TRAMag  = handles.TARGET.STD_B0_SHIMMED(:,:,handles.cTRASlice);       

        %% Scale to range from 0..1
        B0Max = str2num(get(handles.editB0Max,'String'));
        B0Min = -B0Max;
        SAGMag  = (SAGMag(:)-B0Min)/(B0Max-B0Min);
        CORMag  = (CORMag(:)-B0Min)/(B0Max-B0Min);
        TRAMag  = (TRAMag(:)-B0Min)/(B0Max-B0Min);

        %% Convert
        handles.DISPLAY.ShimmedMap   =   zeros([handles.settings.NPixel, handles.settings.NPixel, handles.settings.NPixel,3]); 
        Index = round(SAGMag *(length(ConversionMap)-1))+1;
        Index(Index==0) = 2;
        Index(isnan(Index)) = 1;
        Index(Index>cmap_size) = cmap_size;
        Index(Index<0) = 2;
        handles.DISPLAY.ShimmedMap(handles.cSAGSlice,:,:,:) = reshape(ConversionMap(Index(:),:),[1, handles.settings.NPixel ,handles.settings.NPixel, 3]);

        Index = round(CORMag *(length(ConversionMap)-1))+1;
        Index(Index==0) = 2;
        Index(isnan(Index)) = 1;
        Index(Index>cmap_size) = cmap_size;
        Index(Index<0) = 2;
        handles.DISPLAY.ShimmedMap(:,handles.cCORSlice,:,:) = reshape(ConversionMap(Index(:),:),[handles.settings.NPixel, 1 ,handles.settings.NPixel, 3]);

        Index = round(TRAMag *(length(ConversionMap)-1))+1;
        Index(Index==0) = 2;
        Index(isnan(Index)) = 1;
        Index(Index>cmap_size) = cmap_size;
        Index(Index<0) = 2;
        handles.DISPLAY.ShimmedMap(:,:,handles.cTRASlice,:,:) = reshape(ConversionMap(Index(:),:),[handles.settings.NPixel ,handles.settings.NPixel, 1, 3]);

        %% Update handles structure
        guidata(hObject, handles);
    end