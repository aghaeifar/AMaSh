function UpdateB0Map(hObject)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);  
    try
        %% Convert to truecolor
        cmap_size = 1024;
        ConversionMap = bipolar(cmap_size);
        ConversionMap = [0 0 0;ConversionMap];    

        %% Masked or not?
        B0SAG  = (handles.B0_OBJ.STD_MAP_B0(handles.cSAGSlice,:,:));
        B0COR  = (handles.B0_OBJ.STD_MAP_B0(:,handles.cCORSlice,:));
        B0TRA  = (handles.B0_OBJ.STD_MAP_B0(:,:,handles.cTRASlice));
        if handles.flags.enMask
            B0SAG  = B0SAG.*handles.B0_OBJ.STD_Mask(handles.cSAGSlice,:,:);
            B0COR  = B0COR.*handles.B0_OBJ.STD_Mask(:,handles.cCORSlice,:);
            B0TRA  = B0TRA.*handles.B0_OBJ.STD_Mask(:,:,handles.cTRASlice);
        end
        %% Scale to range from 0..1
        B0Max = str2num(get(handles.editB0Max,'String'));
        B0Min = -B0Max;
        B0SAG   = (B0SAG-B0Min)/(B0Max - B0Min );
        B0COR   = (B0COR-B0Min)/(B0Max - B0Min );
        B0TRA   = (B0TRA-B0Min)/(B0Max - B0Min );

        %% Convert
        handles.DISPLAY.B0Map   =   zeros([handles.settings.NPixel, handles.settings.NPixel, handles.settings.NPixel, 3]);

        Index = round(B0SAG *(length(ConversionMap)-1))+1;
        Index(Index<=0) = 2;
        Index(isnan(Index)) = 1;
        Index(Index>cmap_size) = cmap_size;
        handles.DISPLAY.B0Map(handles.cSAGSlice,:,:,:) = reshape(ConversionMap(Index(:),:),[1,handles.settings.NPixel ,handles.settings.NPixel,3]);

        Index = round(B0COR *(length(ConversionMap)-1))+1;
        Index(Index<=0) = 2;
        Index(isnan(Index)) = 1;
        Index(Index>cmap_size) = cmap_size;
        handles.DISPLAY.B0Map(:,handles.cCORSlice,:,:) = reshape(ConversionMap(Index(:),:),[handles.settings.NPixel,1 ,handles.settings.NPixel,3]);

        Index = round(B0TRA *(length(ConversionMap)-1))+1;
        Index(Index<=0) = 2;
        Index(isnan(Index)) = 1;
        Index(Index>cmap_size) = cmap_size;
        handles.DISPLAY.B0Map(:,:,handles.cTRASlice,:) = reshape(ConversionMap(Index(:),:),[handles.settings.NPixel ,handles.settings.NPixel,1,3]);
    
        %% Update handles structure
        guidata(hObject, handles);
    end