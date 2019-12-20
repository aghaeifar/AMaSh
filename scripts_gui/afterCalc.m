% ------------------------------------------------------------------
% -------------------  What to do after shimming calculations  -----
function afterCalc(hObject)
    handles = guidata(hObject);     
    
    coor_x = []; coor_y = []; coor_z = [];    
    temp_STD_B0_SHIMMED = [];
    if handles.flags.isShimBoxSODA2D   
        % remove slice over sampling and return to original SODA
        CoordsBckup = handles.TARGET.ShimBox.Coords;
        handles.TARGET.ShimBox.Thickness(:) = handles.TARGET.ShimBox.Thickness / (2*handles.settings.adjacent_slices+1);
        handles.TARGET.ShimBox.NPixelSlice(:) = handles.TARGET.ShimBox.NPixelSlice - 4*handles.settings.adjacent_slices;
        handles.TARGET.ShimBox = handles.TARGET.ShimBox.calcPixLoc;
    end
    for cSlc = 1:handles.TARGET.ShimBox.NSlices
        data_temp = handles.TARGET.B0_SHIMMED{cSlc}(:); % contains NAN for anything ourside mask (mask of brain + mask of basis-maps)
        if handles.flags.isShimBoxSODA2D  % Calculation of B0 deviation without slice oversampling when SODA is selected as shimBOX
            data_temp  = reshape(handles.TARGET.B0_SHIMMED{cSlc}, size(CoordsBckup{cSlc}.SAG));
            ind = ceil(size(data_temp, 3)/2); % we do this since we need just the real slice which will be measured, not the neighborhood
            data_temp  = data_temp(:,:,ind-1:ind+1);                                      
        end            
        handles.TARGET.MEAN(2,cSlc) = mean(data_temp(:), 'omitnan'); % voxels outside of the mask are NAN (mask of brain + mask of basis-maps)
        handles.TARGET.STD(2,cSlc)  = std( data_temp(:), 'omitnan');
        % I have to make a single vector for coordinate, otherwise "InterpolateToStandardFoVMex" must be called per slice which make the whole calculation too slow and also adjucent slices will overlap
        coor_x = cat(1, coor_x, handles.TARGET.ShimBox.Coords{cSlc}.SAG(:));
        coor_y = cat(1, coor_y, handles.TARGET.ShimBox.Coords{cSlc}.COR(:));
        coor_z = cat(1, coor_z, handles.TARGET.ShimBox.Coords{cSlc}.TRA(:));
        temp_STD_B0_SHIMMED = cat(1, temp_STD_B0_SHIMMED, data_temp(:));                                                             
    end
    handles.TARGET.STD_B0_SHIMMED = InterpolateToStandardFoVMex(temp_STD_B0_SHIMMED, coor_x, coor_y, coor_z, handles.settings.FOV, handles.settings.NPixel, 1);
    handles.TARGET.STD_B0_SHIMMED(handles.TARGET.STD_B0_SHIMMED==0)=nan;
                                 
    guidata(hObject, handles);
%     shimbox2soda(hObject);
    updateShimInfo(hObject);   
    UpdateShimmedField(hObject); 
    UpdateB0Map(hObject);
    UpdateMainImages(hObject);   
    setStatus(hObject, 0);
end


