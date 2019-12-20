% ------------------------------------------------------------------
% --------------      Interpolate to standard FOV     --------------
function myInterpolate(hObject, enShimMapIntplt, adjSlc)
    % enShimMapIntplt = apply interpolation for basis-maps
    % adjSlc = number of adjacent pixels to employ
    %    
    handles = guidata(hObject);     
    handles.TARGET.MAP_Shims = [];
    handles.TARGET.MAP_B0 = [];
    handles.TARGET.Mask = [];
    %% Select shim box
    %
    handles.flags.isShimBoxSODA2D = false;
    switch handles.settings.shim_scope.active
        case handles.settings.shim_scope.options{1} % Global shimming, use ADJVOL, green box
            ShimBox = handles.TARGET.ADJVOL;   
            
        case handles.settings.shim_scope.options{2} % Slice-Wise shimming, use SODA, yellow box            
            ShimBox = handles.TARGET.SLCVOL;
            if handles.TARGET.SLCVOL.Dim == 2
                ShimBox.Thickness(:) = (2*adjSlc+1)*ShimBox.Thickness;  % re-adjust slice thickness | 2*adjSlc -> adjacent pixels from both sides 
                ShimBox.NPixelSlice(:) = 4*adjSlc+ShimBox.NPixelSlice;  % 4*adjSlc+3 -> Each pixel includes 3 position values (start, middle, end) and we have 2*adjSlc+1 pixels
                ShimBox = ShimBox.calcPixLoc;                           % update SODA for slice over sampling
                handles.flags.isShimBoxSODA2D = true;
            end
            
        case handles.settings.shim_scope.options{3} % Group-Wise shimming, use NAVVOL, blue box
            ShimBox = handles.TARGET.NAVVOL; 
            if ShimBox.NSlices > handles.TARGET.SLCVOL.NSlices
                h = msgbox('Number of shimming regions can''t be more than number of slices/slabs', 'Warning','Warn'); set(h, 'WindowStyle', 'modal');
                return;
            end
    end
    %% Initiate variables
    %     
    handles.TARGET.MAP_B0 = cell(ShimBox.NSlices, 1);
    handles.TARGET.Image  = cell(ShimBox.NSlices, 1);
    handles.TARGET.Mask   = cell(ShimBox.NSlices, 1);
    if ~handles.flags.isDuplicateSession % otherwise STD&MEAN before shim will be set to zero
        handles.TARGET.MEAN = zeros(3, ShimBox.NSlices); % row1: before shim; row2: simulation; row3: after shim
        handles.TARGET.STD  = zeros(3, ShimBox.NSlices);  
        handles.TARGET.ShimBox = ShimBox;
    end
    
    %% organize slices/slabs to be shimmed
    %
%     handles.TARGET.indShimBox = 1:ShimBox.NSlices;
%     handles.TARGET.indSODA    = (1:handles.TARGET.SODA.NSlices) + ~handles.flags.isShimBoxSODA*ShimBox.NSlices;    
    [X,Y,Z] = ndgrid(linspace(-handles.settings.FOV/2, handles.settings.FOV/2, handles.settings.NPixel));  
    coor_x = []; coor_y = []; coor_z = [];
    coor_sz= zeros(ShimBox.NSlices, 1);
    % shim box coordinate [make a single vector for coordinate, otherwise "interpn" must be called per slab which make the whole calculation too slow. Slabs don't have same size necessarily]
    for cSlc=1:ShimBox.NSlices % # of slices or groups or adjustment volume
        coor_x = cat(1, coor_x, ShimBox.Coords{cSlc}.SAG(:));  
        coor_y = cat(1, coor_y, ShimBox.Coords{cSlc}.COR(:));
        coor_z = cat(1, coor_z, ShimBox.Coords{cSlc}.TRA(:));
        coor_sz(cSlc) = numel(ShimBox.Coords{cSlc}.SAG(:));     % save slabs sample size
    end
    coor_sz = [0; cumsum(coor_sz)]; % we need cumulative sum of the slabs size for the future indexing

    %% Interpolate data and calculate std
    %
    if handles.flags.isDataLoaded % Is B0 map imported?
        cprintf('text', 'Interpolating B0Map & Maks...\n') 
        % Interpolate B0 map to all of the slabs that are stored in the previous section
        MAP_B0 = interpn(X,Y,Z, handles.B0_OBJ.STD_MAP_B0, coor_x, coor_y, coor_z);
        Image  = interpn(X,Y,Z, handles.B0_OBJ.STD_Image , coor_x, coor_y, coor_z); % is required later? maybe
        if handles.flags.enMask
            Mask = interpn(X,Y,Z, handles.B0_OBJ.STD_Mask, coor_x, coor_y, coor_z);
            Mask(Mask<0.5) = nan;
            Mask(Mask>=0.5) = 1;
        else
            Mask = ones(size(coor_x));
        end        
                
        for cSlc=1:ShimBox.NSlices % # of slices or groups or adjustment volume
            ind = coor_sz(cSlc)+1:coor_sz(cSlc+1); % select indices for corresponding slice
            sz  = size(ShimBox.Coords{cSlc}.SAG);
            handles.TARGET.MAP_B0{cSlc} = reshape(MAP_B0(ind), sz);
            handles.TARGET.Image{cSlc}  = reshape(Image(ind), sz);
            handles.TARGET.Mask{cSlc}   = reshape(Mask(ind), sz);
            data_temp = MAP_B0(ind) .* Mask(ind); % to exclude voxels outside the mask
            if handles.flags.isShimBoxSODA2D  % Calculation of B0 deviation without slice oversampling when SODA is selected as shimBOX                
                data_temp  = reshape(data_temp, size(ShimBox.Coords{cSlc}.SAG));
                ind = ceil(size(data_temp, 3)/2); % we do this since we need just the real slice which will be measured, not the neighborhood
                data_temp  = data_temp(:,:,ind-1:ind+1);
            end
            % Attention: here STD doesn't include mask of basis-maps (but later for simulation result does)
            handles.TARGET.MEAN(2*handles.flags.isDuplicateSession + 1,cSlc) = mean(data_temp(:), 'omitnan'); % 1 when is not duplicated session and 3 when it is
            handles.TARGET.STD( 2*handles.flags.isDuplicateSession + 1,cSlc) = std( data_temp(:), 'omitnan'); % Calculate B0 deviation
        end           
        
        cprintf('green', 'Done\n')
    end
    
    if handles.flags.isShimMapsLoaded && enShimMapIntplt
        cprintf('text', 'Interpolating Shim basis-sets...\n');
        MAP_Shims = zeros([numel(coor_x), handles.SHIM.NShimCha]);
        for cCha=1:size(handles.SHIM.STD_B0, 4)
            MAP_Shims(:,cCha) = interpn(X,Y,Z, handles.SHIM.STD_B0(:,:,:,cCha) .* handles.SHIM.STD_Mask, coor_x, coor_y, coor_z);
        end
        
        handles.TARGET.MAP_Shims = cell(ShimBox.NSlices, 1);
        for cSlc=1:ShimBox.NSlices % # of slices or groups or adjustment volume
            ind = coor_sz(cSlc)+1:coor_sz(cSlc+1);
            handles.TARGET.MAP_Shims{cSlc} = MAP_Shims(ind,:);
        end
        cprintf('green', 'Done\n')
    end 
       
    guidata(hObject, handles)