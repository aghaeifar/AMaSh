classdef SODA_OBJ
    %SODA Slice orientation date
    %   All values are either in units of mm or radians if not specified
    %   otherwise
    % Example for usge: 
    % mrprot = readVB17Header();
    % SODA = SODA_OBJ('mrprot',mrprot);
    
    properties
            Dim
            MainOrientation
            
            ReadoutFOV
            PhaseFOV
            Thickness
            
            Position           
            Normal             
            
            NPixelReadout
            NPixelPhase
            NPixelSlice
            
            NPixelReadoutInterp = 0
            NPixelPhaseInterp   = 0
            NPixelSliceInterp   = 0
            
            PixelSizeReadout
            PixelSizePhase
            PixelSizeSlice   
            
            PixelSizeReadoutInterp = 0
            PixelSizePhaseInterp   = 0
            PixelSizeSliceInterp   = 0 
            
            NSlices
            Coords
            RotMatrix            
            Sliceshift
            InplaneRot
            InPlaneRotMatrix

            LarmorConst = 42.5756 % [MHz/T]
            DefaultFileName = '\\10.41.60.40\upload\SODA_ADJ_DATA\SODA.ini';
    end
    
    methods
        function obj = SODA_OBJ(varargin)
            
            if(mod(nargin,2))
                error('Number of arguments must be even: "File",File or "mrprot", mrprot')
            end
            
            filename = [];
            mrprot = [];
            ext = [];
       
            for i=1:2:nargin
                switch varargin{i}
                    case 'File'
                        [~,~,ext] = fileparts(varargin{i+1});
                        filename = varargin{i+1};
                    case 'mrprot'
                        mrprot = varargin{i+1};
                end
            end
                
            if(~isempty(mrprot))     
                    obj = extract(obj,mrprot); 
                    obj = calcPixLoc(obj);
                    for i=1:obj.NSlices
                        obj.Sliceshift{i} = dot( obj.Normal{i},obj.Position{i}); % [mm]
                    end   
            else
                 if(strcmp(ext,'.dat'))
                     mrprot = readVB17Header(filename);
                     obj = extract(obj,mrprot); 
                     obj = calcPixLoc(obj);
                     for i=1:obj.NSlices
                        obj.Sliceshift{i} = dot( obj.Normal{1},obj.Position{i}); % [mm]
                     end
                 elseif (strcmp(ext,'.dcm')||strcmp(ext,'.ima')||strcmp(ext,'.DCM')||strcmp(ext,'.IMA'))
                     obj = readSODAfromDICOM(obj,filename);
                     obj = calcPixLoc(obj);
                 else
                     error('Wrong file extension. Expecting "*.hdr" or "*.ini"')
                 end
            end   
        end
                             
        function obj = readSODAfromDICOM(obj, filename)
            info1   = dicominfo(filename);
            info2   = SiemensInfo(info1); % an external function
            obj.Dim = 2 + fix(info2.sKSpace.ucDimension/4); % 2 -> 2 or 4 -> 3
            % Slice
            obj.NSlices = info2.sSliceArray.lSize;
            
            % FoV
            obj.ReadoutFOV	= info2.sSliceArray.asSlice(1).dReadoutFOV;
            obj.PhaseFOV    = info2.sSliceArray.asSlice(1).dPhaseFOV;
            obj.Thickness   = info2.sSliceArray.asSlice(1).dThickness;
            
            % Matrix Dimension
            obj.NPixelReadout	= info2.sKSpace.lBaseResolution;
            obj.NPixelPhase     = info2.sKSpace.lPhaseEncodingLines;
            if(obj.Dim==2)
                obj.NPixelSlice   = 3;
            else
                obj.NPixelSlice   = info2.sKSpace.lPartitions;
            end            
            
            % Pixelsize            
            obj.PixelSizeReadout   = info2.sSliceArray.asSlice(1).dReadoutFOV/info2.sKSpace.lBaseResolution;
            obj.PixelSizePhase     = info2.sSliceArray.asSlice(1).dPhaseFOV/info2.sKSpace.lPhaseEncodingLines;
            if (obj.Dim==2)
                obj.PixelSizeSlice     = info2.sSliceArray.asSlice{cSlice}.dThickness;
            else
                obj.PixelSizeSlice     = info2.sSliceArray.asSlice(1).dThickness/info2.sKSpace.lPartitions;
            end
            
            obj.NSlices = length(info2.sSliceArray.asSlice);    
            for cSli = 1:obj.NSlices
                % In Plane Rotation
                try % this field doesn't exist in 3D sequence 
                    obj.InplaneRot{cSli} = info2.sSliceArray.asSlice(cSli).dInPlaneRot;
                catch
                    obj.InplaneRot{cSli} = 0;
                end
                % Normal
                try
                    NSag = info2.sSliceArray.asSlice(cSli).sNormal.dSag;
                catch
                    NSag = 0;
                end
                try
                    NTra = info2.sSliceArray.asSlice(cSli).sNormal.dTra;
                catch
                    NTra = 0;
                end
                try
                    NCor = info2.sSliceArray.asSlice(cSli).sNormal.dCor;
                catch
                    NCor = 0;
                end                
                obj.Normal{cSli} = [NSag;NCor;NTra];
                % Position
                obj.Position{cSli} =    cell2mat(struct2cell(info2.sSliceArray.asSlice(cSli).sPosition));             
            end
        end       
      
        function obj = extract(obj,mrprot)
            obj.NSlices = length(mrprot.MeasYaps.sSliceArray.asSlice);            
            for cSlice =1:obj.NSlices
                %% Determine normal
                try
                    NSag = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sNormal.dSag;
                catch
                    NSag = 0;
                end
                try
                    NTra = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sNormal.dTra;
                catch
                    NTra = 0;
                end
                try
                    NCor = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sNormal.dCor;
                catch
                    NCor = 0;
                end

                obj.Normal{cSlice} = [NSag; NCor; NTra];
 
                %% Determine centre of slice 
                try
                    if ischar(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dSag)
                        dSag = str2double(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dSag);
                    else
                        dSag = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dSag;
                    end
                catch
                    dSag = 0;
                end
                try
                    if ischar(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dCor)
                        dCor = str2double(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dCor);
                    else
                        dCor = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dCor;
                    end
                catch
                    dCor = 0;
                end
                try
                    if ischar(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dTra)
                        dTra = str2double(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dTra);
                    else
                        dTra = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dTra;
                    end
                catch
                    dTra = 0;
                end
                obj.Position{cSlice} = [dSag;dCor;dTra];

                %% Inplane rotation
                try
                    if ischar(mrprot.MeasYaps.sSliceArray.asSlice{1, 1}.dInPlaneRot)
                        dRot = str2double(mrprot.MeasYaps.sSliceArray.asSlice{1, 1}.dInPlaneRot);
                    else
                        dRot = mrprot.MeasYaps.sSliceArray.asSlice{1, 1}.dInPlaneRot;
                    end
                catch
                    dRot = 0;
                end
                obj.InplaneRot{cSlice} = dRot;

                %% Slice dimensions
                SliOS = 0;
                try
                    if(mrprot.MeasYaps.sKSpace.dSliceOversamplingForDialog>0)
                        warning('Slice oversamplign detected: Check results carefully')
                        SliOS = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dThickness*mrprot.MeasYaps.sKSpace.dSliceOversamplingForDialog;
                    end
                end
                obj.Thickness{cSlice}  = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dThickness+SliOS;
                obj.PhaseFOV{cSlice}   = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dPhaseFOV;
                obj.ReadoutFOV{cSlice} = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dReadoutFOV;

                %% Determine dimension & Pixels & Pixels size
                obj.NPixelPhase   = mrprot.MeasYaps.sKSpace.lPhaseEncodingLines;
                obj.NPixelReadout = mrprot.MeasYaps.sKSpace.lBaseResolution;
                
                obj.PixelSizePhase   = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dPhaseFOV/obj.NPixelPhase;
                obj.PixelSizeReadout = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dReadoutFOV/obj.NPixelReadout;
                
                if(strcmp(mrprot.MeasYaps.sKSpace.ucDimension,'0x2') || sum(mrprot.MeasYaps.sKSpace.ucDimension == 2)) % 2D
                    obj.Dim            = 2;
                    obj.NPixelSlice    = 3; % aa: we use 3 since we need to have volume 
                    obj.PixelSizeSlice = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dThickness;
                elseif(strcmp(mrprot.MeasYaps.sKSpace.ucDimension,'0x4') || sum(mrprot.MeasYaps.sKSpace.ucDimension == 4)) % 3D
                    obj.Dim            = 3;
                    obj.NPixelSlice    = mrprot.MeasYaps.sKSpace.lPartitions;
                    obj.PixelSizeSlice = (mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dThickness+SliOS)/obj.NPixelSlice;
                else
                    error('Unknown dimension');
                end                
            end            
        end
        
        
        function obj = calcPixLoc(obj)                           
            % Loop over slices
            for cSlice =1:obj.NSlices                
                % Determine dominant obj.MainOrientation
                [~,obj.MainOrientation{cSlice}] = max(abs(obj.Normal{cSlice}));                
                % Calculate pixel positions for slice in logical coordinate system [PHASE, READ, SLICE]
                pFOV = obj.PhaseFOV{cSlice} -obj.PhaseFOV{cSlice}/  obj.NPixelPhase;
                rFOV = obj.ReadoutFOV{cSlice} - obj.ReadoutFOV{cSlice}/obj.NPixelReadout;
                sFOV = obj.Thickness{cSlice} -obj.Thickness{cSlice}/ obj.NPixelSlice;
                if obj.Dim==2 % 2D
                    sFOV = obj.Thickness{cSlice};
                end
                
                [PHASE, READ, SLICE] = ndgrid(linspace( -pFOV/2, pFOV/2, obj.NPixelPhase  ),... 
                                              linspace(  rFOV/2,-rFOV/2, obj.NPixelReadout),...
                                              linspace( -sFOV/2, sFOV/2, obj.NPixelSlice  ));                                                     

                % Transformation from logical to patient coordinate system           
                if(obj.MainOrientation{cSlice}==1) % SAG
                    obj.InPlaneRotMatrix{cSlice} = [    0                    0                     1; ... 
                                                        cos(obj.InplaneRot{cSlice})  sin(obj.InplaneRot{cSlice})   0;...
                                                        -sin(obj.InplaneRot{cSlice})  cos(obj.InplaneRot{cSlice})   0];
                    initNormal = [1; 0; 0];  
                elseif(obj.MainOrientation{cSlice}==2) %COR
                    obj.InPlaneRotMatrix{cSlice} = [    cos(obj.InplaneRot{cSlice})  sin(obj.InplaneRot{cSlice})   0;...
                                                        0                    0                     1; ... 
                                                        sin(obj.InplaneRot{cSlice})  -cos(obj.InplaneRot{cSlice})  0];
                    initNormal = [0; 1; 0];  
                elseif(obj.MainOrientation{cSlice}==3) %TRA
                    obj.InPlaneRotMatrix{cSlice} = [    sin(obj.InplaneRot{cSlice})  -cos(obj.InplaneRot{cSlice})  0;...
                                                        cos(obj.InplaneRot{cSlice})  sin(obj.InplaneRot{cSlice})   0; ... 
                                                        0                    0                     1];
                    initNormal = [0; 0; 1];                       
                end
                
                SAG = obj.InPlaneRotMatrix{cSlice}(1,1)*PHASE + obj.InPlaneRotMatrix{cSlice}(1,2)*READ + obj.InPlaneRotMatrix{cSlice}(1,3)*SLICE;
                COR = obj.InPlaneRotMatrix{cSlice}(2,1)*PHASE + obj.InPlaneRotMatrix{cSlice}(2,2)*READ + obj.InPlaneRotMatrix{cSlice}(2,3)*SLICE;
                TRA = obj.InPlaneRotMatrix{cSlice}(3,1)*PHASE + obj.InPlaneRotMatrix{cSlice}(3,2)*READ + obj.InPlaneRotMatrix{cSlice}(3,3)*SLICE;
                
                % Compute rotation matrix to align the logical coordinate normal vector with normal vector in the patient coordiate system [SAG, COR, TRA]   
                if(norm(obj.Normal{cSlice})>1.001 || norm(obj.Normal{cSlice})<0.999 )
                   error('SODA.Normal must be unit vector') 
                end

                v = cross(initNormal,obj.Normal{cSlice});
                s = norm(v);     % sine of angle
                c = dot(initNormal,obj.Normal{cSlice});   % cosine of angle

                if(s == 0)
                   obj.RotMatrix{cSlice} = eye(3)*c;
                else
                   V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                   obj.RotMatrix{cSlice} = eye(3)  + V + V*V*(1-c)/s^2;
                end

                % Apply rotation matrix, add shift and set READ as first dimension
                obj.Coords{cSlice}.SAG = permute(obj.RotMatrix{cSlice}(1,1)*SAG + obj.RotMatrix{cSlice}(1,2)*COR + obj.RotMatrix{cSlice}(1,3)*TRA + obj.Position{cSlice}(1),[2 1 3]);
                obj.Coords{cSlice}.COR = permute(obj.RotMatrix{cSlice}(2,1)*SAG + obj.RotMatrix{cSlice}(2,2)*COR + obj.RotMatrix{cSlice}(2,3)*TRA + obj.Position{cSlice}(2),[2 1 3]);
                obj.Coords{cSlice}.TRA = permute(obj.RotMatrix{cSlice}(3,1)*SAG + obj.RotMatrix{cSlice}(3,2)*COR + obj.RotMatrix{cSlice}(3,3)*TRA + obj.Position{cSlice}(3),[2 1 3]);

                obj.Coords{cSlice}.SAG = obj.Coords{cSlice}.SAG(end:-1:1,:,:);
                obj.Coords{cSlice}.COR = obj.Coords{cSlice}.COR(end:-1:1,:,:);
                obj.Coords{cSlice}.TRA = obj.Coords{cSlice}.TRA(end:-1:1,:,:);
            end         
        end
    end
end

