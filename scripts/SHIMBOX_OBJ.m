classdef SHIMBOX_OBJ
    %
    %  Read shim region
    %
    
    properties
        Dim
        MainOrientation
        
        FOV
        ReadoutFOV
        PhaseFOV
        Thickness
        
        GroupInd
        IsInterleaved = 1
        Position
        Normal
        
        NPixel
        NPixelReadout
        NPixelPhase
        NPixelSlice
        
        NSlices
        NGroups = 1
        Coords
        RotMatrix
        Sliceshift
        InplaneRot
        InPlaneRotMatrix
        
        LarmorConst = 42.5756 % [MHz/T]
        DefaultFileName = '\\10.41.60.40\upload\SODA_ADJ_DATA\SODA.ini';
    end
    
    methods
        function obj = SHIMBOX_OBJ(varargin)
            
            if nargin == 1
                addr = varargin{1};
                if exist(addr, 'file') ~= 2
                    error(['File path is not correct: ' addr])
                end
                obj = readBoxfromFile(obj,addr);
            elseif nargin == 2
                mrprot = varargin{1};
                boxName = varargin{2};
                obj = readBoxfromMrprot(obj, mrprot, boxName);
            end
            
            obj = calcPixLoc(obj);
            
            for i=1:obj.NSlices
                obj.Sliceshift{i} = dot(obj.Normal(i,:), obj.Position(i,:)); % [mm]
            end
        end
        
        %% % Get data from .ini file
        function obj = readBoxfromMrprot(obj, mrprot, boxName)
            if strcmp(boxName, 'SLCVOL')
                obj.NSlices = length(mrprot.MeasYaps.sSliceArray.asSlice);
                if(strcmp(mrprot.MeasYaps.sKSpace.ucDimension,'0x2') || sum(mrprot.MeasYaps.sKSpace.ucDimension == 2)) % 2D
                    obj.Dim = 2;
                elseif(strcmp(mrprot.MeasYaps.sKSpace.ucDimension,'0x4') || sum(mrprot.MeasYaps.sKSpace.ucDimension == 4)) % 3D
                    obj.Dim = 3;
                end
                
                for cSlc =1:obj.NSlices
                    % Slice dimensions
                    SliOS = 0;
                    if isfield(mrprot.MeasYaps.sKSpace, 'dSliceOversamplingForDialog')
                        if(mrprot.MeasYaps.sKSpace.dSliceOversamplingForDialog>0)
                            warning('Slice oversamplign detected: Check results carefully')
                            SliOS = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.dThickness*mrprot.MeasYaps.sKSpace.dSliceOversamplingForDialog;
                        end
                    end
                    obj.ReadoutFOV(cSlc) = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.dReadoutFOV;
                    obj.PhaseFOV(cSlc)   = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.dPhaseFOV;
                    obj.Thickness(cSlc)  = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.dThickness+SliOS;
                    obj.FOV(cSlc, :)     = [obj.ReadoutFOV(cSlc), obj.PhaseFOV(cSlc), obj.Thickness(cSlc)];
                    % Determine dimension & Pixels & Pixels size
                    obj.NPixelReadout(cSlc) = mrprot.MeasYaps.sKSpace.lBaseResolution;
                    obj.NPixelPhase(cSlc)   = mrprot.MeasYaps.sKSpace.lPhaseEncodingLines;
                    obj.NPixelSlice(cSlc)   = mrprot.MeasYaps.sKSpace.lPartitions;
                    obj.NPixel(cSlc,:)      = [obj.NPixelReadout(cSlc), obj.NPixelPhase(cSlc), obj.NPixelSlice(cSlc)]; 
                    if obj.Dim == 2
                        obj.NPixelSlice(cSlc)   = 3;
                    end
                    % In Plane Rotation
                    obj.InplaneRot(cSlc) = 0;
                    if isfield(mrprot.MeasYaps.sSliceArray.asSlice{cSlc}, 'dInPlaneRot')
                        obj.InplaneRot(cSlc) = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.dInPlaneRot;
                    end
                    % Normal
                    obj.Normal(cSlc,:) = [0,0,0];
                    if isfield(mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sNormal, 'dSag')
                        obj.Normal(cSlc,1) = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sNormal.dSag;
                    end
                    if isfield(mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sNormal, 'dCor')
                        obj.Normal(cSlc,2) = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sNormal.dCor;
                    end
                    if isfield(mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sNormal, 'dTra')
                        obj.Normal(cSlc,3) = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sNormal.dTra;
                    end
                    % Position
                    obj.Position(cSlc,:) = [0,0,0];
                    if isfield(mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sPosition, 'dSag')
                        obj.Position(cSlc,1) = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sPosition.dSag;
                    end
                    if isfield(mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sPosition, 'dCor')
                        obj.Position(cSlc,2) = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sPosition.dCor;
                    end
                    if isfield(mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sPosition, 'dTra')
                        obj.Position(cSlc,3) = mrprot.MeasYaps.sSliceArray.asSlice{cSlc}.sPosition.dTra;
                    end
                end
            elseif strcmp(boxName, 'ADJVOL')
				obj.Dim = 3;
				obj.NSlices = 1;
				%Slice dimensions
                obj.Thickness  = mrprot.MeasYaps.sAdjData.sAdjVolume.dThickness;
                obj.PhaseFOV   = mrprot.MeasYaps.sAdjData.sAdjVolume.dPhaseFOV;
                obj.ReadoutFOV = mrprot.MeasYaps.sAdjData.sAdjVolume.dReadoutFOV;    
                obj.FOV        = [obj.ReadoutFOV, obj.PhaseFOV, obj.Thickness];
                % Pixels
                obj.NPixelPhase   = obj.PhaseFOV;
                obj.NPixelReadout = obj.ReadoutFOV;
                obj.NPixelSlice   = obj.Thickness;
                obj.NPixel        = [obj.NPixelReadout, obj.NPixelPhase, obj.NPixelSlice]; 
                % Determine normal
				obj.Normal = [0,0,0];
                if isfield(mrprot.MeasYaps.sAdjData.sAdjVolume.sNormal, 'dSag')
                    obj.Normal(1) = mrprot.MeasYaps.sAdjData.sAdjVolume.sNormal.dSag;
                end
                if isfield(mrprot.MeasYaps.sAdjData.sAdjVolume.sNormal, 'dCor')
                    obj.Normal(2) = mrprot.MeasYaps.sAdjData.sAdjVolume.sNormal.dCor;
                end
                if isfield(mrprot.MeasYaps.sAdjData.sAdjVolume.sNormal, 'dTra')
                    obj.Normal(3) = mrprot.MeasYaps.sAdjData.sAdjVolume.sNormal.dTra;
                end
                % Determine centre of slice
                obj.Position = [0,0,0];
                if isfield(mrprot.MeasYaps.sAdjData.sAdjVolume.sPosition, 'dTra')
                    obj.Position(1) = mrprot.MeasYaps.sAdjData.sAdjVolume.sPosition.dSag;
                end
                if isfield(mrprot.MeasYaps.sAdjData.sAdjVolume.sPosition, 'dCor')
                    obj.Position(2) = mrprot.MeasYaps.sAdjData.sAdjVolume.sPosition.dCor;
                end
                if isfield(mrprot.MeasYaps.sAdjData.sAdjVolume.sPosition, 'dTra')
                    obj.Position(3) = mrprot.MeasYaps.sAdjData.sAdjVolume.sPosition.dTra;
                end
                % Inplane rotation
                obj.InplaneRot = 0;
                if isfield(mrprot.MeasYaps.sAdjData.sAdjVolume, 'dInPlaneRot')
                    obj.InplaneRot = mrprot.MeasYaps.sAdjData.sAdjVolume.dInPlaneRot;
                end                  
            elseif strcmp(boxName, 'NAVVOL')
                
            end
        end
        %% % Get data from .ini file
        function obj = readBoxfromFile(obj,filename)
            if(isempty(filename))
                warning(['Reading from ',obj.DefaultFileName])
                filename = obj.DefaultFileName;
            end
            
            % Dimension
            obj.Dim = cell2mat(inifile(filename,'read',{ 'Parameters','','Dim','i','none'}));
            % Slice
            obj.NSlices         = cell2mat(inifile(filename,'read',{ 'Parameters','','NSlices','i','none'}));
            obj.NGroups         = cell2mat(inifile(filename,'read',{ 'Parameters','','NGroups','i','none'}));
            obj.IsInterleaved   = cell2mat(inifile(filename,'read',{ 'Parameters','','Interleaved','i','none'}));
            
            % Matrix Dimension
            % obj.NPixelReadout	= cell2mat(inifile(filename,'read',{ 'MatrixSize','','Read','d','none'}));
            % obj.NPixelPhase     = cell2mat(inifile(filename,'read',{ 'MatrixSize','','Phase','d','none'}));
            % obj.NPixelSlice     = cell2mat(inifile(filename,'read',{ 'MatrixSize','','Slice','d','none'}));
            
            for cSlc = 1:obj.NSlices
                % FoV [FoVx FoVy Thickness]
                obj.FOV(cSlc,:) = cell2mat(inifile(filename,'read',{ 'FOV','',['FOV[',num2str(cSlc-1),']'],'d','none'}));
                obj.ReadoutFOV(cSlc)= obj.FOV(cSlc,1);
                obj.PhaseFOV(cSlc)  = obj.FOV(cSlc,2);
                obj.Thickness(cSlc) = obj.FOV(cSlc,3);
                % Matrix Dimension
                obj.NPixel(cSlc,:) = cell2mat(inifile(filename,'read',{ 'MatrixSize','',['MatSz[',num2str(cSlc-1),']'],'d','none'}));
                obj.NPixelReadout(cSlc)= obj.NPixel(cSlc,1);
                obj.NPixelPhase(cSlc)  = obj.NPixel(cSlc,2);
                obj.NPixelSlice(cSlc)  = obj.NPixel(cSlc,3);
                % In Plane Rotation
                obj.InplaneRot(cSlc) = cell2mat(inifile(filename,'read',{ 'InPlaneRot','',['PRot[',num2str(cSlc-1),']'],'d','none'}));
                % Normal
                obj.Normal(cSlc,:)   = cell2mat(inifile(filename,'read',{ 'Normal','',['NORM[',num2str(cSlc-1),']'],'d','none'}));
                % Position
                obj.Position(cSlc,:) = cell2mat(inifile(filename,'read',{ 'Position','',['POS[',num2str(cSlc-1),']'],'d','none'}));
                % Group
                obj.GroupInd{cSlc} = cell2mat(inifile(filename,'read',{ 'Group','',['GRP[',num2str(cSlc-1),']'],'i','none'})) + 1; % index starts from zero
            end
        end
        
        %% Calculate pixels location
        function obj = calcPixLoc(obj)
            % Loop over slices
            for cSlc =1:obj.NSlices
                % Determine dominant obj.MainOrientation
                [~,obj.MainOrientation(cSlc)] = max(abs(obj.Normal(cSlc)));
                % Calculate pixel positions for slice in logical coordinate system [PHASE, READ, SLICE]
                
                pFOV = obj.PhaseFOV(cSlc)   - obj.PhaseFOV(cSlc)  / obj.NPixelPhase(cSlc);
                rFOV = obj.ReadoutFOV(cSlc) - obj.ReadoutFOV(cSlc)/ obj.NPixelReadout(cSlc);
                sFOV = obj.Thickness(cSlc)  - obj.Thickness(cSlc) / obj.NPixelSlice(cSlc);
                if obj.Dim==2 % 2D
                    sFOV = obj.Thickness(cSlc);
                end
                
                [PHASE, READ, SLICE] = ndgrid(linspace( -pFOV/2, pFOV/2, obj.NPixelPhase(cSlc)  ),...
                    linspace(  rFOV/2,-rFOV/2, obj.NPixelReadout(cSlc)),...
                    linspace( -sFOV/2, sFOV/2, obj.NPixelSlice(cSlc)  ));
                
                % Transformation from logical to patient coordinate system
                if(obj.MainOrientation(cSlc)==1) % SAG
                    obj.InPlaneRotMatrix{cSlc} = [    0                    0                     1; ...
                        cos(obj.InplaneRot(cSlc))  sin(obj.InplaneRot(cSlc))   0;...
                        -sin(obj.InplaneRot(cSlc))  cos(obj.InplaneRot(cSlc))   0];
                    initNormal = [1; 0; 0];
                elseif(obj.MainOrientation(cSlc)==2) %COR
                    obj.InPlaneRotMatrix{cSlc} = [    cos(obj.InplaneRot(cSlc))  sin(obj.InplaneRot(cSlc))   0;...
                        0                    0                     1; ...
                        sin(obj.InplaneRot(cSlc))  -cos(obj.InplaneRot(cSlc))  0];
                    initNormal = [0; 1; 0];
                elseif(obj.MainOrientation(cSlc)==3) %TRA
                    obj.InPlaneRotMatrix{cSlc} = [    sin(obj.InplaneRot(cSlc))  -cos(obj.InplaneRot(cSlc))  0;...
                        cos(obj.InplaneRot(cSlc))  sin(obj.InplaneRot(cSlc))   0; ...
                        0                    0                     1];
                    initNormal = [0; 0; 1];
                end
                
                SAG = obj.InPlaneRotMatrix{cSlc}(1,1)*PHASE + obj.InPlaneRotMatrix{cSlc}(1,2)*READ + obj.InPlaneRotMatrix{cSlc}(1,3)*SLICE;
                COR = obj.InPlaneRotMatrix{cSlc}(2,1)*PHASE + obj.InPlaneRotMatrix{cSlc}(2,2)*READ + obj.InPlaneRotMatrix{cSlc}(2,3)*SLICE;
                TRA = obj.InPlaneRotMatrix{cSlc}(3,1)*PHASE + obj.InPlaneRotMatrix{cSlc}(3,2)*READ + obj.InPlaneRotMatrix{cSlc}(3,3)*SLICE;
                
                % Compute rotation matrix to align the logical coordinate normal vector with normal vector in the patient coordiate system [SAG, COR, TRA]
                if(norm(obj.Normal(cSlc,:))>1.001 || norm(obj.Normal(cSlc,:))<0.999 )
                    error('SODA.Normal must be unit vector')
                end
                
                v = cross(initNormal,obj.Normal(cSlc, :));
                s = norm(v);     % sine of angle
                c = dot(initNormal,obj.Normal(cSlc, :));   % cosine of angle
                
                if(s == 0)
                    obj.RotMatrix{cSlc} = eye(3)*c;
                else
                    V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                    obj.RotMatrix{cSlc} = eye(3)  + V + V*V*(1-c)/s^2;
                end
                
                % Apply rotation matrix, add shift and set READ as first dimension
                obj.Coords{cSlc}.SAG = permute(obj.RotMatrix{cSlc}(1,1)*SAG + obj.RotMatrix{cSlc}(1,2)*COR + obj.RotMatrix{cSlc}(1,3)*TRA + obj.Position(cSlc, 1),[2 1 3]);
                obj.Coords{cSlc}.COR = permute(obj.RotMatrix{cSlc}(2,1)*SAG + obj.RotMatrix{cSlc}(2,2)*COR + obj.RotMatrix{cSlc}(2,3)*TRA + obj.Position(cSlc, 2),[2 1 3]);
                obj.Coords{cSlc}.TRA = permute(obj.RotMatrix{cSlc}(3,1)*SAG + obj.RotMatrix{cSlc}(3,2)*COR + obj.RotMatrix{cSlc}(3,3)*TRA + obj.Position(cSlc, 3),[2 1 3]);
                
                obj.Coords{cSlc}.SAG = obj.Coords{cSlc}.SAG(end:-1:1,:,:);
                obj.Coords{cSlc}.COR = obj.Coords{cSlc}.COR(end:-1:1,:,:);
                obj.Coords{cSlc}.TRA = obj.Coords{cSlc}.TRA(end:-1:1,:,:);
            end  % for
        end  % function
    end % method
end % class

