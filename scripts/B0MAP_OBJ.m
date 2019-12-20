classdef B0MAP_OBJ
%  B0MAP_OBJ:
%
%  A class to reconstruct B0 map 
%  External Dependenscies: uipickfiles 
%  Authors: Ali Aghaeifar (ali.aghaeifar@tuebingen.mpg.de), Christian Mirkes
%
    
    properties      
        filename = []
        mrprot
        SODA
        twix
        
        TE   
        mask_volume = 0; % mL
        
        RAW_complex= [];
        RAW_Image  = [];
        RAW_Mask   = [];
        MAP_B0_rad = [];
        MAP_B0  = [];
        MAP_B0_Wrapped = [];
        MAP_B0_Unwrapped = [];        
        
        STD_FOV      
        STD_NPixel   
        MaxSeparation ;
        
        STD_Image = []; 
        STD_Mask = [];
        STD_MAP_B0 = [];
        STD_MAP_B0_Wrapped = [];
        STD_MAP_B0_Unwrapped = [];
        
        unwrapping_mode
        reco_mode
        mask_mode
        mask_erode
        NShim_coils = 20;
        SHIM
    end
    
    methods
        function this = B0MAP_OBJ(varargin)            
            if(mod(nargin,2))
                error('Number of arguments must be even')
            end
            this = parse_inputs(this, varargin{:});                                    
            this = recoRAWDATA(this);  
            if strcmp(this.twix.image.softwareVersion, 'vb')
                this.SHIM = getShimSetting(this.twix, 36); % we have 36 (4+32) channels
            end

            if isempty(this.MaxSeparation)
                this.MaxSeparation = max([this.SODA.PixelSizeReadout,this.SODA.PixelSizePhase,this.SODA.PixelSizeSlice])/2;                
            end

            this = unwrapB0(this);   
            this = CalcMask(this, this.mask_erode);
            this.STD_Image = single(Interpolate2StdFov(this, this.RAW_Image));               
        end        
        
        %% --------------------------------------------------------------------------------------- %%
        function this = recoRAWDATA(this)
            if(isempty(this.filename))
                this.filename = uipickfiles('FilterSpec', pwd, 'Type', {'*.dat' 'Siemens rawdata file'});
                if isempty(this.filename)
                    error('Reconstruction was aborted by user');
                end
            end
            
            this.TE = [];
            raw_image = [];
            for cEco=1:numel(this.filename)
                [temp, this.mrprot, this.SODA, this.twix] = RecoImage(this.filename{cEco}); 
                checkConsistency(this); 
                this.TE = cat(2, this.TE, cell2mat(this.mrprot.MeasYaps.alTE(1:this.twix.image.NEco)) * 1e-6);
                raw_image = cat(8, raw_image, temp); % Concatenate along eco dimension
            end
            clear temp % save memory
            
            if numel(this.TE) == 1
                error('More than one echo is expected in rawfile')
            end            
            if strcmp(this.twix.hdr.MeasYaps.sKSpace.ucDimension, '0x2') || isequal(this.twix.hdr.MeasYaps.sKSpace.ucDimension, 2)% measurment is 2D, swap dimension of slices with partitions
                raw_image = permute(raw_image, [1 2 5 4 3 6:ndims(raw_image)]); % Col,Lin,Par,Cha,Slc/Slb,... -> Col,Lin,Slc/Slb,Cha,Par,...
            end
%             raw_image = applyHannFilt(raw_image);
            % Combine the coils
            disp('Combining the coils...');
            tempB0 = zeros(size(raw_image));
            if strcmpi(this.reco_mode, 'sos') % sum of squares
                for cEco=1:numel(this.TE) - 1
                    tempB0(:,:,:,:,:,:,:,cEco) = abs(raw_image(:,:,:,:,:,:,:,cEco+1)) .* raw_image(:,:,:,:,:,:,:,cEco+1) .* conj(raw_image(:,:,:,:,:,:,:,cEco)) .* abs(raw_image(:,:,:,:,:,:,:,cEco)); % weighting each echo by it's magnitude
                end
                this.RAW_Image = sum(raw_image(:,:,:,:,1,1,1,1).*conj(raw_image(:,:,:,:,1,1,1,1)), 4).^0.4; % sum over channels
            elseif strcmpi(this.reco_mode, 'adapt') % adaptive combine
                ceoff = [3 3 3];
                if size(raw_image, 3) == 1 % in the case of 2D single slice 
                    ceoff = [3 3];
                elseif size(raw_image, 3)<5 % Before was [3 3 3], for par<5 didn't work
                    ceoff = [1 1 1];
                end                
                [~,wfull,cmap] = openadapt_3D(permute(raw_image(:,:,:,:,1),[4 1 2 3]), ceoff, true);
                wfull = permute(wfull,[2 3 4 1]);
                for cEco = 1:numel(this.TE) - 1
                    tempB0(:,:,:,:,:,:,:,cEco) = bsxfun(@times, raw_image(:,:,:,:,:,:,:,cEco+1),wfull) .* conj(bsxfun(@times, raw_image(:,:,:,:,:,:,:,cEco),wfull));
                end 
                this.RAW_complex   = abs(sum(bsxfun(@times, raw_image(:,:,:,:,1,1,1,:),wfull),4) .* squeeze(sum(abs(cmap))).^2); % pick the first echo as magnitude image
                this.RAW_Image     = abs(sum(bsxfun(@times, raw_image(:,:,:,:,1,1,1,1),wfull),4) .* squeeze(sum(abs(cmap))).^2); % pick the first echo as magnitude image
            end            
            this.MAP_B0_rad = squeeze(angle(sum(tempB0, 4)));   
            this.MAP_B0_rad = this.MAP_B0_rad(:,:,:,1:end-1);
        end
        
        %% --------------------------------------------------------------------------------------- %%
        function DataSet_STD = Interpolate2StdFov(this, DataSet) % DataSet: [Col, Lin, Par/Slc/Slb]        
            % Compute maximal sepration of data points in the original image
            if this.SODA.Dim==2  % Stupid 2D data
                for i=1:this.SODA.NSlices
                    SAG(:,:,i) = this.SODA.Coords{i}.SAG(:,:,2); % 2 = center of the slice
                    COR(:,:,i) = this.SODA.Coords{i}.COR(:,:,2);
                    TRA(:,:,i) = this.SODA.Coords{i}.TRA(:,:,2);
                end
            else
                SAG = this.SODA.Coords{1}.SAG;
                COR = this.SODA.Coords{1}.COR;
                TRA = this.SODA.Coords{1}.TRA;
            end
%             DataSet_STD = raw2stdfov(DataSet, SAG, COR, TRA, this.STD_FOV, this.STD_NPixel);
            DataSet_STD = InterpolateToStandardFoV(DataSet, SAG, COR, TRA, this.STD_FOV, this.STD_NPixel, this.MaxSeparation);% InterpolateToStandardFoV needs 3D input
        end
        
        %% --------------------------------------------------------------------------------------- %%
        function this = unwrapB0(this, shift_2pi)
            if nargin < 2
                shift_2pi = 0;
            end
            % Calculate phase difference (B0 Map)    
            this.MAP_B0_Wrapped = this.MAP_B0_rad(:,:,:,1) / (this.TE(2) - this.TE(1)) / (2*pi);
            if strcmpi(this.unwrapping_mode, 'none')                
                this.MAP_B0_Unwrapped = this.MAP_B0_Wrapped;
            elseif strcmpi(this.unwrapping_mode, 'spatial')
                disp('Unwrapping B0Map (spatial)...');
                if isequal(this.twix.hdr.MeasYaps.sSliceArray.lSize, 1) && (strcmp(this.twix.hdr.MeasYaps.sKSpace.ucDimension, '0x2') || isequal(this.twix.hdr.MeasYaps.sKSpace.ucDimension,2)) % have to use isequal
                    this.MAP_B0_Unwrapped = this.MAP_B0_Wrapped;
                    warning('Unwrapping is not supported for single slice.');
                else
                    this.MAP_B0_Unwrapped = double(UnWrap_mex(single(this.MAP_B0_rad(:,:,:,1)+ shift_2pi*2*pi)) / (this.TE(2) - this.TE(1)) / (2*pi));
                end
            elseif strcmpi(this.unwrapping_mode, 'temporal')
                disp('Unwrapping B0Map (temporal)...');
                if numel(this.TE) < 3
                    error('Temporal unwrapping needs more than two echos');
                end
                this.MAP_B0_Unwrapped = unwrap_b0_temporal(this.MAP_B0_rad, this.TE); 
                % under construction
                % this.MAP_B0_Unwrapped = this.MAP_B0_rad(:,:,:,1,1,1,2) .* conj(this.MAP_B0_rad(:,:,:,1,1,1,1)) / (this.TE(3)+this.TE(1)-2*this.TE(2)) / (2*pi);
            end 
            this.MAP_B0 = this.MAP_B0_Wrapped;    
             % Interpolates
            this.STD_MAP_B0_Wrapped      = Interpolate2StdFov(this, this.MAP_B0_Wrapped);   this.STD_MAP_B0_Wrapped(this.STD_MAP_B0_Wrapped==0)=nan;
            this.STD_MAP_B0_Unwrapped    = Interpolate2StdFov(this, this.MAP_B0_Unwrapped); this.STD_MAP_B0_Unwrapped(this.STD_MAP_B0_Unwrapped==0)=nan;
            this.STD_MAP_B0              = this.STD_MAP_B0_Wrapped;
        end
        
        %% --------------------------------------------------------------------------------------- %%
        function this = CalcMask(this, ErodeKernelSize)
            this.mask_erode = ErodeKernelSize;
            if strcmpi(this.mask_mode, 'bet')
                if isunix
                    this.RAW_Mask = double(imerode(betFSL(abs(this.RAW_Image)), strel('sphere', ErodeKernelSize)));
                else
                    this.RAW_Mask = imerode(bet2(single(abs(this.RAW_Image))), strel('sphere', ErodeKernelSize));
                end
            elseif strcmpi(this.mask_mode, 'threshold')
                this.RAW_Mask=sqrt(sum(abs(this.RAW_Image).^2,4));
                this.RAW_Mask = this.RAW_Mask/max(this.RAW_Mask(:));
                T = 3*mean(col(this.RAW_Mask(1:4,1:4,1:4)));
                this.RAW_Mask(this.RAW_Mask<T) = 0;
                this.RAW_Mask(this.RAW_Mask>=T) = 1;
                this.RAW_Mask = imerode(this.RAW_Mask, ones(ErodeKernelSize,ErodeKernelSize,ErodeKernelSize));
                this.RAW_Mask = padarray(this.RAW_Mask,[1 1 1]);
                this.RAW_Mask = double(bwareaopen(this.RAW_Mask, double(round(sum(col(this.RAW_Mask))/2))));
                this.RAW_Mask = imdilate(this.RAW_Mask, ones(3,3,3));
                this.RAW_Mask = imfill(this.RAW_Mask, 'holes');
                this.RAW_Mask = imerode(this.RAW_Mask, ones(2,2,2));
                this.RAW_Mask = this.RAW_Mask(2:end-1, 2:end-1, 2:end-1);
            elseif strcmpi(this.mask_mode, 'none')
                this.RAW_Mask = ones(size(this.MAP_B0));
            end
            this.STD_Mask = Interpolate2StdFov(this, this.RAW_Mask);
            this.STD_Mask(this.STD_Mask<0.5)  = nan;
            this.STD_Mask(this.STD_Mask>=0.5) = 1;
            this.RAW_Mask(this.RAW_Mask<0.5)  = nan;
            this.RAW_Mask(this.RAW_Mask>=0.5) = 1;
            this.mask_volume = nansum(this.STD_Mask(:)) * (this.STD_FOV/this.STD_NPixel)^3 / 1000;
        end
        
        %% --------------------------------------------------------------------------------------- %%
        function checkConsistency(this)
            if this.twix.image.dataSize(8) < 2
                error('Less than 2 echos was found.');
            end
%             if this.mrprot.MeasYaps.sGroupArray.lSize ~= 1
%                 error('Multi group is not supported.');
%             end
            if this.twix.hdr.MeasYaps.sSliceArray.lSize ~= 1 && strcmp(this.twix.hdr.MeasYaps.sKSpace.ucDimension, '0x4')
                error('Multi slab in 3D acquisiton is not supported.');
            end
            
            % A warning message for user
            if (strcmp(this.twix.hdr.MeasYaps.sKSpace.ucDimension,'0x4') || sum(this.twix.hdr.MeasYaps.sKSpace.ucDimension == 4)) ...     % is 3D
               && strcmp(this.twix.hdr.MeasYaps.sTXSPEC.ucExcitMode, '0x2') ... % is non-selective excitation
               && this.twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness <100 % Volume thickness less than 100mm
               h = msgbox('Small volume is excited with a non-selective pulse! B0-map might not be correct.', 'Warning','Warn'); 
               set(h, 'WindowStyle', 'modal')
            end              
        end
        
        %% --------------------------------------------------------------------------------------- %%
        function this = parse_inputs(this, varargin)   
            parser = inputParser;
            parser.PartialMatching = false;

            addParameter(parser, 'file', [], @(x) assert(exist(x, 'file') ~= 0, 'Input file is not valid'));
            addParameter(parser, 'files', [], @(x) assert(iscell(x), 'Input files is not a cell'));
            addParameter(parser, 'fov', 300, @(x) assert(isnumeric(x) && (x>0), 'Invalid FOV'));
            addParameter(parser, 'nst', 201, @(x) assert(isnumeric(x) && (x>0), 'Invalid NST'));     
            addParameter(parser, 'maxsep', [], @(x) assert(isnumeric(x) && (x>0), 'Invalid Max Separation'));  
            addParameter(parser, 'unwrapping', 'spatial', @(x) assert(ismember(x, ["none","spatial","temporal"]), 'Invalid unwrapping method'));   
            addParameter(parser, 'coilcombine', 'adapt', @(x) assert(ismember(x, ["sos","adapt"]), 'Invalid coil combining method'));  
            addParameter(parser, 'mask', 'none', @(x) assert(ismember(x, ["none","bet","threshold"]), 'Invalid masking method'));     
            addParameter(parser, 'mask_erode', 1, @(x) assert(isnumeric(x) && (x>0), 'Invalid mask erode size'));  
            
            parse(parser, varargin{:});
            
            this.filename        = cell(0);
            if ~isempty(parser.Results.file)
                this.filename    = cellstr(parser.Results.file);
            elseif ~isempty(parser.Results.files)
                this.filename    = parser.Results.files;                
            end
            this.STD_FOV         = parser.Results.fov; 
            this.STD_NPixel      = parser.Results.nst;
            this.MaxSeparation   = parser.Results.maxsep;
            this.unwrapping_mode = parser.Results.unwrapping;
            this.reco_mode       = parser.Results.coilcombine;      
            this.mask_mode       = parser.Results.mask;  
            this.mask_erode      = parser.Results.mask_erode;  
        end        
   
    end % Methods  
end % classdef

