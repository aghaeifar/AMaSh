function varargout = BasisMapsDisp(varargin)
% BASISMAPSDISP MATLAB code for BasisMapsDisp.fig
%      BASISMAPSDISP, by itself, creates a new BASISMAPSDISP or raises the existing
%      singleton*.
%
%      H = BASISMAPSDISP returns the handle to a new BASISMAPSDISP or the handle to
%      the existing singleton*.
%
%      BASISMAPSDISP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BASISMAPSDISP.M with the given input arguments.
%
%      BASISMAPSDISP('Property','Value',...) creates a new BASISMAPSDISP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BasisMapsDisp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BasisMapsDisp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BasisMapsDisp

% Last Modified by GUIDE v2.5 07-Dec-2016 17:35:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BasisMapsDisp_OpeningFcn, ...
                   'gui_OutputFcn',  @BasisMapsDisp_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BasisMapsDisp is made visible.
function BasisMapsDisp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BasisMapsDisp (see VARARGIN)

% Choose default command line output for BasisMapsDisp
handles.output = hObject;
handles.DATA = varargin{1};
sz = [size(handles.DATA), 1];
if numel(varargin) > 1
    handles.DISPLAY.Overlay = varargin{2};
else
    handles.DISPLAY.Overlay = ones(sz(1:3));
end


handles.DISPLAY.Alpha = zeros(size(handles.DISPLAY.Overlay));
if all(handles.DISPLAY.Overlay(:)) && ~sum(isnan(handles.DISPLAY.Overlay(:)))
    handles.DISPLAY.Alpha(handles.DISPLAY.Overlay>0) = 0;
else
    handles.DISPLAY.Alpha(handles.DISPLAY.Overlay>0) = 0.4;
end

% guidata(hObject, handles);
% handles = guidata(hObject);

handles.cSAGSlice = ceil((sz(1)+1)/2);
handles.cCORSlice = ceil((sz(2)+1)/2);
handles.cTRASlice = ceil((sz(3)+1)/2);

set(handles.sliderChannel,'Value',1)
set(handles.sliderChannel,'Min', 1)
set(handles.sliderChannel,'Max', size(handles.DATA,4))

try
    set(handles.sliderChannel,'SliderStep', [1/(sz(4)-1),1/(sz(4)-1)])
catch
    set(handles.sliderChannel,'SliderStep', [1,1])    
end
if sz(4) == 1
   set(handles.sliderChannel,'Enable', 'off'); 
end
set(handles.editChannel,'String',num2str(get(handles.sliderChannel,'Value')))

min_max = handles.DATA(~isnan(handles.DATA));
set(handles.editHzMax,'String', num2str(max(floor( mean(min_max(:)) + std(min_max(:)/2)), abs(ceil ( mean(min_max(:)) - std(min_max(:))/2)))));

axes(handles.axesSAG3)
imagesc(zeros(sz(1),sz(2)),[0 1]); axis off, axis xy;
axis off
colormap(gray)

axes(handles.axesCOR3)
imagesc(zeros(sz(2),sz(3)),[0 1]); axis off, axis xy;
axis off

axes(handles.axesTRA3)
imagesc(zeros(sz(1), sz(3)),[0 1]); axis off, axis xy;
axis off

guidata(hObject, handles);
UpdateShimBasisSets(hObject,str2num(get(handles.editChannel,'String')));
% handles = guidata(hObject);
UpdateMainImages(hObject);
        
% Update handles structure
% guidata(hObject, handles);

% UIWAIT makes BasisMapsDisp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BasisMapsDisp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editShimChannel_Callback(hObject, eventdata, handles)
% hObject    handle to editShimChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editShimChannel as text
%        str2double(get(hObject,'String')) returns contents of editShimChannel as a double



function editChannel_Callback(hObject, eventdata, handles)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);
    
    %% Check input
    if(str2num(get(handles.editChannel,'String')) > size(handles.DATA, 4))
        set(handles.editChannel,'String',num2str(size(handles.DATA, 4)))
    end

    if(str2num(get(handles.editChannel,'String')) < 1 )
        set(handles.editChannel,'String',num2str(1))
    end

    %% Update GUI
    set(handles.sliderChannel,'Value',str2num(get(handles.editChannel,'String')))
    
    %% Update basis shims
    guidata(hObject, handles);
    UpdateShimBasisSets(hObject,str2num(get(handles.editChannel,'String')));
    handles = guidata(hObject);

    %% Update main images
    guidata(hObject, handles);
    UpdateMainImages(hObject);
    handles = guidata(hObject);

    %% Update handles
    guidata(hObject, handles);


% --- Executes on slider movement.
function sliderChannel_Callback(hObject, eventdata, handles)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);

    %% Update GUI
    set(handles.sliderChannel,'Value', ceil(get(handles.sliderChannel,'Value')))
    set(handles.editChannel,'String',num2str(get(handles.sliderChannel,'Value')));

    %% Update basis sets
    guidata(hObject, handles);
    UpdateShimBasisSets(hObject,str2num(get(handles.editChannel,'String')));
    handles = guidata(hObject);

    %% Update main images
    guidata(hObject, handles);
    UpdateMainImages(hObject);
    handles = guidata(hObject);

    %% Update handles
    guidata(hObject, handles);


function editHzMax_Callback(hObject, eventdata, handles)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);   
    %% Update basis shims
    UpdateShimBasisSets(hObject,str2num(get(handles.editChannel,'String')));
    %% Update main images
    UpdateMainImages(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateMainImages(hObject)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);
    
    %% Update images for shim basis sets
    cla(handles.axesSAG3,'reset')
    imagesc(abs(permute(squeeze(handles.DISPLAY.ShimBasis(handles.cSAGSlice,:,end:-1:1,:)),[2 1 3])), 'Parent', handles.axesSAG3);
    hold(handles.axesSAG3,'on')
    I=imagesc(permute(squeeze(handles.DISPLAY.Overlay(handles.cSAGSlice,:,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesSAG3);
    set(handles.axesSAG3,'xtick',[],'ytick',[])
    set(handles.axesSAG3,'XColor',[0.941176470588235 0.941176470588235 0.941176470588235])
    set(handles.axesSAG3,'YColor',[0.941176470588235 0.941176470588235 0.941176470588235])
    set(I, 'AlphaData', permute(squeeze(handles.DISPLAY.Alpha(handles.cSAGSlice,:,end:-1:1,:)),[2 1 3]));

    cla(handles.axesCOR3,'reset')
    imagesc(permute(squeeze(handles.DISPLAY.ShimBasis(:,handles.cCORSlice,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesCOR3);
    hold(handles.axesCOR3,'on')
    I2=imagesc(permute(squeeze(handles.DISPLAY.Overlay(:,handles.cCORSlice,end:-1:1,:)),[2 1 3]), 'Parent', handles.axesCOR3);
    set(handles.axesCOR3,'xtick',[],'ytick',[])
    set(handles.axesCOR3,'XColor',[0.941176470588235 0.941176470588235 0.941176470588235])
    set(handles.axesCOR3,'YColor',[0.941176470588235 0.941176470588235 0.941176470588235])
    set(I2, 'AlphaData', permute(squeeze(handles.DISPLAY.Alpha(:,handles.cCORSlice,end:-1:1,:)),[2 1 3]));

    cla(handles.axesTRA3,'reset')
    imagesc(permute(squeeze(handles.DISPLAY.ShimBasis(:,:,handles.cTRASlice,:)),[2 1 3]), 'Parent', handles.axesTRA3);
    hold(handles.axesTRA3,'on')
    I3=imagesc(permute(squeeze(handles.DISPLAY.Overlay(:,:,handles.cTRASlice,:)),[2 1 3]), 'Parent', handles.axesTRA3);
    set(handles.axesTRA3,'xtick',[],'ytick',[])
    set(handles.axesTRA3,'XColor',[0.941176470588235 0.941176470588235 0.941176470588235])
    set(handles.axesTRA3,'YColor',[0.941176470588235 0.941176470588235 0.941176470588235])
    set(I3, 'AlphaData', permute(squeeze(handles.DISPLAY.Alpha(:,:,handles.cTRASlice,:)),[2 1 3]));
 
    %% Update handles structure
    guidata(hObject, handles);
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function UpdateShimBasisSets(hObject, Channel)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);

    cmap_size = 512;
    sz = size(handles.DATA);
    %% Convert to truecolor
    ConversionMap = bipolar(cmap_size);
    ConversionMap = [0 0 0;ConversionMap];

    %% Masked or not?
    SAGMag  = handles.DATA(handles.cSAGSlice,:,:,Channel);
    CORMag  = handles.DATA(:,handles.cCORSlice,:,Channel);
    TRAMag  = handles.DATA(:,:,handles.cTRASlice,Channel);

    %% Scale to range from 0..1
    Hzmax = str2num(get(handles.editHzMax,'String'));
    Hzmin = -Hzmax;
    SAGMag  = (SAGMag(:)-Hzmin)/(Hzmax-Hzmin);
    CORMag  = (CORMag(:)-Hzmin)/(Hzmax-Hzmin);
    TRAMag  = (TRAMag(:)-Hzmin)/(Hzmax-Hzmin);

    %% Convert
    handles.DISPLAY.ShimBasis   =   zeros([sz(1:3),3]); 
    Index = round(SAGMag *(length(ConversionMap)-1))+1;
    Index(Index==0) = 2;
    Index(isnan(Index)) = 1;
    Index(Index>cmap_size) = cmap_size;
    Index(Index<0) = 2;
    handles.DISPLAY.ShimBasis(handles.cSAGSlice,:,:,:) = reshape(ConversionMap(Index(:),:),[1,sz(2) ,sz(3),3]);

    Index = round(CORMag *(length(ConversionMap)-1))+1;
    Index(Index==0) = 2;
    Index(isnan(Index)) = 1;
    Index(Index>cmap_size) = cmap_size;
    Index(Index<0) = 2;
    handles.DISPLAY.ShimBasis(:,handles.cCORSlice,:,:) = reshape(ConversionMap(Index(:),:),[sz(1),1 ,sz(3),3]);

    Index = round(TRAMag *(length(ConversionMap)-1))+1;
    Index(Index==0) = 2;
    Index(isnan(Index)) = 1;
    Index(Index>cmap_size) = cmap_size;
    Index(Index<0) = 2;
    handles.DISPLAY.ShimBasis(:,:,handles.cTRASlice,:,:) = reshape(ConversionMap(Index(:),:),[sz(1),sz(2),1,3]);

    %% Update handles structure
    guidata(hObject, handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
    %% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);

    %% SAG   
    mousePoint = get(handles.axesSAG3, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    xrange = range(get(handles.axesSAG3, 'Xlim'));
    yrange = range(get(handles.axesSAG3, 'Ylim'));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        handles.cTRASlice = yrange-mouseY+1;
        handles.cCORSlice = mouseX;
        guidata(hObject, handles);
        UpdateShimBasisSets(hObject,str2num(get(handles.editChannel,'String')));
        UpdateMainImages(hObject);
        handles = guidata(hObject);
        guidata(hObject, handles);
    end
    
    %% COR    
    mousePoint = get(handles.axesCOR3, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    xrange = range(get(handles.axesCOR3, 'Xlim'));
    yrange = range(get(handles.axesCOR3, 'Ylim'));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        handles.cTRASlice = yrange-mouseY+1;
        handles.cSAGSlice = mouseX;
        guidata(hObject, handles);
        UpdateShimBasisSets(hObject,str2num(get(handles.editChannel,'String')));
        UpdateMainImages(hObject);
        handles = guidata(hObject);
        guidata(hObject, handles);
    end
    
    %% TRA  
    mousePoint = get(handles.axesTRA3, 'CurrentPoint');
    mouseX = round(mousePoint(1,1));
    mouseY = round(mousePoint(1,2));
    xrange = range(get(handles.axesTRA3, 'Xlim'));
    yrange = range(get(handles.axesTRA3, 'Ylim'));
    if(mouseX>0 && mouseX<xrange && mouseY>0 && mouseY<yrange)
        handles.cCORSlice = mouseY;
        handles.cSAGSlice = mouseX;
        guidata(hObject, handles);
        UpdateShimBasisSets(hObject,str2num(get(handles.editChannel,'String')));
        UpdateMainImages(hObject);
        handles = guidata(hObject);
        guidata(hObject, handles);
    end 
    
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
