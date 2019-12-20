function recoB0Map(hObject, dataType)
% Get correct handles (handles provided by the function arguments are not working)
    handles = guidata(hObject);

    if nargin <2
       dataType = 'rawdata'; 
    end
    % Get raw data file
    info = 'Please select binary file to read';
    [filename, pathname]=uigetfile(fullfile(handles.settings.path2data,'*.dat'),info);
   
    if(pathname == 0)
        return;
    end
    handles.settings.path2data = pathname;
    guidata(hObject, handles)
    setStatus(hObject, 1);
    set(handles.edit_shimInfo, 'String', 'Please wait...');
    drawnow;    
    
    %% Reconstruct B0 map
    reco_mode = {'sos', 'adapt'};  % the orders must be identical with setting GUI
    unwrapping_method = {'none', 'spatial', 'temporal'};
    handles.B0_OBJ = B0MAP_OBJ('File',[pathname,filename],'NST',handles.settings.NPixel,'FOV',handles.settings.FOV, 'coilcombine', reco_mode(handles.settings.reco_mode), 'unwrapping', unwrapping_method(handles.settings.unwrapping_method));
    handles.flags.isDataLoaded = true;
    fclose all;          
        
    %% Convert Background to truecolor (grayscale)
    handles.B0_OBJ.STD_Image = abs(handles.B0_OBJ.STD_Image(:,:,:,end))/max(abs(handles.B0_OBJ.STD_Image(:)));
    if isequal(size(handles.B0_OBJ.STD_Mask), size(handles.DISPLAY.Mask)) % use old mask
        handles.B0_OBJ.STD_Mask = handles.DISPLAY.Mask; 
    else
        handles.DISPLAY.Mask = ones(size(handles.B0_OBJ.STD_Image));
    end
    
    ConversionMap = [0,0,0;8,8,8;15,15,15;22,22,22;30,30,30;37,37,37;44,44,44;52,52,52;59,59,59;66,66,66;74,74,74;81,81,81;88,88,88;...
                     96,96,96;103,103,103;110,110,110;118,118,118;125,125,125;132,132,132;140,140,140;147,147,147;154,154,154;162,162,162;...
                     169,169,169;176,176,176;184,184,184;191,191,191;198,198,198;206,206,206;213,213,213;220,220,220;228,228,228;235,235,235;...
                     242,242,242;250,250,250;257,257,257;264,264,264;271,271,271;279,279,279;286,286,286;293,293,293;301,301,301;308,308,308;...
                     315,315,315;323,323,323;330,330,330;337,337,337;345,345,345;352,352,352;359,359,359;367,367,367;374,374,374;381,381,381;...
                     389,389,389;396,396,396;403,403,403;411,411,411;418,418,418;425,425,425;433,433,433;440,440,440;447,447,447;455,455,455;...
                     462,462,462;469,469,469;477,477,477;484,484,484;495,495,495;497,497,497;500,500,500;503,503,503;505,505,505;508,508,508;...
                     511,511,511;513,513,513;516,516,516;519,519,519;522,522,522;524,524,524;527,527,527;530,530,530;532,532,532;535,535,535;...
                     538,538,538;540,540,540;543,543,543;546,546,546;548,548,548;551,551,551;554,554,554;557,557,557;559,559,559;562,562,562;...
                     565,565,565;567,567,567;570,570,570;573,573,573;575,575,575;578,578,578;581,581,581;583,583,583;586,586,586;589,589,589;...
                     591,591,591;594,594,594;597,597,597;600,600,600;602,602,602;605,605,605;608,608,608;610,610,610;613,613,613;616,616,616;...
                     618,618,618;621,621,621;624,624,624;626,626,626;629,629,629;632,632,632;635,635,635;637,637,637;640,640,640;643,643,643;...
                     645,645,645;648,648,648;651,651,651;653,653,653;656,656,656;659,659,659;661,661,661;664,664,664;667,667,667;670,670,670;...
                     672,672,672;675,675,675;678,678,678;680,680,680;683,683,683;686,686,686;688,688,688;691,691,691;694,694,694;696,696,696;...
                     699,699,699;702,702,702;705,705,705;707,707,707;710,710,710;713,713,713;715,715,715;718,718,718;721,721,721;723,723,723;...
                     726,726,726;729,729,729;731,731,731;734,734,734;737,737,737;739,739,739;742,742,742;745,745,745;748,748,748;750,750,750;...
                     753,753,753;756,756,756;758,758,758;761,761,761;764,764,764;766,766,766;769,769,769;772,772,772;774,774,774;777,777,777;...
                     780,780,780;783,783,783;785,785,785;788,788,788;791,791,791;793,793,793;796,796,796;799,799,799;801,801,801;804,804,804;...
                     807,807,807;809,809,809;812,812,812;815,815,815;818,818,818;820,820,820;823,823,823;826,826,826;828,828,828;831,831,831;...
                     834,834,834;836,836,836;839,839,839;842,842,842;844,844,844;847,847,847;850,850,850;853,853,853;855,855,855;858,858,858;...
                     861,861,861;863,863,863;866,866,866;869,869,869;871,871,871;874,874,874;877,877,877;879,879,879;882,882,882;885,885,885;...
                     887,887,887;890,890,890;893,893,893;896,896,896;898,898,898;901,901,901;904,904,904;906,906,906;909,909,909;912,912,912;...
                     914,914,914;917,917,917;920,920,920;922,922,922;925,925,925;928,928,928;931,931,931;933,933,933;936,936,936;939,939,939;...
                     941,941,941;944,944,944;947,947,947;949,949,949;952,952,952;955,955,955;957,957,957;960,960,960;963,963,963;966,966,966;...
                     968,968,968;971,971,971;974,974,974;976,976,976;979,979,979;982,982,982;984,984,984;987,987,987;990,990,990;992,992,992;...
                     995,995,995;998,998,998;1000,1000,1000]/1000;
                 
    Index = round(handles.B0_OBJ.STD_Image*length(ConversionMap));
    Index(Index==0) = 1;
    Index(isnan(Index)) = 1;
    
    handles.DISPLAY.Image = ConversionMap(Index(:),:);
    handles.DISPLAY.Image = reshape(handles.DISPLAY.Image, handles.settings.NPixel ,handles.settings.NPixel ,handles.settings.NPixel, 3);
    handles.DISPLAY.Min = 0;
    handles.DISPLAY.Max = max(handles.DISPLAY.Image(:));              
    % Set filename
    [~, Filename] = fileparts(handles.B0_OBJ.filename{1});
    set(handles.mainFigure,'Name', ['B0GUI - ', Filename]);
    guidata(hObject, handles);
    
    UpdateShimmedField(hObject); % update will be applied inside next line: checkbox_unwrap_Callback();
    
    % check unwrap status
    if handles.flags.enUnwrap
        handles.B0_OBJ.MAP_B0 = handles.B0_OBJ.MAP_B0_Unwrapped;
        handles.B0_OBJ.STD_MAP_B0 = handles.B0_OBJ.STD_MAP_B0_Unwrapped;
    else 
        handles.B0_OBJ.MAP_B0 = handles.B0_OBJ.MAP_B0_Wrapped;
        handles.B0_OBJ.STD_MAP_B0 = handles.B0_OBJ.STD_MAP_B0_Wrapped;
    end
    guidata(hObject,handles );
    
    % Calc Init mean
    if handles.flags.SODALoaded
        myInterpolate(hObject, false, handles.settings.adjacent_slices);
    end      
    
    handles = guidata(hObject); % Since hObject has been modified in [checkbox_unwrap_Callback, UpdateShimmedField, myInterpolate]
    if isfield(handles, 'TARGET') && isfield(handles.TARGET, 'oldSHIM') && ~isequal(handles.TARGET.oldSHIM.SH, handles.B0_OBJ.SHIM.SH)
        h = msgbox('Scanner''s shim setting is changed!', 'Warning','Warn'); 
        set(h, 'WindowStyle', 'modal')
    end
    handles.TARGET.oldSHIM = handles.B0_OBJ.SHIM;
    guidata(hObject, handles);
    
    UpdateShimmedField(hObject);
    UpdateB0Map(hObject);
    UpdateMainImages(hObject);
    updateShimInfo(hObject);
    setStatus(hObject, 0);    