% ------------------------------------------------------------------
% -------------------        Write coefficient to file    ----------

function writeShimCoef_Callback(hObject)

handles = guidata(hObject);
savepath = handles.settings.path2save_coef ;
fid = fopen([savepath, '/coef.txt'], 'w');
if fid == -1
    f_path = uigetdir;
    if ~f_path        
        return
    end
    handles.settings.path2save_coef  = f_path;
    fid = fopen(fullfile(handles.settings.path2save_coef , '/coef.txt'), 'w');
end

ShimCoef = handles.TARGET.NewShimSettings;

% if ~handles.flags.isShimBoxSODA2D || handles.TARGET.SLCVOL.Dim == 3 % if sequence is 3D or greenBox is selected
%     ShimCoef = handles.TARGET.NewShimSettings;
% else
%     sodaPos = handles.TARGET.SLCVOL.Position;
%     sodaSlc = handles.TARGET.SLCVOL.NSlices;
%     shimSlc = handles.TARGET.ShimBox.NSlices;
%     ShimBox = handles.TARGET.ShimBox;
%     ShimCoef = zeros(sodaSlc, handles.SHIM.NShimCha);
%     % The criteria is that center of slice stays in a volume, then shim coefficient of that volume will be used for this slice
%     for j=1:sodaSlc
%         for i=1:shimSlc
%             mnmx_x = [maxmin(ShimBox.Coords{i}.SAG); maxmin(ShimBox.Coords{i}.COR); maxmin(ShimBox.Coords{i}.TRA)];
%             if all(mnmx_x(:,1) <= sodaPos(j,:)') && all(sodaPos(j,:) <= mnmx_x(:,2)')
%                 ShimCoef(j,:) = handles.TARGET.NewShimSettings(i,:); 
%                 break;
%             end
%         end
%     end
% end

for i=1: size(ShimCoef, 1)
    freq = sprintf('%.0f ', ShimCoef(i,1));
    grad = sprintf('%.1f ', ShimCoef(i,2:4));
    mc   = sprintf('%.3f ', ShimCoef(i,5:end)); mc(end) = [];
    if numel(ShimCoef(i,5:end)) == 16
        mc = sprintf('%.3f ',[ShimCoef(i,5:12), zeros(1,16), ShimCoef(i,13:end)]);
%         mc = sprintf('%.3f ', [ShimCoef(i,5:end), zeros(1,16)]);
        mc(end) = [];
    end
    fprintf(fid, '%s%s%s', freq, grad, mc);
    fprintf(fid, '\n');
end    
% for i=1: size(handles.TARGET.NewShimSettings, 1)
%     freq = num2str(handles.TARGET.NewShimSettings(i,1), '%.0f');
%     grad = num2str(handles.TARGET.NewShimSettings(i,2:4), '%.1f ');
%     mc   = num2str(handles.TARGET.NewShimSettings(i,5:end), '%.3f ');
%     fprintf(fid, '%s %s %s', freq, grad, mc);
%     fprintf(fid, '\n');
% end
fclose(fid);
guidata(hObject, handles);
end % function

%%
function mnx = maxmin(data)
    mnx(1) = min(data(:));
    mnx(2) = max(data(:));
end