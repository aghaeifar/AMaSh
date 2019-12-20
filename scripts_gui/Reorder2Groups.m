function Reorder2Groups(hObject)

handles = guidata(hObject); 
if handles.TARGET.SLCVOL.Dim == 3
   return; 
end


% we do this since we need just the real slice which will be measured, not the neighborhood
ind = round(size(handles.TARGET.MAP_B0,3)/2 + 0.5);

for i=1:handles.TARGET.SLCVOL.NGroups   
    MAP_B0_Group = handles.TARGET.MAP_B0(:,:,groups == i,:);
end

ind =  2*handles.flags.isDuplicateSession + 1;% 1 when is not duplicated session and 3 when it is
handles.TARGET.MAP_B0(:,:,thickRange,:).*handles.TARGET.Mask(:,:,thickRange,:)
handles.TARGET.MAP_B0(:,:,thickRange,:).*handles.TARGET.Mask(:,:,thickRange,:)


handles.TARGET.MAP_B0 = interpn(X,Y,Z, handles.B0_OBJ.STD_MAP_B0, X_dst, Y_dst, Z_dst);
handles.TARGET.Image  = interpn(X,Y,Z, handles.B0_OBJ.STD_Image, X_dst, Y_dst, Z_dst); % is required later?