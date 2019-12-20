function newname = siemensRawRenamer(oldname)

if nargin < 1
    return;    
end
%%
if isfolder(oldname)
    list = dir(fullfile(oldname, '*.dat'));
    for i = 1:numel(list)
        [pathstr,name,ext] = fileparts(fullfile(oldname, list(i).name));
        k = strfind(name, '_');
        if numel(k) <3
           warning(['non-valid rawfile name: ' name]);
           continue;
        end
        name(k(end):end)=[];
        name(1:k(2)) = [];
        newname = fullfile(pathstr, [name,ext]);
        movefile(fullfile(oldname, list(i).name), newname);
    end
else
    [pathstr,name,ext] = fileparts(oldname);
    k = strfind(name, '_');
    if numel(k) <3
       error('non-valid rawfile name') ;
    end
    name(k(end):end)=[];
    name(1:k(2)) = [];
    newname = fullfile(pathstr, [name,ext]);
end
