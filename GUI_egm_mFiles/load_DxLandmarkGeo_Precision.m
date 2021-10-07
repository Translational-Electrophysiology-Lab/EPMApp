function [ParamOut] = load_DxLandmarkGeo_Precision(filename,do_plot);
% 
% filename = 'G:\Others\Velocity_Export\study_dwsG600573_2018_02_13_13_10_26\2018_03_16_17_57_16\DxLandmarkGeo.xml';
% [ParamOut] = load_DxLandmarkGeo_Precision(filename,0);


if nargin<2
    do_plot = 0;
end

ParamOut = [];
fid=fopen(filename);

str = '<Volumes number="';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a=='"');
ParamOut.Volumes_number = a(ii(1)+1:ii(2)-1);
a = fgetl(fid);
ii = find(a=='"');
ParamOut.Volumes_name = a(ii(1)+1:ii(2)-1);
a = fgetl(fid);
a = fgetl(fid);
str = 'Exported from study';
i1 = strfind(a,str)+length(str)+1;
i2 = strfind(a,' -->')-1;
ParamOut.Study_name =  a(i1:i2);

str = '<Vertices number="';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a=='"');
ParamOut.Vertices_number = str2double(a(ii(1)+1:ii(2)-1));
x = fscanf(fid,'%f');
ParamOut.MESH.vertices = reshape(x,[3 size(x,1)/3]).';
clear x ii
% =
str = 'Data values at each vertex of DxL map ';
while 1
    a = fgetl(fid);
        if contains(a,str);break;end
        if a==-1;return;end
end
n = a(findstr(a,str)+length(str) : findstr(a,str)+length(str)+2);
if isequal(n,'P-P')
    ParamOut.Data_type = 'Voltage';
elseif isequal(n,'LAT')
    ParamOut.Data_type = 'LAT';
else
    ParamOut.Data_type = 'NA';
end

%
str = '<Map_data number';
while 1
    a = fgetl(fid);
    if contains(a,str);break;end
    if a==-1;return;end
end
x = fscanf(fid,'%f');
ParamOut.Data = x;
% =
str = '<Map_color number';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
x = fscanf(fid,'%f');
ParamOut.Data_color = x;
a = fgetl(fid);a = fgetl(fid);a = fgetl(fid);
a = fgetl(fid);
ii = find(diff(isspace(a))==-1);
ParamOut.MESH.Data_color_lim(1) = str2double(a(ii(1)+1:ii(2)-1));
ParamOut.MESH.Data_color_lim(2) = str2double(a(ii(2)+1:end));
% =
str = '<Map_status number';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
x = fscanf(fid,'%f');
ParamOut.Data_status = x;
% =
str = '<Normals number';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
x = fscanf(fid,'%f');
ParamOut.MESH.Normals_number = reshape(x,[3 size(x,1)/3]).';
% =
str = '<Polygons number=';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
    
end
x = fscanf(fid,'%f');
ParamOut.MESH.Faces = reshape(x,[3 size(x,1)/3]).';
% =
str = '<Surface_of_origin number';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
x = fscanf(fid,'%f');
ParamOut.Surface_of_origin = x;
% =
str = '<Labels number';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a=='"');
lab_num = str2double(a(ii(1)+1:ii(2)-1));
ParamOut.Label_name = cell(1,lab_num);
ParamOut.Label_position = nan(lab_num,3);
str = 'Label name';
for ilab = 1:lab_num
    while 1
        a = fgetl(fid);
        if ~isempty(strfind(a,str));break;end
        if a==-1;return;end
    end
    ii = find(a=='"');
    ParamOut.Label_name{ilab} = a(ii(1)+1:ii(2)-1);
    a = fgetl(fid);
    ii = find(diff(isspace(a))==-1);
    ParamOut.Label_position(ilab,1) = str2double(a(ii(1)+1:ii(2)-1));
    ParamOut.Label_position(ilab,2) = str2double(a(ii(2)+1:ii(3)-1));
    ParamOut.Label_position(ilab,3) = str2double(a(ii(3)+1:end-1));
end
% =
str = '<Rotation';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii1 = find(a=='>');ii2 = find(a=='<');
a = a(ii1(1)+1:ii2(2)-1);
x = sscanf(a,'%f');
ParamOut.Rotation = reshape(x,[3 3]);
%
a = fgetl(fid);
ii1 = find(a=='>');ii2 = find(a=='<');
a = a(ii1(1)+1:ii2(2)-1);
ParamOut.Translation = sscanf(a,'%f');
%
a = fgetl(fid);
ii1 = find(a=='>');ii2 = find(a=='<');
a = a(ii1(1)+1:ii2(2)-1);
ParamOut.Scaling = sscanf(a,'%f');


fclose all;

if do_plot
    tri = ParamOut.MESH.Faces;
    xyz = ParamOut.MESH.vertices;
    c = ParamOut.Data;
    cnan=c;cnan(c==0)=nan;
    figure;
    ax(1) = subplot(121);
    h = trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3),1);
    
    ax(2) = subplot(122);
    h2 = trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3),cnan);
    set(h2,'edgeColor','none','faceColor','interp','facelighting','gouraud');
    
    linkprop(ax,'view');
    load('E:\UCL\Scripts_all\Scripts_mo\GUI_egm\GUI_egm_mFiles\Geo_Chann\colormap_carto.mat');
    colormap(cmap)
    
    set(h,'edgecolor','none','facecolor',[0 .4 .7],'facelighting','gouraud',...
        'ambientstrength',0.9,'specularcolorreflectance',0.4,'specularexponent',3,'specularstrength',0.3);
    
    l(1) = light('Position',[-1 -1 0],'Style','infinite','parent',ax(1));
    l(2) = light('Position',[-1 -1 0],'Style','infinite','parent',ax(2));
    
    axis(ax,'image','off')

    cl = [prctile(cnan,2) prctile(cnan,98)];
    if cl(1)==cl(2);
        cl = [min(cnan) max(cnan)];
    end
    ax(2).CLim = cl;
    cb = colorbar(ax(2));
    title(cb,ParamOut.Data_type);ylabel(cb,'(ms)');
    
    set(gcf,'units','centimeters','paperunits','centimeters','position',[3 3 20 8],'paperPosition',[0 0 20 8],'papersize',[20 8])
    
    set(cb,'position',[.90 .15 .025 .7])
    set(ax(1),'position',[0 0 .45 .95])
    set(ax(2),'position',[0.45 0 .45 .95])
    linkprop(ax,'view');
end


