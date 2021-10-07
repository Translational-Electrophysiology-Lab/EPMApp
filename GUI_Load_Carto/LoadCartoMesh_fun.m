function [CartoData,fig_h] = LoadCartoMesh_fun(filename,do_plot);

if nargin<2
    do_plot=0;
end
% clear
% %% Type of mesh
% name = '1-Epi.mesh';
% % name = '2-Endo.mesh';
%
% %% Select Index Type
% % Indices_name = {'NodeNumber','Unipolar','Bipolar','LAT','Impedance','A1','A2','A2-A1','SCI','ICL','ACL','Force','Paso'};
% index_type = 'Bipolar';
% %%
%
% %%
% do_save = 1;
%
% %%
% path_name = 'E:\UCL\Data&People\Data_Barts_VT\151211_pt2\Patient 2015_12_11\Study 1\Export_Study-1-12_11_2015-15-51-21\';
% dir_save = 'E:\UCL\Data&People\Data_Barts_VT\151211_pt2\MAT\';
% filename =   [path_name,name];
fid=fopen(filename);
%% Header
seekstr = 'GeneralAttributes';
a=fgetl(fid);
while 1
    a=fgetl(fid);
    if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
    if (strfind(a, seekstr)), break, end
end
celltraces = textscan(fid,'%s%s',100,'delimiter','=');


headname = {'MeshID','MeshName','NumVertex','NumTriangle','TopologyStatus','MeshColor','Matrix','NumVertexColors','ColorsIDs','ColorsNames'};
for i = 1:length(headname)
    ii = find(~cellfun(@isempty,strfind(celltraces{1},headname{i})));
    if ~isempty(ii)
    Header.(headname{i}) = celltraces{2}{ii};
    end
end

%% Verices
fseek(fid,0,-1);
seekstr = 'VerticesSection';
while 1
    a=fgetl(fid);
    if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
    if (strfind(a, seekstr)), break, end
end
h = textscan(fid,'%s %s %s %s %s %s %s %s',1,'collectoutput',1);
h{1}{1} = '# Node';
Vertices_label = h{1};
a=fgetl(fid);
celltraces = textscan(fid,'%f %f %f %f %f %f %f %f', 'delimiter', '=','collectoutput',1);
Vertices = cell2mat(celltraces);

%% TrianglesSection
seekstr = 'TrianglesSection';
while 1
    a=fgetl(fid);
    if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
    if (strfind(a, seekstr)), break, end
end
h = textscan(fid,'%s %s %s %s %s %s %s %s',1,'collectoutput',1);
h{1}{1} = '# Node';
Triangles_label = h{1};
a=fgetl(fid);
celltraces = textscan(fid,'%f %f %f %f %f %f %f %f', 'delimiter', '=','collectoutput',1);
Triangles = cell2mat(celltraces);

%% [VerticesColorsSection]
fseek(fid,0,-1);
a=fgetl(fid);
seekstr = 'VerticesColorsSection';
while ~isequal(a,-1)
    a=fgetl(fid);
    %     if ~ischar(a), display('No VerticesColorsSection'); , end
    if (strfind(a, seekstr)), break, end
end
if ~isequal(a,-1)
    a=fgetl(fid);
    a=fgetl(fid);
    a=fgetl(fid);
    col = Header.ColorsNames;
    ii = find(col==' ');ii2 = find(diff(ii)>1);
    ii3 = [1 ii(ii2) ii(end)];
    for j = 1:length(ii3)-1
        hh = col(ii3(j):ii3(j+1));
        hh(isspace(hh))=[];
        Indices_name{j} = (hh);
    end
    Indices_name= ['NodeNumber',Indices_name];
    h = textscan(fid,repmat('%f ',[1,length(Indices_name)]),'delimiter', '=','collectoutput',1);
    hmat = cell2mat(h);
    hmat(hmat==-10000)=nan;
    
    for i=1:length(Indices_name);
        fname = Indices_name{i};
        fname(fname=='-') = '_';
        try
        Indices.(fname) = hmat(:,i);
        end
    end
else
    Indices = [];
    disp('No VerticesColorsSection')
end

%% [VerticesAttributesSection]
fseek(fid,0,-1);
a=fgetl(fid);
seekstr = 'VerticesAttributesSection';
while ~isequal(a,-1)
    a=fgetl(fid);
%     if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
    if (strfind(a, seekstr)), break, end
end
h = textscan(fid,'%s',5,'delimiter',',','collectoutput',1);
VerticesAttributes.Note = h{1};
a=fgetl(fid);
clear h
h = textscan(fid,'%f %f %f','delimiter', '=','collectoutput',1);
VerticesAttributes.NodeNumber = h{1}(:,1);
VerticesAttributes.Scar = h{1}(:,2);
VerticesAttributes.EML = h{1}(:,3);


%% Mesh and plot
MESH.triangles = Triangles(:,2:4)+1;
MESH.vertices = Vertices(:,2:4);
MESH.vertices_norm = Vertices(:,5:7);


if do_plot & ~isempty(Indices)
    load('colormap_carto.mat')

    index_type_all = {'LAT','Unipolar','Bipolar'};
    
    fv.faces = MESH.triangles;
    fv.vertices =  MESH.vertices;
    
    figure;
    for i = 1:length(index_type_all);
        index_type = index_type_all{i};
        ax(i) = subplot(1,3,i);
        fv.facevertexcdata = Indices.(index_type);
        if mean(isnan(fv.facevertexcdata))==1
            error([index_type,': not available'])
        end
        
        
        hpatch = patch(fv,'Parent',ax(i));
        set(hpatch,'edgecolor','none','facecolor','interp','facelighting','phong',...,
            'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
            'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
        cb(i) = colorbar;
        caxis([prctile(fv.facevertexcdata,[2 98])])
        title(index_type)
    end
    set(gcf,'units','centimeters','position',[2 2 18 8],'paperposition',[0 0 18 8],'papersize',[18 8])
    set(ax(1),'position',[0 0 .3 .9])
    set(ax(2),'position',[0.32 0 .3 .9])
    set(ax(3),'position',[0.64 0 .3 .9])
    
    set(cb(1),'position',[0.28 .1 .02 .8])
    set(cb(2),'position',[0.28+.32 .1 .02 .8])
    set(cb(3),'position',[0.28+.64 .1 .02 .8])
    
    xlabel(cb(1),'(ms)')
    xlabel(cb(2),'(mV)')
    xlabel(cb(3),'(mV)')
    
    axis(ax,'off','equal')
    hlink=linkprop(ax,'view');
    key = 'graphics_linkprop';
    % Store link object on first subplot axes
    setappdata(ax(1),key,hlink); 
            
    colormap(cmap(1:round(size(cmap,1)/8-1):end,:));
    fig_h.ax = ax;
    fig_h.cb = cb;
else
    fig_h.ax = [];
    fig_h.cb = [];
end

CartoData.Indices = Indices;
CartoData.MESH = MESH;
CartoData.Header = Header;
name = filename;
ii = find(filename=='\');
CartoData.name = name(ii(end)+1:end-5);

%
% %% save mat
% filename = [dir_save,name(1:end-5)];
% display(['Saving: ',filename]);
