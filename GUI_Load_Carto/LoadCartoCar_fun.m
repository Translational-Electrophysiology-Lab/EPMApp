function [CAR,Tabale_car,fig_h] = LoadCartoCar_fun(filename,do_plot);

if nargin<2
    do_plot=0;
end

tab_legend = {'Point','Backward Compatibility','Point Index','Backward Compatibility','X','Y','Z','Alpha(angular)','Beta(angular)','Gamma(angular)','Unipolar Value (mV)','Bipolar Value (mV)','LAT (ms)','Impedance (ohms)','PointType','LaberID','Point Label Size','Point Label','Label Color','Catherter ID','Is TE point'};
PointType = {'Normal = 0','Location Only= 1','Scar = 2','Floating = 3','TE= 4'};
LabelID = {'without Tag=–1','none = 4','His = 5','Pacing Site = 6','Double Potential = 7',
    'Fragmented Signal = 8', 'Ablation = 9', 'Scar = 10', 'Location Only = 11', 'TE (Transient Event for XP only) = 12'};


fid=fopen(filename);
a=fgetl(fid);
a=fgetl(fid);
if a==-1
   CAR=[];Tabale_car =[];
   return
end
d=textscan(a,'%s');
fseek(fid,0,-1);
a=fgetl(fid);
h = textscan(fid,[repmat('%s',[1 length(d{1})])],'collectoutput',1);

% h = textscan(fid,repmat('%s ',[1 length(tab_legend)+1]),'collectoutput',1);

Tabale_car = h{1};
xyz = str2double(h{1}(:,5:7));
Index_point = str2double(h{1}(:,3));
Indices.LAT = str2double(h{1}(:,13));
Indices.Unipolar = str2double(h{1}(:,11));
Indices.Bipolar = str2double(h{1}(:,12));
Indices.Catheter_ID = str2double(h{1}(:,20));

CAR.Index_point = Index_point;
CAR.Indices = Indices;
CAR.xyz = xyz;
CAR.tab_legend = tab_legend;
%%
xyz_mesh = xyz;
[A,iu] = unique(xyz,'stable','rows');
iun = setdiff([1:size(xyz,1)],iu);
iiko = sum(isnan(xyz),2)>0 | sum(isinf(xyz),2)>0;
iiko(iun) = 1;
xyz_mesh(iiko,:) = [];

DT = delaunayTriangulation(xyz_mesh);
CAR.MESH.vertices = xyz_mesh;

if ~isempty(DT.ConnectivityList)
    [triangles,v] = convexHull(DT); % mesh_original defines the faces used to do the mesh
    CAR.MESH.triangles= triangles;
else
    CAR.MESH.triangles= [];
end

Indices.LAT(iiko) = [];
Indices.Unipolar(iiko) = [];
Indices.Bipolar(iiko) = [];

if do_plot
    load('colormap_carto.mat')
    
    index_type_all = {'LAT','Unipolar','Bipolar'};
    
    
    %     fv.faces = MESH.triangles;
    fv.faces = DT.ConnectivityList;
    fv.vertices =  CAR.MESH.vertices;
    
    figure;
    for i = 1:length(index_type_all);
        index_type = index_type_all{i};
        ax(i) = subplot(1,3,i);
        fv.facevertexcdata = Indices.(index_type)(~iiko,:);
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

CAR.MESH.Indices = Indices;
fclose(fid);
