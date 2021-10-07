%% combination of surf_index_interp.m and surf_index_v3

function [H,index_val_interp,geo] = patch_index_mo(var1,var2,var3,var4,var5,var6);
%% plot surface
%% If first argument is NOT an axis handle:
% surf_index_mo([1:240],xx(:,k)-min(x(:,k)),'mo_sock2',1,0);
% var1 : number of channels to plot
% var2 : values of channels
% var3 : sock name
% var4 : interpolation; 0=no;1=yes;2=no but better representation of sock
% var5 : allow negative values; 0=no [default]; 1=yes
%% If first argument is  an axis handle:
% surf_index_mo(gca,[1:240],xx(:,k)-min(x(:,k)),'mo_sock2',1,1);
% var1 : axis handle
% var2 : number of channels to plot
% var3 : values of channels
% var4 : sock name
% var5 : interpolation; 0=no;1=yes;2=no but better representation of sock
% var6 : allow negative values; 0=no [default]; 1=yes
%%
if ishandle(var1)
    ax = var1;
    index_chan = var2;
    index_val = var3;
    sock_name = var4;
    if nargin>4
        flag_interp = var5;
    else
        flag_interp = 0;
    end
    if nargin==6
        allow_neg_val = var6;
    else
        allow_neg_val = 0;
    end
else
    ax = gca;
    index_chan = var1;
    index_val = var2;
    sock_name = var3;
    if nargin>3
        flag_interp = var4;
    else
        flag_interp = 0;
    end
    if nargin==5
        allow_neg_val = var5;
    else
        allow_neg_val = 0;
    end
end
clear var*
%%
addpath .\GUI_egm_mFiles\Geo_Chann
load(['ALLgeoDATA_',sock_name,'_mesh']);
load(['ALLgeoDATA_',sock_name]);
%%
% Input is already interp
if length(index_val)==size(MESH.Vertices_Interp,1)
    flag_interp = 3;
end

if flag_interp==1
    %%
    v = nan(size(MESH.Vertices_Original,1),1);
    for i=1:length(index_val)
        a = find(index_chan(i) == channel_num);
        if ~isnan(channel_num(index_chan(i)))
            v(a) = index_val(i);
        end
    end
    iiok = ~isnan(v);
    F = scatteredInterpolant(MESH.Vertices_Original(iiok,1),MESH.Vertices_Original(iiok,2),MESH.Vertices_Original(iiok,3),v(iiok));
    %
    fv.faces = MESH.Faces_OnlyVentricles_Interp;
    fv.vertices = MESH.Vertices_Interp;
    index_val_interp = F([MESH.Vertices_Interp(:,1),MESH.Vertices_Interp(:,2),MESH.Vertices_Interp(:,3)]);
    fv.facevertexcdata = index_val_interp;
    %
    fvtop.faces = MESH.Faces_Top_Interp;
    fvtop.vertices = MESH.Vertices_Interp;
    fvtop.facevertexcdata = 0;
    
    hpatch = patch(fv,'Parent', ax);
    hold(ax,'on')
    hpatch_top = patch(fvtop,'Parent', ax);
    set(hpatch,'edgecolor','none','facecolor','interp','facelighting','phong',...,
        'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
        'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
    set(hpatch_top,'edgecolor','none','facecolor',[1 1 1]*.35,'facelighting','phong',...,
        'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
        'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
    
    set(ax,'clim',[prctile(index_val_interp,5) prctile(index_val_interp,95)])
    
    
elseif flag_interp==0
    index_val_interp = [];
    v = nan(1,length(channel_num));
    v_ind = nan(1,length(channel_num));
    for i= 1:length(index_chan)
        a = find(index_chan(i) == channel_num);
        if isnan(v(a))
            v(a) =index_val(i);
            v_ind(a) =i;
        end
    end
    v(isnan(v)) = -10;
    %
    fv.faces = MESH.Faces_OnlyVentricles;
    fv.vertices = MESH.Vertices_Original;
    fv.facevertexcdata = v(:);
    %
    fvtop.faces = MESH.Faces_Top;
    fvtop.vertices = MESH.Vertices_Original;
    fvtop.facevertexcdata = 0;
    
    hpatch = patch(fv,'Parent', ax);
    hold(ax,'on')
    hpatch_top = patch(fvtop,'Parent', ax);
    set(hpatch,'edgecolor','none','facecolor','flat','facelighting','phong',...,
        'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
        'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
    set(hpatch_top,'edgecolor','none','facecolor',[1 1 1]*.35,'facelighting','phong',...,
        'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
        'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
    
    set(ax,'clim',[prctile(v,5) prctile(v,95)])
elseif flag_interp==2
    %%
    v = nan(size(MESH.Vertices_Original,1),1);
    for i=1:length(index_val)
        a = find(index_chan(i) == channel_num);
        if ~isnan(channel_num(index_chan(i)))
            v(a) = index_val(i);
        end
    end
    iiok = ~isnan(v);
    F = scatteredInterpolant(MESH.Vertices_Original(iiok,1),MESH.Vertices_Original(iiok,2),MESH.Vertices_Original(iiok,3),v(iiok));
    %
    fv.faces = MESH.Faces_OnlyVentricles;
    fv.vertices = MESH.Vertices_Original;
    index_val_interp = F([MESH.Vertices_Original(:,1),MESH.Vertices_Original(:,2),MESH.Vertices_Original(:,3)]);
    fv.facevertexcdata = index_val_interp;
    %
    fvtop.faces = MESH.Faces_Top_Interp;
    fvtop.vertices = MESH.Vertices_Interp;
    fvtop.facevertexcdata = 0;
    
    hpatch = patch(fv,'Parent', ax);
    hold(ax,'on')
    hpatch_top = patch(fvtop,'Parent', ax);
    set(hpatch,'edgecolor','none','facecolor','interp','facelighting','phong',...,
        'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
        'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
    set(hpatch_top,'edgecolor','none','facecolor',[1 1 1]*.35,'facelighting','phong',...,
        'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
        'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
    
    set(ax,'clim',[prctile(index_val_interp,5) prctile(index_val_interp,95)])
    
    
elseif flag_interp==3
    
%     F = scatteredInterpolant(MESH.Vertices_Interp(:,1),MESH.Vertices_Interp(:,2),MESH.Vertices_Interp(:,3),index_val);
%     v = F([MESH.Vertices_Interp(:,1),MESH.Vertices_Interp(:,2),MESH.Vertices_Interp(:,3)]);
    %%
    fv.faces = MESH.Faces_OnlyVentricles_Interp;
    fv.vertices = MESH.Vertices_Interp;
    v = index_val;
    index_val_interp = v;
    fv.facevertexcdata = index_val_interp;
    %
    fvtop.faces = MESH.Faces_Top_Interp;
    fvtop.vertices = MESH.Vertices_Interp;
    fvtop.facevertexcdata = 0;
    
    hpatch = patch(fv,'Parent', ax);
    hold(ax,'on')
    hpatch_top = patch(fvtop,'Parent', ax);
    set(hpatch,'edgecolor','none','facecolor','interp','facelighting','phong',...,
        'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
        'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
    set(hpatch_top,'edgecolor','none','facecolor',[1 1 1]*.35,'facelighting','phong',...,
        'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
        'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
    
    set(ax,'clim',[prctile(index_val_interp,5) prctile(index_val_interp,95)])
    
    
end
H.cbar = colorbar('peer',ax);
H.hpatch = hpatch;
H.hpatch_top = hpatch_top;

set(gcf,'CurrentAxes',ax)
for i=1:size(lad_xyz,1)-1
    hold on, 
    h_lad(i) = line(lad_xyz(i:i+1,1),lad_xyz(i:i+1,2),lad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [0 0 0]);
end
if exist('plad_xyz','var')
    for i=1:size(plad_xyz,1)-1
       hold on, 
       h_plad(i) = line(plad_xyz(i:i+1,1),plad_xyz(i:i+1,2),plad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [1 1 1]);
    end
else
     h_plad = [];
end
if ~isempty(strfind(sock_name,'new_sock4')) | ~isempty(strfind(sock_name,'old_sock4'))
    for i = 1:length(elect_name)
        if elect_name{i}(1)=='g' | elect_name{i}(1)=='K'
            col{i} = 'g';
        elseif elect_name{i}(1)=='p' | elect_name{i}(1)=='A'
            col{i} = 'm';
        elseif elect_name{i}(1)=='y' | elect_name{i}(1)=='B'
            col{i} = 'y';
        elseif elect_name{i}(1)=='r' | elect_name{i}(1)=='C'
            col{i} = 'r';
        end
    end
elseif ~isempty(strfind(sock_name,'mo_sock1'))
    for i = 1:length(elect_name)
        if elect_name{i}(1)=='w'
            col{i} = [1 1 1];
        elseif elect_name{i}(1)=='o'
            col{i} = [1 165/255 0];
        elseif elect_name{i}(1)=='y'
            col{i} = [1 1 0];
        elseif elect_name{i}(1)=='r'
            col{i} = [1 0 0];
        end
    end
elseif ~isempty(strfind(sock_name,'mo_sock2'))
    for i = 1:length(elect_name)
        if elect_name{i}(1)=='w'
            col{i} = [1 1 1];
        elseif elect_name{i}(1)=='g'
            col{i} = [0 1 0];
        elseif elect_name{i}(1)=='b'
            col{i} = [0 0 1];
        elseif elect_name{i}(1)=='r'
            col{i} = [1 0 0];
        end
    end
elseif ~isempty(strfind(sock_name,'new_sock'))
    for i = 1:length(elect_name)
        if elect_name{i}(1)=='w'
            col{i} = [1 1 1];
        elseif elect_name{i}(1)=='g'
            col{i} = [0 1 0];
        elseif elect_name{i}(1)=='y'
            col{i} = [1 1 0];
        elseif elect_name{i}(1)=='r'
            col{i} = [1 0 0];
        end
    end
else
    col(1:length(elect_name)) = {'k'};
end

H.cross = zeros(1,size(xyz,1));
H.points = zeros(1,size(xyz,1));
chanTot = [1:size(xyz,1)];
for i=1:size(xyz,1)
    if sum(chanTot==i)
        hold on
        H.points(i) = plot3(ax,xyz(i,1),xyz(i,2),xyz(i,3),'marker','.','markersize',8,'color',col{i});
        if isnan(v(i))|v(i)==-10
            % mark electrodes with nan
            H.cross(i) = plot3(ax,xyz(i,1),xyz(i,2),xyz(i,3),'marker','x','markersize',10,'color','k','linewidth',2);
        end
    end
end

H.cross(H.cross==0)=[];
H.points(H.points==0)=[];

dcm_obj = datacursormode(gcf);
if ~isempty(strfind(sock_name,'old_sock4'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mrc_oldsock4)
elseif  ~isempty(strfind(sock_name,'old_sock6'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mrc_oldsock6)
elseif  ~isempty(strfind(sock_name,'sock6_newPCB'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mrc_oldsock6_newPCB)
elseif  ~isempty(strfind(sock_name,'new_sock4'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mrc_newsock4)
elseif  ~isempty(strfind(sock_name,'new_sock'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mrc_newsock)
elseif  ~isempty(strfind(sock_name,'mo_sock1'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mo_sock1)
elseif  ~isempty(strfind(sock_name,'mo_sock2'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mo_sock2)
end

set(ax,'view',[100 -10])
H.geo.lad_xyz =lad_xyz;
H.geo.plad_xyz =plad_xyz;
H.geo.xyz =xyz;
H.ax = ax;
H.h_lad = h_lad;
H.h_plad = h_plad;

geo = load(['ALLgeoDATA_',sock_name]);
geo.Mesh = MESH;
