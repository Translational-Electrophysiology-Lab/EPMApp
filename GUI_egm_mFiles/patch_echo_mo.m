function [H,index_val_interp] = patch_echo_mo(index_val,sock_name,nanvalues);

if nargin<3
    nanvalues = min(index_val)-abs(min(index_val)/10);
end
% function [H,index_val_interp] = patch_echo_mo(var1,var2,var3,var4,var5,var6);
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
% if ishandle(var1)
%     ax = var1;
%     index_chan = var2;
%     index_val = var3;
%     sock_name = var4;
%     if nargin>4
%         flag_interp = var5;
%     else
%         flag_interp = 0;
%     end
%     if nargin==6
%         allow_neg_val = var6;
%     else
%         allow_neg_val = 0;
%     end
% else
%     ax = gca;
%     index_chan = var1;
%     index_val = var2;
%     sock_name = var3;
%     if nargin>3
%         flag_interp = var4;
%     else
%         flag_interp = 0;
%     end
%     if nargin==5
%         allow_neg_val = var5;
%     else
%         allow_neg_val = 0;
%     end
% end
% clear var*
%%
ax = gca;
if ~isequal(sock_name,'mo_sock2')
    error('Only with sock geo = mo_sock2')
end
if length(index_val)~=17
   error('index_val mush be a vector of 17 elements') 
end

addpath E:\UCL\Scripts_mo\GUI_egm\GUI_egm_mFiles\Geo_Chann
load(['ALLgeoDATA_',sock_name,'_mesh'])
load(['ALLgeoDATA_',sock_name]);
%%
ii_septal = [2,3,8,9,14];
values = nan(256,1);
for i = 1:length(index_val)
   if sum(ii_septal==i)==0
       values(Bullseye{i})=index_val(i);
   else
       values(Bullseye{i})=nanvalues;
   end
end
iiok = ~isnan(values);
F = scatteredInterpolant(MESH.Vertices_Original(iiok,1),MESH.Vertices_Original(iiok,2),MESH.Vertices_Original(iiok,3),values(iiok));
fv.faces = MESH.Faces_OnlyVentricles_Interp;
fv.vertices = MESH.Vertices_Interp;
index_val_interp = F([MESH.Vertices_Interp(:,1),MESH.Vertices_Interp(:,2),MESH.Vertices_Interp(:,3)]);
fv.facevertexcdata = index_val_interp;
%
fvtop.faces = MESH.Faces_Top_Interp;
fvtop.vertices = MESH.Vertices_Interp;
fvtop.facevertexcdata = 0;

% figure
hpatch = patch(fv);
hpatch_top = patch(fvtop);
set(hpatch,'edgecolor','none','facecolor','interp','facelighting','phong',...,
    'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
    'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
set(hpatch_top,'edgecolor','none','facecolor',[1 1 1]*.35,'facelighting','phong',...,
    'edgelighting','none','ambientstrength',0.6,'diffusestrength',0.6,...,
    'specularstrength',0.2,'specularexponent',5,'specularcolorreflectance',0.5,'marker','none','markersize',6,'markerfacecolor','k');
H.cbar = colorbar('peer',ax);
H.hpatch = hpatch;
H.hpatch_top = hpatch_top;

set(gcf,'CurrentAxes',ax)
for i=1:size(lad_xyz,1)-1
    hold on, line(lad_xyz(i:i+1,1),lad_xyz(i:i+1,2),lad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [0 0 0])
end
if exist('plad_xyz','var')
    for i=1:size(plad_xyz,1)-1
        hold on, line(plad_xyz(i:i+1,1),plad_xyz(i:i+1,2),plad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [1 1 1])
    end
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

chanTot = [1:size(xyz,1)];
for i=1:size(xyz,1)
    if sum(chanTot==i)
        hold on
        H.points(i) = plot3(ax,xyz(i,1),xyz(i,2),xyz(i,3),'marker','.','markersize',8,'color',col{i});
        if isnan(values(i))|values(i)==nanvalues
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
H.ax = ax;