%% combination of surf_index_interp.m and surf_index_v3

function [H,index_val_interp] = surf_index_mo(var1,var2,var3,var4,var5,var6);
%% plot surface
%% If first argument is an axis handle:
% surf_index_mo([1:240],xx(:,k)-min(x(:,k)),'mo_sock2',1,0);
% var1 : number of channels to plot
% var2 : values of channels 
% var3 : sock name
% var4 : interpolation; 0=no;1=yes;2=no but better representation of sock 
% var5 : allow negative values; 0=no [default]; 1=yes
%% If first argument is NOT an axis handle:
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
% addpath .\GUI_egm_mFiles\Geo_Chann
% addpath .\Geo_Chann\
addpath E:\UCL\Scripts_all\Scripts_mo\GUI_egm\GUI_egm_mFiles\Geo_Chann
addpath E:\UCL\Scripts_all\Scripts_mo\GUI_egm\GUI_egm_mFiles\
load(['ALLgeoDATA_',sock_name])
%%


if flag_interp==1 
    %%
    x = xyz(:,1);
    y = xyz(:,2);
    z = xyz(:,3);
    v = nan(size(x));
    for i=1:length(index_val)
        a = find(index_chan(i) == channel_num);
        if ~isnan(channel_num(index_chan(i)))
            v(a) = index_val(i);
        end
    end
    iiok = ~isnan(v);
    F = scatteredInterpolant(x(iiok),y(iiok),z(iiok),v(iiok));
    vtot = F([x,y,z]);
    F2 = scatteredInterpolant(x,y,z,vtot,'natural');
    vtot2 = F2([surf_x(:),surf_y(:),surf_z(:)]);

    if ~allow_neg_val
        vtot2(vtot2<0) = 0;
    end
    index_val_interp = nan(1,240);
    index_val_interp(channel_num(~isnan(channel_num))) = vtot(~isnan(channel_num));
    
    surf_c = reshape(vtot2(:),size(surf_x,1),size(surf_x,2));

%     surf_c = vtot(closest_index);
%     surf_c(isnan(surf_c))=-10;
%     surf_c = reshape(surf_c(:),size(surf_x,1),size(surf_x,2));
    %%
    surf_x2 = surf_x;
    surf_y2 = surf_y;
    surf_z2 = surf_z;
    surf_c2 = surf_c;
    if sum(isnan(surf_x(:)))==0
        for i=1:size(surf_x,2)
            P = [surf_x(end,i) surf_y(end,i) surf_z(end,i)];
            D = sum( (xyz - repmat(P,[size(xyz,1),1])).^2,2);
            [~,im] = min(D);
            zth(i) = xyz(im,3)+1;
        end
        z = [zth zth zth];
        L = 20;
        zthsm = filtfilt(hamming(L),sum(hamming(L)),z);
        zthsm = zthsm(length(zth)+1:end-length(zth));
        for i=1:size(surf_x,2)
            iiko = surf_z(:,1)>zthsm(i);
            surf_x2(iiko,i) = nan;
            surf_y2(iiko,i) = nan;
            surf_z2(iiko,i) = nan;
            surf_c2(iiko,i) = nan;
        end
    end
    surf_c2(isnan(surf_c2))=-10;
    % ----
    hold(ax,'off')
    surf(ax,surf_x,surf_y,surf_z,surf_c2,'FaceColor','flat','EdgeColor','none','FaceLighting','phong'),
    H.cbar = colorbar('peer',ax);
    
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
    surf_c = v(closest_index);
    surf_c(isnan(surf_c))=-10;
    surf_c = reshape(surf_c(:),size(surf_x,1),size(surf_x,2));
    hold(ax,'off')
    surf(ax,surf_x,surf_y,surf_z,surf_c,'FaceColor','interp','EdgeColor','none','FaceLighting','phong'),

elseif flag_interp==2
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
    surf_c = v(closest_index);
    surf_c(isnan(surf_c))=-10;
    surf_c = reshape(surf_c(:),size(surf_x,1),size(surf_x,2));
    %
    surf_x2 = surf_x;
    surf_y2 = surf_y;
    surf_z2 = surf_z;
    surf_c2 = surf_c;
    if sum(isnan(surf_x(:)))==0
        for i=1:size(surf_x,2)
            P = [surf_x(end,i) surf_y(end,i) surf_z(end,i)];
            D = sum( (xyz - repmat(P,[size(xyz,1),1])).^2,2);
            [~,im] = min(D);
            zth(i) = xyz(im,3)+1;
        end
        z = [zth zth zth];
        L = 20;
        zthsm = filtfilt(hamming(L),sum(hamming(L)),z);
        zthsm = zthsm(length(zth)+1:end-length(zth));
        for i=1:size(surf_x,2)
            iiko = surf_z(:,1)>zthsm(i);
            surf_x2(iiko,i) = nan;
            surf_y2(iiko,i) = nan;
            surf_z2(iiko,i) = nan;
            surf_c2(iiko,i) = nan;
        end
    end
    surf_c2(isnan(surf_c2))=-10;
    % ----
    hold(ax,'off')
    surf(ax,surf_x,surf_y,surf_z,surf_c2,'FaceColor','flat','EdgeColor','none','FaceLighting','phong'),
    H.cbar = colorbar('peer',ax);

elseif flag_interp==3 
    %%
    x = xyz(:,1);
    y = xyz(:,2);
    z = xyz(:,3);
    v = nan(size(x));
    for i=1:length(index_val)
        a = find(index_chan(i) == channel_num);
        if ~isnan(channel_num(index_chan(i)))
            v(a) = index_val(i);
        end
    end
    iiok = ~isnan(v);
    F = scatteredInterpolant(x(iiok),y(iiok),z(iiok),v(iiok));
    vtot = F([x,y,z]);
    
    if ~allow_neg_val
        vtot2(vtot<0) = 0;
    end
    index_val_interp = nan(1,240);
    index_val_interp(channel_num(~isnan(channel_num))) = vtot(~isnan(channel_num));
    
    x = interp1([1:length(x)],x,linspace(1,length(x),2000));
    y = interp1([1:length(y)],y,linspace(1,length(y),2000));
    z = interp1([1:length(z)],z,linspace(1,length(z),2000));
    %     F2 = scatteredInterpolant(x(:),y(:),z(:),vtot(:),'natural');
    vtot = F([x(:),y(:),z(:)]);

    tri = delaunay(x,y,z);
    trisurf(tri,x,y,z,vtot, 'LineStyle', 'none');shading interp
    H.cbar = colorbar('peer',ax);

end
H.cbar = colorbar('peer',ax);

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

chanTot = [1:size(xyz,1)];
for i=1:size(xyz,1)
    if sum(chanTot==i)
        hold on
        H.points(i) = plot3(ax,xyz(i,1),xyz(i,2),xyz(i,3),'marker','.','markersize',8,'color',col{i});
        if isnan(v(i))
            % mark electrodes with nan
            H.cross(i) = plot3(ax,xyz(i,1),xyz(i,2),xyz(i,3),'marker','x','markersize',10,'color','k','linewidth',2);
        end
    end
end

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

H.ax = ax;
H.h_lad = h_lad;
H.h_plad = h_plad;