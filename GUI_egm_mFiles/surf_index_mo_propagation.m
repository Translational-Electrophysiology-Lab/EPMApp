%% combination of surf_index_interp.m and surf_index_v3

function [H,index_val_interp] = surf_index_mo_propagation(index_chan,index_val,sock_name,type_flag,intervals);
%
%% Ex (1): surf_index_mo_propagation([1:length(dataDT)],dataDT,'new_sock4',1,6);
%% Ex (2): surf_index_mo_propagation([1:length(dataDT)],dataDT,sockname,2,[0 8 20 60 66 83 100]);

%% INPUT
% index_chan :  channels for which we have a value 
% index_val :   value at a given channel
% sock_name:    name of the sock
% type_flag:    ==0: no 3D interpolation (raw data);
%               ==1 3D interpolation (channel for which there is a nan are assigned a value by interpolation);
%               ==2 no 3D interpolation but a better representation of the sock
% intervals:        (a) Number of frames (scalar) or (b) percentiles in which index_val (or its interpolated version) will be divided 

%%
addpath .\GUI_egm_mFiles\Geo_Chann
addpath .\Geo_Chann\
load(['ALLgeoDATA_',sock_name])
%%



if type_flag==1
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
%     F = scatteredInterpolant(x(iiok),y(iiok),z(iiok),v(iiok));
    F = scatteredInterpolant(x(iiok),y(iiok),z(iiok),v(iiok),'natural');

    vtot = F([x,y,z]);
    index_val_interp = nan(1,240);
    index_val_interp(channel_num(~isnan(channel_num))) = vtot(~isnan(channel_num));
    
    surf_c = vtot(closest_index);
    surf_c(isnan(surf_c))=-10;
    surf_c = reshape(surf_c(:),size(surf_x,1),size(surf_x,2));
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
    
    %%
    if length(intervals)==1
        int_time = prctile(surf_c2(:),linspace(0,100,intervals+1));
        Nplot = intervals;
    else
        int_time = prctile(surf_c2(:),intervals);
        Nplot = length(intervals)-1;
    end
    
    figure
    for jj = 1:Nplot
        D = surf_c2;
%         D([surf_c2>=int_time(jj)&surf_c2<=int_time(jj+1)]) = 10;
%         D([surf_c2<int_time(jj)|surf_c2>int_time(jj+1)]) = 0;
        D([surf_c2<=int_time(jj+1)]) = 10;
        D([surf_c2>int_time(jj+1)]) = 0;
        D(isnan(D))=-10;
        ax(jj)=subplot(1,Nplot,jj);
        surf(ax(jj),surf_x,surf_y,surf_z,D,'FaceColor','flat','EdgeColor','none','FaceLighting','phong'),
        title(['[',num2str(int_time(jj),3),' - ',num2str(int_time(jj+1),3),']'])
        for i=1:size(lad_xyz,1)-1
            hold on,
            H.lad(jj) = line(lad_xyz(i:i+1,1),lad_xyz(i:i+1,2),lad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [0 0 0]);
        end
        for i=1:size(plad_xyz,1)-1
            hold on,
            H.plad(jj) = line(plad_xyz(i:i+1,1),plad_xyz(i:i+1,2),plad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [1 1 1]);
        end
        
    end
    colormap([.7 .7 .7;0 0 1;1 0 0])
    shading interp
    %%
elseif type_flag==0
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
    %     surf_c(isnan(surf_c))=-10;
    surf_c = reshape(surf_c(:),size(surf_x,1),size(surf_x,2));
    
    % -
   if length(intervals)==1
        int_time = prctile(surf_c(:),linspace(0,100,intervals+1));
        Nplot = intervals;
    else
        int_time = prctile(surf_c(:),intervals);
        Nplot = length(intervals)-1;
    end
    
    figure
    for jj = 1:Nplot
        D = surf_c;
        D([surf_c>=int_time(jj)&surf_c<=int_time(jj+1)]) = 10;
        D([surf_c<int_time(jj)|surf_c>int_time(jj+1)]) = 0;
        D(isnan(D))=-10;
        ax(jj)=subplot(1,Nplot,jj);
        surf(ax(jj),surf_x,surf_y,surf_z,D,'FaceColor','flat','EdgeColor','none','FaceLighting','phong'),
        title(['[',num2str(int_time(jj),3),' - ',num2str(int_time(jj+1),3),']'])
        for i=1:size(lad_xyz,1)-1
            hold on,
            H.lad(jj) = line(lad_xyz(i:i+1,1),lad_xyz(i:i+1,2),lad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [0 0 0]);
        end
        for i=1:size(plad_xyz,1)-1
            hold on,
            H.plad(jj) = line(plad_xyz(i:i+1,1),plad_xyz(i:i+1,2),plad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [1 1 1]);
        end
    end
    colormap([.7 .7 .7;0 0 1;1 0 0])
    shading interp
    
%%    
elseif type_flag==2
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
%     surf_c(isnan(surf_c))=-10;
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
    
   if length(intervals)==1
        int_time = prctile(surf_c2(:),linspace(0,100,intervals+1));
        Nplot = intervals;
    else
        int_time = prctile(surf_c2(:),intervals);
        Nplot = length(intervals)-1;
    end
    
    figure
    for jj = 1:Nplot
        D = surf_c2;
        D([surf_c2>=int_time(jj)&surf_c2<=int_time(jj+1)]) = 10;
        D([surf_c2<int_time(jj)|surf_c2>int_time(jj+1)]) = 0;
        D(isnan(D))=-10;
        ax(jj)=subplot(1,Nplot,jj);
        surf(ax(jj),surf_x,surf_y,surf_z,D,'FaceColor','flat','EdgeColor','none','FaceLighting','phong'),
        title(['[',num2str(int_time(jj),3),' - ',num2str(int_time(jj+1),3),']'])
        for i=1:size(lad_xyz,1)-1
            hold on, 
            H.lad(jj) = line(lad_xyz(i:i+1,1),lad_xyz(i:i+1,2),lad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [0 0 0]);
        end
        for i=1:size(plad_xyz,1)-1
            hold on, 
            H.plad(jj) = line(plad_xyz(i:i+1,1),plad_xyz(i:i+1,2),plad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [1 1 1])
        end
    end
    colormap([.7 .7 .7;0 0 1;1 0 0])
    shading interp
    
    %     surf_c2(isnan(surf_c2))=-10;
    %     % ----
    %     hold(ax,'off')
    %     surf(ax,surf_x,surf_y,surf_z,surf_c2,'FaceColor','flat','EdgeColor','none','FaceLighting','phong'),
    %     H.cbar = colorbar('peer',ax);
    %
end

% set(gcf,'CurrentAxes',ax)
% for i=1:size(lad_xyz,1)-1
%     hold on, line(lad_xyz(i:i+1,1),lad_xyz(i:i+1,2),lad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [0 0 0])
% end
% for i=1:size(plad_xyz,1)-1
%     hold on, line(plad_xyz(i:i+1,1),plad_xyz(i:i+1,2),plad_xyz(i:i+1,3),'LineWidth',5.5,'Color', [1 1 1])
% end
%
% if ~isempty(strfind(sock_name,'new_sock4')) | ~isempty(strfind(sock_name,'old_sock4'))
%     for i = 1:length(elect_name)
%         if elect_name{i}(1)=='g' | elect_name{i}(1)=='K'
%             col{i} = 'g';
%         elseif elect_name{i}(1)=='p' | elect_name{i}(1)=='A'
%             col{i} = 'm';
%         elseif elect_name{i}(1)=='y' | elect_name{i}(1)=='B'
%             col{i} = 'y';
%         elseif elect_name{i}(1)=='r' | elect_name{i}(1)=='C'
%             col{i} = 'r';
%         end
%     end
% elseif ~isempty(strfind(sock_name,'mo_sock1'))
%     for i = 1:length(elect_name)
%         if elect_name{i}(1)=='w'
%             col{i} = [1 1 1];
%         elseif elect_name{i}(1)=='o'
%             col{i} = [1 165/255 0];
%         elseif elect_name{i}(1)=='y'
%             col{i} = [1 1 0];
%         elseif elect_name{i}(1)=='r'
%             col{i} = [1 0 0];
%         end
%     end
% elseif ~isempty(strfind(sock_name,'mo_sock2'))
%     for i = 1:length(elect_name)
%         if elect_name{i}(1)=='w'
%             col{i} = [1 1 1];
%         elseif elect_name{i}(1)=='g'
%             col{i} = [0 1 0];
%         elseif elect_name{i}(1)=='b'
%             col{i} = [0 0 1];
%         elseif elect_name{i}(1)=='r'
%             col{i} = [1 0 0];
%         end
%     end
% else
%     col(1:length(elect_name)) = {'k'};
% end

% chanTot = [1:size(xyz,1)];
% for i=1:size(xyz,1)
%     if sum(chanTot==i)
%         hold on
%         H.points(i) = plot3(ax,xyz(i,1),xyz(i,2),xyz(i,3),'marker','.','markersize',8,'color',col{i});
%         if isnan(v(i))
%             % mark electrodes with nan
%             H.cross(i) = plot3(ax,xyz(i,1),xyz(i,2),xyz(i,3),'marker','x','markersize',10,'color','k','linewidth',2);
%         end
%     end
% end

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
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mrc_new_sock)
elseif  ~isempty(strfind(sock_name,'mo_sock1'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mo_sock1)
elseif  ~isempty(strfind(sock_name,'mo_sock2'))
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_mo_sock2)
end

H.ax = ax;
hlink = linkprop(ax,'view')
% hlink = getappdata(gca,'graphics_linkprop');
setappdata(H.ax(end),'graphics_linkprop',hlink);

L = nan;
if Nplot == 5
    L = 0.15; W = 0.88; DL = 0.05;
    set(gcf,'units','centimeters','paperUnits','centimeters','position',[2  2  18  7],'papersize',[18 7],'paperposition',[0 0 18 7])
elseif Nplot == 4
    L = 0.2; W = 0.88; DL = 0.05;
    set(gcf,'units','centimeters','paperUnits','centimeters','position',[2  2  14  7],'papersize',[14 7],'paperposition',[0 0 16 7])
elseif Nplot ==3
    L = 0.25; W = 0.88; DL = 0.05;
    set(gcf,'units','centimeters','paperUnits','centimeters','position',[2  2  14  7],'papersize',[14 7],'paperposition',[0 0 14 7])
elseif Nplot ==6
    L = 0.12; W = 0.88; DL = 0.04;
    set(gcf,'units','centimeters','paperUnits','centimeters','position',[2  2  18  7],'papersize',[18 7],'paperposition',[0 0 18 7])
end

if ~isnan(L)
    for i = 1:Nplot
        axis(ax(i),'tight','off')
        set(ax(i),'position',[.025+(L+DL)*(i-1) 0.05 L W])
    end
end


