dir_save = 'E:\UCL\Scripts_mo\GUI_egm\GUI_egm_mFiles\Geo_Chann\';
addpath E:\UCL\Scripts_mo\GUI_egm\GUI_egm_mFiles\Geo_Chann
addpath E:\UCL\Scripts_mo\GUI_egm\GUI_egm_mFiles\
do_fig = 0;
%%

sockname_tot = {'new_sock4','old_sock4','old_sock6','old_sock6_newPCB','mo_sock1','mo_sock2'};
Dist_max = 12.5;

for isocks = 6%1:length(sockname_tot)
%%
sockname = sockname_tot{isocks};
display([' * ',sockname])
filename = [dir_save,'ALLgeoDATA_',sockname,'.mat'];
load(filename)
chan_closest_chan = cell(1,240);
chan_closest_chan_D = cell(1,240);
chan_closest_xyz = cell(1,size(xyz,1));
chan_closest_xyz_D = cell(1,size(xyz,1));

for i=1:size(xyz,1)
    
    Dtot = sqrt(sum((xyz - repmat(xyz(i,:),[size(xyz,1),1])).^2,2));
    [d,id] = sort(Dtot,'ascend');
    id(1)=[]; % this is equal to i
    d(1) = [];
    xyzs = xyz(id,:);
    
    % to balance
    s = xyzs - repmat(xyz(i,:),[size(xyz,1)-1 1]);
    %     [th,ph,r] = cart2sph(s(2:10,1),s(2:10,2),s(2:10,3));
    is_diff_direction = sum(diff(sign(s)),2)~=0;
    
    chan_closest_xyz{i} = id(d(1:end-1)<=Dist_max);
    chan_closest_xyz_D{i} = d(d(1:end-1)<=Dist_max);
    iix = channel_num(chan_closest_xyz{i});
    dixx = chan_closest_xyz_D{i};
    dixx(isnan(iix)) =[];
    iix(isnan(iix)) = [];
    
    jj = channel_num(i);
    if ~isnan(jj)
        chan_closest_chan{jj} =iix;
        chan_closest_chan_D{jj} = dixx;
    end
    
    if do_fig
        figure(1),
        hold off
        ax(1) = subplot(121);
        H = surf_index_mo([1:240],[ones(1,240)],sockname,0); %
        set(H.cross,'visible','off')
        hold on,plot3(xyz(i,1),xyz(i,2),xyz(i,3),'xr','markersize',10,'linewidth',3)
        hold on,plot3(xyz(chan_closest_xyz{i},1),xyz(chan_closest_xyz{i},2),xyz(chan_closest_xyz{i},3),'xk','markersize',10,'linewidth',3)
        title(num2str(i))
        
        ax(2) = subplot(122);
        hold off,plot3(xyz(i,1),xyz(i,2),xyz(i,3),'xr','markersize',10,'linewidth',3)
        hold on,plot3(xyz(chan_closest_xyz{i},1),xyz(chan_closest_xyz{i},2),xyz(chan_closest_xyz{i},3),'xk','markersize',10,'linewidth',3)
        title(num2str(i))
        
        set(ax,'xlim',[-35 35],'ylim',[-35 35],'zlim',[0 60])
        linkprop(ax,'view')
        pause
    end
    
end

channels_closes.chan = chan_closest_chan;
channels_closes.chan_D = chan_closest_chan_D;
channels_closes.xyz = chan_closest_xyz;
channels_closes.xyz_D = chan_closest_xyz_D;
channels_closes.Dist_max = Dist_max;

channels_closes.ref_bipol = nan(length(channels_closes.chan),1);
channels_closes.ref_bipol_name = cell(length(channels_closes.chan),1);
for ic = 1:length(channels_closes.chan)
    if ~isempty(channels_closes.chan{ic})
        channels_closes.ref_bipol(ic) = channels_closes.chan{ic}(1);
        channels_closes.ref_bipol_name{ic} = elect_name{channel_num==channels_closes.chan{ic}(1)};
    end
end

display(['* Adding to ',filename])
save(filename,'channels_closes','-append')

clearvars -except sockname_tot do_fig dir_save Dist_max
end

%%
%
% sockname = 'old_sock4';
% filename = [dir_save,'ALLgeoDATA_',sockname,'.mat'];
% load(filename)
% chan_closest_chan = cell(1,240);
% chan_closest_chan_D = cell(1,240);
% chan_closest_xyz = cell(1,size(xyz,1));
% chan_closest_xyz_D = cell(1,size(xyz,1));
% for i=1:size(xyz,1)
%
%     Dtot = sqrt(sum((xyz - repmat(xyz(i,:),[size(xyz,1),1])).^2,2));
%     [d,id] = sort(Dtot,'ascend');
%     id(1)=[]; % this is equal to i
%     d(1) = [];
%     xyzs = xyz(id,:);
%
%     % to balance
%     s = xyzs - repmat(xyz(i,:),[size(xyz,1)-1 1]);
%     %     [th,ph,r] = cart2sph(s(2:10,1),s(2:10,2),s(2:10,3));
%     is_diff_direction = sum(diff(sign(s)),2)~=0;
%
%     %     chan_closest{i} = id(d(1:end-1)<=10 & is_diff_direction(:));
%     chan_closest_xyz{i} = id(d(1:end-1)<=10);
%     chan_closest_xyz_D{i} = d(d(1:end-1)<=10);
%     iix = channel_num(chan_closest_xyz{i});
%     dixx = chan_closest_xyz_D{i};
%     dixx(isnan(iix)) =[];
%     iix(isnan(iix)) = [];
%
%     jj = channel_num(i);
%     if ~isnan(jj)
%         chan_closest_chan{jj} =iix;
%         chan_closest_chan_D{jj} = dixx;
%     end
%
%     if do_fig
%         figure(1),
%         hold off
%         ax(1) = subplot(121);
%         H = surf_index_mo([1:240],[ones(1,240)],sockname,0); %
%         set(H.cross,'visible','off')
%         hold on,plot3(xyz(i,1),xyz(i,2),xyz(i,3),'xr','markersize',10,'linewidth',3)
%         hold on,plot3(xyz(chan_closest_xyz{i},1),xyz(chan_closest_xyz{i},2),xyz(chan_closest_xyz{i},3),'xk','markersize',10,'linewidth',3)
%         title(num2str(i))
%
%         ax(2) = subplot(122);
%         hold off,plot3(xyz(i,1),xyz(i,2),xyz(i,3),'xr','markersize',10,'linewidth',3)
%         hold on,plot3(xyz(chan_closest_xyz{i},1),xyz(chan_closest_xyz{i},2),xyz(chan_closest_xyz{i},3),'xk','markersize',10,'linewidth',3)
%         title(num2str(i))
%
%         set(ax,'xlim',[-35 35],'ylim',[-35 35],'zlim',[0 60])
%         linkprop(ax,'view')
%         pause
%     end
% end
% channels_closes.chan = chan_closest_chan;
% channels_closes.chan_D = chan_closest_chan_D;
% channels_closes.xyz = chan_closest_xyz;
% channels_closes.xyz_D = chan_closest_xyz_D;
% display(['* Adding to ',filename])
% save(filename,'channels_closes','-append')

%%



