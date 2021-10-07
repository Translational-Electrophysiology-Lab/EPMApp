function output_txt = myfunctioncursor_mrc(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
%output_txt = {['X: ',num2str(pos(1),4)],...
%    ['Y: ',num2str(pos(2),4)]};

% load('E:\UCL\Scripts_mo\MRC_updated\Geo_Chann\new_sock4_geo_new.mat')
% load('E:\UCL\Scripts_mo\MRC_updated\Geo_Chann\new_sock4_chann.mat')
load([pwd,'\Geo_Chann\new_sock4_geo_new.mat'])
load([pwd,'\Geo_Chann\new_sock4_chann.mat'])

% ie = find([pos(1)==xyz(:,1) & pos(2)==xyz(:,2) & pos(3)==xyz(:,3)]);
[~,ie] = min([abs(pos(1)-xyz(:,1)) + abs(pos(2)-xyz(:,2)) + abs(pos(3)-xyz(:,3))]);
output_txt= { [elect_name{ie}],
              ['#',num2str(channel_num(ie))]
              };

