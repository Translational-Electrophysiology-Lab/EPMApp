function output_txt = myfunctioncursor_mrc_newsock4(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');



load('E:\UCL\Scripts_all\Scripts_mo\GUI_egm\GUI_egm_mFiles\Geo_Chann\ALLgeoDATA_new_sock4.mat')

% ie = find([pos(1)==xyz(:,1) & pos(2)==xyz(:,2) & pos(3)==xyz(:,3)]);
[~,ie] = min([abs(pos(1)-xyz(:,1)) + abs(pos(2)-xyz(:,2)) + abs(pos(3)-xyz(:,3))]);

output_txt= { [elect_name{ie}],
              ['#',num2str(channel_num(ie))]
              };

