function output_txt = myfunctioncursor_segment_remove(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

t = get(event_obj,'Target');
ax = get(t,'parent');
x = event_obj.Position;
hold(ax,'on');
aa = plot3(ax,x(1),x(2),x(3),'or','markerfacecolor','k');
aa.Tag = 'segment';

cc = findobj(ax,'tag','segment');
output_txt= { ['X']};
    