function [Markers2] = remove_all_markers(Markers,ii_chan_ko,beat_n);

if nargin<2
    ii_chan_ko = 1:size(Markers.dt,2);
end

if nargin<3
    beat_n = 1:size(Markers.dt,1);
end
Markers2 = Markers;

vv = fieldnames(Markers);
vv(strcmp(vv,'iiTwpos'))=[];
for ic = 1:length(vv);
    if ~isstruct(Markers.(vv{ic})) & ~iscell(Markers.(vv{ic}))
    Markers2.(vv{ic})(beat_n,ii_chan_ko) = nan;
    end
end