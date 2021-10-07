function myfun_CellSelectionCallback_ATRT(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
%handles=guidata(source);



ii = eventdata.Indices;
A = guidata(hObject);

ii_manually_deleted = [A.tag_table.Data{2:end,5}];
ii_manually_checked = [A.tag_table.Data{2:end,6}];

ichan = A.dat{ii(1),2};
ff = findobj(0,'name','check channels');
if isempty(ff)
    hfig = figure('name','check channels');
    hfig.Children = gca;
else
    hfig = ff(1);
end
clear ff ff2
cc = hfig.Children;
if isempty(cc)
    ax = axes(hfig);
else
    ax = cc(1);
end
hold(ax,'off');
aa = plot(ax,A.data.S(:,ichan));
xlabel(ax,'Time (samples)');
ylabel(ax,'mV');

tt = A.data.ParamSig.Label{ichan};tt(tt=='_')='-';
if ii_manually_deleted(ichan);
    title(ax,[tt,' DELETED']);
    ax.XColor = 'r';
    ax.YColor = 'r';
    aa.Color = 'k';
else
    title(ax,[tt,' OK']);
end


set(A.tag_n_bad_chan,'string',num2str(sum(ii_manually_deleted)));
set(A.tag_n_ok_chan,'string',num2str(sum(~ii_manually_deleted)));
set(A.tag_n_check,'string',num2str(sum(ii_manually_checked)));

hold(ax,'on');
if isfield(A.data,'Markers');
    x = round(A.data.Markers.dt(:,ichan)/1000*A.data.ParamSig.frequency);
    x(isnan(x)|x>size(A.data.S,1)) = [];
    y = A.data.S(x,ichan);
    plot(ax,x,y,'o','markerfacecolor',[.2 .2 .9],'color','k');
    
    x = round(A.data.Markers.rt_Wyatt(:,ichan)/1000*A.data.ParamSig.frequency);
    x(isnan(x)) = [];
    y = A.data.S(x,ichan);
    plot(ax,x,y,'x','markerfacecolor',[.2 .2 .8],'color','k','marker','x','markersize',10,'linewidth',2);
end

xlim(ax,[0 size(A.data.S,1)])
hold(ax,'off')

