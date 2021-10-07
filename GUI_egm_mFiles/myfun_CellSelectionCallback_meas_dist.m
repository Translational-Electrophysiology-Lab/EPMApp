function myfun_CellSelectionCallback_meas_dist(hObject, eventdata)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
%handles=guidata(source);



ii = eventdata.Indices;
if size(ii,1)==0
    return
end
if ii(2)==4
    return
end


ff = findobj(0,'name','GUI_measure_distance');
if ~isempty(ff)
    hfig = ff(1);
else
    return
end
hD = guidata(hfig);
hD.ii_select = [eventdata.Indices];
guidata(hfig,hD)
