function output_txt = myfunctioncursor_markers(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
ax = get(get(event_obj,'target'),'parent');
handles = guidata(gcf);
iic = get(handles.tag_List,'value');


if ax == handles.tag_ax_sig
    
    xx = pos(1)-handles.spikes;xx(xx<0)=[];
    [~,im] = min(xx);
    tt = (pos(1)-handles.spikes(im));
    ip = round(pos(1)/1000*handles.ParamSig.frequency);
    if get(handles.tag_show_filter,'value')
        if get(handles.tag_norm,'value')
            signal = handles.signals_proc(ip,iic)./max(abs(handles.signals_proc(:,iic)));
            jj = find(signal==pos(2));
        else
            jj = find(handles.signals_proc(ip,iic)==pos(2));
        end
    else
        if get(handles.tag_norm,'value')
            signal = handles.S(ip,iic)./max(abs(handles.S(:,iic)));
            jj = find(signal==pos(2));
        else
            jj = find(handles.S(ip,iic)==pos(2));
        end
    end
    
    if ~isempty(jj)
        output_txt= { ['IC',num2str(iic(jj(1)))],
            ['t =',num2str(tt,3),' ms']
            ['x =',num2str(pos(1)),' samp']};
    end
    
    
    
elseif ax == handles.tag_ax_spikes
    
    if get(handles.tag_show_RT,'value')&(get(handles.tag_method,'value')==2)
        rti = handles.MarkersC.rt_Wyatt(:,iic) - repmat(handles.spikes(:),[1 length(iic)]);
        [ii,jj] = find(handles.MarkersC.rt_Wyatt(:,iic)==pos(1)&rti == pos(2));
        if ~isempty(ii)
            output_txt= { ['IC',num2str(iic(jj(1)))],
                ['RT=',num2str(rti(ii(1),jj(1)),3)]
                };
        end
    end
    
    if get(handles.tag_show_RT,'value')&(get(handles.tag_method,'value')==3)
        rti = handles.MarkersC.rt_Alternative(:,iic) - repmat(handles.spikes(:),[1 length(iic)]);
        [ii,jj] = find(handles.MarkersC.rt_Alternative(:,iic)==pos(1)&rti == pos(2));
        
        if ~isempty(ii)
            output_txt= { ['IC',num2str(iic(jj(1)))],
                ['RT=',num2str(rti(ii(1),jj(1)),3)]
                };
        end
    end
    
    if get(handles.tag_show_DT,'value')
        dti = handles.MarkersC.dt(:,iic) - repmat(handles.spikes(:),[1 length(iic)]);
        [ii,jj] = find(handles.MarkersC.dt(:,iic)==pos(1)&dti == pos(2));
        
        if ~isempty(ii)
            output_txt= { ['IC',num2str(iic(jj(1)))],
                ['AT=',num2str(dti(ii(1),jj(1)),3)]
                };
        end
    end
else
    return
end

