function output_txt = myfunctioncursor_markers_CARTO(~,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).



t = get(event_obj,'Target');
ax = get(t,'parent');
handles = guidata(get(ax,'parent'));
ib = str2double(get(handles.tag_beat_number,'string'));
iic = handles.matfile.ParamSig.Label;
%
%
if ax==handles.tag_mesh;
    delete(findall(handles.tag_egm,'Type','hggroup'));
%         bb = findobj(handles.tag_mesh,'markersize',18,'linewidth',1,'color','k','marker','x');
%         delete(bb);
%     
        bb = findobj(handles.tag_mesh,'marker','x');
        delete(bb);
        
    %     aa = findobj(handles.tag_mesh,'tag','datatip');
    %     delete(aa)
    %     clear aa
    %         output_txt= {['x=',num2str(event_obj.Position(1))]};
    
    
    pos = get(event_obj,'Position');
    
    %     if ~handles.tag_project.Value
    
    [~,im] = min(sum(abs(handles.matfile.geo.xyz  - repmat(pos,[size(handles.matfile.geo.xyz,1),1])),2));
    %     if D>5
    %         disp('no point found');
    %         return
    %     end
    dt = handles.matfile.Markers.dt(ib,im) - handles.matfile.spikes(ib);
    rt = handles.matfile.Markers.rt_Wyatt(ib,im) - handles.matfile.spikes(ib);
    
    
    
    
    hold(handles.tag_mesh,'on');
    plot3(handles.tag_mesh,handles.matfile.geo.xyz(im,1),handles.matfile.geo.xyz(im,2),handles.matfile.geo.xyz(im,3),'xk','linewidth',1,'markersize',18);
    
    output_txt= { ['Node #',num2str(im)],['AT=',num2str(dt,'%1.0f'),' ms'],['ARI=',num2str(rt-dt,'%1.0f'),' ms'],['RT=',num2str(rt,'%1.0f'),' ms'],
        };
    %     end
    
    
elseif get(get(event_obj,'Target'),'parent')==handles.tag_egm;
    delete(findall(handles.tag_mesh,'Type','hggroup'));
    pos = get(event_obj,'Position');
    samp = round([pos(1)+handles.matfile.spikes(1)]/1000*handles.matfile.ParamSig.frequency);
    %     im = find(handles.matfile.signals(samp,:)==pos(2));
    if handles.tag_menu.Value == 1
        if isfield(handles.matfile,'signals_proc_Ensite_AT')
            S = handles.matfile.signals_proc_Ensite_AT(samp,:);
        else
            S = handles.matfile.signals_proc(samp,:);
        end
    elseif handles.tag_menu.Value == 2
        if isfield(handles.matfile,'signals_proc_Ensite_RT')
            S = handles.matfile.signals_proc_Ensite_RT(samp,:);
        else
            S = handles.matfile.signals_proc(samp,:);
        end
    else
        S = handles.matfile.signals(samp,:);
    end
    
    im = find(S==pos(2));
    
    if ~isempty(im)
        t0 =  pos(1)- handles.matfile.spikes(ib) + handles.matfile.spikes(1);
        dt = handles.matfile.Markers.dt(ib,im) - handles.matfile.spikes(ib);
        rt = handles.matfile.Markers.rt_Wyatt(ib,im) - handles.matfile.spikes(ib);
        
        output_txt= { ['Node #',num2str(im)],['t=',num2str(t0,'%1.0f'),' ms'],['AT=',num2str(dt,'%1.0f'),' ms'],['ARI=',num2str(rt-dt,'%1.0f'),' ms'],['RT=',num2str(rt,'%1.0f'),' ms'],
            };
    else
        t0 =  pos(1)- handles.matfile.spikes(ib) + handles.matfile.spikes(1);
        output_txt= {['x = ',num2str(pos(1),'%1.0f'),' ms'];['t = ',num2str(t0,'%1.0f'),' ms']};
    end
    
end

