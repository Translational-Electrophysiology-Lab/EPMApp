function [M,H] = video_sock_ECHO_beat2beat_fun(Values,Beats,Clim,Title,CbarLab,View0,do_annotation,do_save_pics);

if nargin<7
    do_annotation = 0;
end
if nargin<6
    View0 = 120;
end
if nargin<5
    CbarLab = '';
end
if nargin<4
    Title = '';
end
if nargin<3
    Clim = [];
end
if nargin<2
    Beats = [1:size(Values,1)]';
end
%%
addpath E:\UCL\Scripts_mo\GUI_egm\GUI_egm_mFiles
addpath E:\UCL\Scripts_mo\GUI_egm\GUI_egm_mFiles\Geo_Chann\
load AllgeoDATA_mo_sock2_mesh
load AllgeoDATA_mo_sock2
[~,ia,ib] = intersect(channel_num,[1:240]);

if ~isnan(Clim)
    if numel(Clim)==2
        Clim = repmat(Clim(:)',[size(Values,1) 1]);
    end
else
    Clim = [prctile(Values,5,2) prctile(Values,95,2)];
end

if size(Beats,2)>2
    nbeats = Beats(:);
    tsec = [];
elseif size(Beats,2)==2
    nbeats = Beats(:,1);
    tsec = Beats(:,2);
elseif numel(Beats)==1
    nbeats = [];
    tsec = [];
else
    nbeats = Beats(:,1);
    tsec = [];
end

do_interp =1;
figure
set(gcf,'units','pixels','position',[50 50 800 300],'paperunits','centimeters','paperposition',[0 0 21 8.5],'papersize',[21 8.5],'InvertHardCopy','off','color','white')
% load colormap_carto
% colormap(cmap)
cmap = jet;
cmap(1,:) = [.8 .8 .8];
colormap(cmap)

nanvalues = min(Values(:))-abs(min(Values(:))/10);
        
for n = 1:size(Values,1)
    values = Values(n,:);
    
    if ~isempty(nbeats)
        Tb = ['#Beat = ',num2str(nbeats(n))];
    else
        Tb = [];
    end
    if isempty(tsec)
        Tt = [];
    else
        Tt= ['t = ',num2str(round(tsec(n))),'sec'];
    end
    
    title_plot = [Title,'  ',Tb,'  ',Tt];
    aa=annotation('textbox',[.01 .9 .9 .1],'string',title_plot,'linestyle','none','fontsize',14);
    %%
    if n==1
        for ipos = 1:5
            ax(ipos) = subplot(1,6,ipos);
            [H(ipos),data_interp] = patch_echo_mo(values(:),'mo_sock2',nanvalues);
            delete(H(ipos).cross)
            delete(H(ipos).points)
            set(H(ipos).cbar,'visible','off')
        end
        
        set(H(5).cbar,'visible','on','position',[.90 .1 .02 .8])
        xlabel(H(end).cbar,CbarLab,'fontsize',12)
        for i = 1:4
            set(ax(i),'view',[View0+360/4*(i-1) -10])
        end
        set(ax(end),'view',[90 -90])
        %
        set(ax(1:4),'xlim',[-30 30],'ylim',[-30 30],'zlim',[0 66])
        set(ax(5),'xlim',[-35 35],'ylim',[-35 35],'zlim',[0 66])
        
        WW = 0.18;DW = 0.002;HH=0.8;
        set(ax(1),'position',[0.01 0.1 WW HH])
        set(ax(2),'position',[0.01+(WW+DW)-0.01 0.1 WW HH])
        set(ax(3),'position',[0.01+2*(WW+DW) 0.1 WW HH])
        set(ax(4),'position',[0.01+3*(WW+DW)-0.01 0.1 WW HH])
        set(ax(5),'position',[0.01+4*(WW+DW)+DW*2 0.15 WW 0.6])
        set(H(5).cbar,'position',[.95 0.15 0.02 .7])
        
        axis(ax,'off')
        if do_annotation
            A(1) = annotation(gcf,'textbox',[.15 .05 .1 .1],'string','LV','horizontalalignment','center','linestyle','none');
            arr(1) = annotation(gcf,'arrow','x',[.18 .18-.02],'y',[.16 .32],'HeadStyle','none');
            arr(2) = annotation(gcf,'arrow','x',[.18 .18+.02],'y',[.16 .32],'HeadStyle','none');
            
            A(2) = annotation(gcf,'textbox',[.52 .05 .1 .1],'string','RV','horizontalalignment','center','linestyle','none');
            arr(3) = annotation(gcf,'arrow','x',[.55 .55-0.02],'y',[.16 .32],'HeadStyle','none');
            arr(4) = annotation(gcf,'arrow','x',[.55 .55+0.02],'y',[.16 .32],'HeadStyle','none');
            
            A(3)= annotation(gcf,'textbox',[.85 .75 .1 .1],'string','LV','horizontalalignment','center','linestyle','none');
            A(4) = annotation(gcf,'textbox',[.77 .75 .1 .1],'string','RV','horizontalalignment','center','linestyle','none');
            set(A,'fontsize',13)
        end
        
    else
        ii_septal = [2,3,8,9,14];
        v = nan(256,1);
        for i = 1:length(values)
            if sum(ii_septal==i)==0
                v(Bullseye{i})=values(i);
            else
                v(Bullseye{i})=nanvalues;
            end
        end
        iiok = ~isnan(v);
        F = scatteredInterpolant(MESH.Vertices_Original(iiok,1),MESH.Vertices_Original(iiok,2),MESH.Vertices_Original(iiok,3),v(iiok));
        fv.faces = MESH.Faces_OnlyVentricles_Interp;
        fv.vertices = MESH.Vertices_Interp;
        data_interp = F([MESH.Vertices_Interp(:,1),MESH.Vertices_Interp(:,2),MESH.Vertices_Interp(:,3)]);
        % fv.facevertexcdata = index_val_interp;
        % %
        % fvtop.faces = MESH.Faces_Top_Interp;
        % fvtop.vertices = MESH.Vertices_Interp;
        % fvtop.facevertexcdata = 0;
        
        
        
        %         data_interp = F([MESH.Vertices_Interp(:,1),MESH.Vertices_Interp(:,2),MESH.Vertices_Interp(:,3)]);
        po = findobj(gcf,'type','patch','facecolor','interp');
        set(po,'FaceVertexCData',data_interp)
        
    end
    set(ax,'clim',[Clim(n,1) Clim(n,2)])
    
    M(n) = getframe(gcf);
    

    if do_save_pics
        if n==1
            RR = num2str(randi(100));
        end
        display(['- Saving as: ',pwd,'\video',RR,'_pic_N_',num2str(n)])
        print(gcf,['video',RR,'_pic_N_',num2str(n)],'-djpeg','-r600')
    end
    
    if n~=size(Values,1)
        delete(aa)
    end
end


