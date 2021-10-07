[filename, pathname] = uigetfile('*.txt','Select a txt file from Bard');
filename = [pathname,filename];
do_only_header = 1;
[~,Param] = LoadBardSignals_moCathLab(filename,do_only_header);

[filenameMAT, pathnameMAT] = uigetfile('*.mat','Select a txt file from Bard',filename);
filenameMAT = [pathnameMAT,filenameMAT];

V = load([filenameMAT],'signals*');

param_all = {'signals_raw','signals','signals_proc'};
for jsig = 1:length(param_all)
    
    if isfield(V,param_all{jsig})
        x = nan(size(V.(param_all{jsig})));
        if size(V.signals_raw,2)==length(Param.range)
            
            for i = 1:length(Param.range)
                Aold = power(10,Param.scale(i))*1000;
                x(:,i) = V.(param_all{jsig})(:,i)/Aold/(2^15)*Param.range(i);
            end
            Vnew.(param_all{jsig}) = x;
        else
            warning('No correction for signals_raw');
        end
    end
end

