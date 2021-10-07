%% Function to export data from ENSITE file.
% Only suitable for St. Jude Medical. File Revision : Velocity 2.0
%% MO 04/2014
function [Data]=load_data_from_ensite_mo(filename1)

% [Data]=load_data_from_ensite_mo(filename1)
% INPUT     
%   - filenemame1: string with path+name+format
% OUTPUT 
%   - Data: strunct array with fields:
%       Type : type of data (ECG, EGMs, Virtual)
%       wavenames : name of waves
%       traces : data (as many columns as length(wavenames))
%       Times : matrices with time info
%       Time_name : type of time variable
%%

hb = waitbar(0,'Importing Data ...');
fid=fopen(filename1);

%% open file. check if this is ENSITE file. get format type
% i.e. find "St. Jude Medical Virtual  data export; file format revision 3"
% in first row of the file
a=fgetl(fid);
%keyboard
if (strfind(a, 'Virtual  data export'))
    format = 0; % old format (Rev 3)
    
elseif (strfind(a, 'Waveform  data export'))
    format = 1;  %new format (Rev 4)
    
elseif (strfind(a, 'Velocity_1.0'))
    format = 2; % Velocity export
    
elseif (strfind(a, 'St. Jude Medical. File Revision : Velocity 2.0'))
    
    format = 2; %Velocity export
    disp('this is velocity 2.0 file')
    
elseif (strfind(a, 'St. Jude Medical. File Revision : Velocity 3.0'))
    
    format = 4; %Velocity export
    disp('this is velocity 3.0 file')
    
elseif (strfind(a, 'MC_DataTool ASCII conversion'))
    
    format = 5; %MC rack export
    disp('this is a MC rack file')
    
elseif (strfind(a, '[Header]'))
    if strfind(fgetl(fid), 'File Type: 1') && strfind(fgetl(fid), 'Version: 1')
        format = 3; % Bard Format
    else
        disp('Unknown/unsuported Bard format');
        return;
    end
else
    disp('Unknown/unsuported ENSITE/Bard format');
    return;
end

%%
if format == 2 %Velocity format
    seekstr = 't_dws,t_secs,t_usecs,t_ref';  %channel names header line - next line start of the actual data
    a=fgetl(fid);
    Data.Type = a(findstr(a,':')+2:end);
    while 1
        a=fgetl(fid);  
        if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
        if (strfind(a, seekstr)), break, end
    end
    
    %% 
    a(end+1) = ',';
    wavenames_temp = line2cellArray(a); % Name of waveform
    dataposition = ftell(fid); 
    tempdata = line2cellArray(fgetl(fid)); % read the first line of data to get the number of columns
    
%     fseek(fid,0,'bof');
%     seekstr = 'Export Start Time (secs usecs)';
%     while 1
%         a=fgetl(fid);  % "fgetLINE" not "fgetONE"
%         if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
%         if (strfind(a, seekstr)), break, end
%     end
%     %get starttime
%     jj = findstr(a,':')+1;
%     starttime = str2num(a(jj:end));
%     starttime = starttime(1)*1000 + starttime(2)/1000
%     %get endtime
%     fseek(fid,0,'bof');
%     %     seekstr = 'Export End Time (secs usecs)';
%     seekstr = 'Export End Time(secs usecs)';  % mo
%     
%     while 1
%         a=fgetl(fid);  % "fgetLINE" not "fgetONE"
%         if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
%         if (strfind(a, seekstr)), break, end
%     end
%     jj = findstr(a,':')+1;
%     endtime = str2num(a(31:end)); endtime = endtime(1)*1000 + endtime(2)/1000;
%     %get endtime
%     timepoints = (endtime - starttime)*2; %in secs - number of samples
%     
    %go back to the data position (end of the header
    fseek(fid,dataposition,'bof'); 
    %makeformat string, put %s as number of data channels
    formatstr = '%*s '; % skip the first one
    ind=0;
    Data.wavenames = cell(1);
    for kk=2:length(tempdata)
        if ~isempty(regexp(wavenames_temp{kk},'_ds'))|~isempty(regexp(wavenames_temp{kk},'_ps'))
            formatstr = [formatstr,'%*f '];
        else
            ind=ind+1;
            Data.wavenames{ind} = wavenames_temp{kk};
            formatstr = [formatstr,'%f '];
        end
    end
    waitbar(1/3,hb)
    %read all data into the one large cellarray
    celltraces = textscan(fid, formatstr, 'delimiter', ',','collectoutput',1);
    %     wavenames = wavenames(1:size(celltraces,2));
    %     ar(1/3,hb);
    %     iiok = cellfun(@isnumeric,celltraces);
    %     traces = cell2mat(celltraces(iiok));
    waitbar(2/3,hb)
    %     wavenames = wavenames(2:length(celltraces));
%     frequency = 2000;%round(1/(traces(2,1) - traces(1,1)));
    
    %%
    Data.Times = celltraces{1}(:,strncmp(Data.wavenames,'t_',2));%clear celltraces
    Data.Times_name = Data.wavenames(strncmp(Data.wavenames,'t_',2));%clear celltraces
    Data.traces = celltraces{1}(:,~strncmp(Data.wavenames,'t_',2));%clear celltraces
    Data.wavenames(strncmp(Data.wavenames,'t_',2))=[];%clear celltraces
    Data.traces(end,:)=[]; % remove nan due to 'eof'
    Data.Times(end,:)=[]; % remove nan due to 'eof'
    %     wavename2 = wavenames(~strncmp(wavenames,'t_',2));
    %     traces2 = traces2(:,~ [~cellfun(@isempty,regexp(wavename2,'_ds')) | ~cellfun(@isempty,regexp(wavename2,'_ps'))]);
    %     wavename2 = wavename2(~ [~cellfun(@isempty,regexp(wavename2,'_ds')) | ~cellfun(@isempty,regexp(wavename2,'_ps'))]);
    
    %     Data.traces = traces2;
    %     Data.wavename = wavename2;
    waitbar(1,hb)
    
    if isequal(Data.Type,'EPcathBIO_RAW')
        fseek(fid,0,'bof');
        seekstr1 = 'used channels';  % mo
        
        while 1
            a=fgetl(fid);  % "fgetLINE" not "fgetONE"
            if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
            if (strfind(a, seekstr1)), break, end
        end
        ind = 0;
        while ~isequal(a,'')&ind<100
            ind = ind+1;
            a = fgetl(fid);
            Data.Info_channels{ind} = a;
        end
    end
    
else
    error('Format not recognized, ask Michele')
    
    
end


%% Close the file
fclose(fid);

%now a function to help clean up any NaN's in the file - as with velocity
%with different (?) export function - for Justine;

disp('Check import')
close(hb)


%% some help functions
function output=line2cellArray(line) %create

nchannels=strfind(line,',');
previndex = 1;

for kk=1:length(nchannels)
    output{kk} = line(previndex:nchannels(kk)-1);
    previndex = nchannels(kk)+1;
end

function output=line2array(line)

nchannels=strfind(line,',');
previndex = 1;

for kk=1:length(nchannels)
    if kk == 1
        output(kk) = convert2seconds(line(previndex:nchannels(kk)-1));
    else
        output(kk) = str2double(line(previndex:nchannels(kk)-1));
    end
    previndex = nchannels(kk)+1;
end

function finaloutput=convert2seconds(tvector)

% Example
% With "jump"                   Without jump
% 22.42.05.655,4.67923
% 22.42.05.656,4.67972  |    22.42.05.656,4.67972
% 22.42.05.661,0        |    22.42.05.656,4.68021
% 22.42.05.662,0.000491 |    22.42.05.657,4.6807

finaloutput = [];
for i=1:length(tvector)
    text1 = tvector{i};
    ndots = strfind(text1,'.');
    output = str2double(text1(1:ndots(1)-1))*60*60;
    output = output + str2double(text1(ndots(1)+1:ndots(2)-1))*60;
    output = output + str2double(text1(ndots(2)+1:end));
    finaloutput = [finaloutput;output];
end





%% this is example of file format Velocity_1 (a: Virtuals)

% St. Jude Medical,Virtual, data export; file format revision,Velocity_1.0
% t_dws: timepoint as seen in DWS in HH.MM.SS.mSmSmS format
% t_ref: timepoint relative to first timepoint in seconds
% V-En: enGuide virtual@active electrode
% V#-##: virtual bank# channel ##
%
% t_dws,t_ref,V-En,V1-1,V1-2,V1-3,V1-4,V1-5,V1-6,V1-7,V1-8,V1-9,V1-10,V1-11,V1-12,V1-13,V1-14,V1-15,V1-16,V2-1,V2-2,V2-3,V2-4,V2-5,V2-6,V2-7,V2-8,V2-9,V2-10,V2-11,V2-12,V2-13,V2-14,V2-15,V2-16,
% 22.42.00.976,0,-1.29441,-0.046046,0.197334,0.238233,-0.0164239,0,0,0,0,0,0,0,0,0,0,0,0,-0.226746,-0.0657985,0.0170073,-3.82667,0,0,0,0,0,0,0,0,0,0,0,0,
% 22.42.00.977,0.000491,-1.30211,-0.044075,0.208643,0.241844,-0.01533,0,0,0,0,0,0,0,0,0,0,0,0,-0.22511,-0.063012,0.0169207,-3.83384,0,0,0,0,0,0,0,0,0,0,0,0,
% ...
% 22.43.01.976,0.47627,-0.330662,-0.160564,-0.0520099,-0.0206661,0.0283873,0,0,0,0,0,0,0,0,0,0,0,0,-0.213691,-0.190305,-0.25384,-0.114204,0,0,0,0,0,0,0,0,0,0,0,0,
% Number of Timepoints: 124106
% Time Interval       : 22.43.01.500-22.43.01.976
% Time Differential   : 18.00.00.000
% EOF

%% this is example of file format Velocity_1 (b: ECG)
% St. Jude Medical,ECGWaveform, data export; file format revision,Velocity_1.0
% t_dws: timepoint as seen in DWS in HH.MM.SS.mSmSmS format
% t_ref: timepoint relative to first timepoint in seconds
%
% t_dws,t_ref,I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6,
% 22.42.00.976,0,0.128,0.234,0.107,-0.181,0.011,0.17,0.062,0.38,0.491,0.709,0.541,0.401,
% ...
% 22.43.02.840,1.07971,0.019,0.022,0.002,-0.021,0.009,0.012,0.008,0.267,0.37,0.201,0.214,0.126,
% Number of Timepoints: 125866
% Time Interval       : 22.43.01.761-22.43.02.841
% Time Differential   : 18.00.00.000
% EOF

%% this is example of file format Velocity_1 (c: ep cath bio)
% St. Jude Medical,EPcathBIOWaveform, data export; file format revision,Velocity_1.0
% Number of Catheters = 2
% Catheter[0](name, num electrodes): RV,10
% 	Electrode[0](name, channel): D,20
% 	Electrode[1](name, channel): 2,21
% 	Electrode[2](name, channel): 3,22
% 	Electrode[3](name, channel): 4,23
% 	Electrode[4](name, channel): 5,24
% 	Electrode[5](name, channel): 6,25
% 	Electrode[6](name, channel): 7,26
% 	Electrode[7](name, channel): 8,27
% 	Electrode[8](name, channel): 9,28
% 	Electrode[9](name, channel): 10,29
% Catheter[1](name, num electrodes): CS,8
% 	Electrode[0](name, channel): D,0
% 	Electrode[1](name, channel): 2,1
% 	Electrode[2](name, channel): 3,2
% 	Electrode[3](name, channel): 4,3
% 	Electrode[4](name, channel): 5,4
% 	Electrode[5](name, channel): 6,5
% 	Electrode[6](name, channel): 7,6
% 	Electrode[7](name, channel): 8,7
%
% used channels:
% channel,catheter name,electrode name
% 0,CS,D
% 1,CS,2
% 2,CS,3
% 3,CS,4
% 4,CS,5
% 5,CS,6
% 6,CS,7
% 7,CS,8
% 20,RV,D
% 21,RV,2
% 22,RV,3
% 23,RV,4
% 24,RV,5
% 25,RV,6
% 26,RV,7
% 27,RV,8
% 28,RV,9
% 29,RV,10
%
% t_dws: timepoint as seen in DWS in HH.MM.SS.mSmSmS format
% t_ref: timepoint relative to first timepoint in seconds
% c###: channel ###
%
% t_dws,t_ref,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39,c40,c41,c42,c43,c44,c45,c46,c47,c48,c49,c50,c51,c52,c53,c54,c55,c56,c57,c58,c59,c60,c61,c62,c63,c64,c65,c66,c67,c68,c69,c70,c71,c72,c73,c74,c75,c76,c77,c78,c79,c80,c81,c82,c83,c84,c85,c86,c87,c88,c89,c90,c91,c92,c93,c94,c95,c96,c97,c98,c99,c100,c101,c102,c103,c104,c105,c106,c107,c108,c109,c110,c111,c112,c113,c114,c115,c116,c117,c118,c119,c120,c121,c122,c123,c124,c125,c126,c127,c128,c129,c130,c131,
% 22.56.27.546,0,-13.41,-14.047,-14.586,-15.562,-14.741,-14.921,-14.541,-14.747,-16.638,-15.858,-15.062,-14.769,-15.539,-16.674,-50,-50,-16.525,-14.429,-16.93,-16.897,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-13.135,-13.575,-12.854,-12.032,-13.145,-13.422,-12.262,-13.405,-13.07,-13.86,-12.31,-12.118,-12.097,-13.691,-13.69,-13.31,-13.17,-12.485,-12.698,-11.842,-15.4,-11.573,-13.592,-11.569,-12.117,-12.66,-12.425,-11.137,-13.052,-14.683,-14.002,-13.512,-15.685,-14.216,-14.166,-10.621,-12.799,-10.762,
% ...
% 22.56.45.543,2.29297,6.881,6.066,5.889,4.978,5.225,5.235,5.627,5.343,6.546,5.541,5.786,6.116,5.647,4.681,0.062,0.02,4.81,6.375,4.224,5.026,0.949,-0.098,-50,-50,-18.551,-14.843,-15.38,-15.106,-21.688,-26.13,4.82,3.915,4.87,5.223,4.143,3.97,5.014,4.058,4.233,3.621,5.153,5.578,5.365,3.944,3.831,3.993,4.204,4.844,4.513,5.1,1.958,5.624,3.725,5.747,5.236,4.705,-3.005,-2.139,-5.242,-5.776,-6.398,-4.542,-8.364,-6.707,2.781,6.344,4.353,6.434,
% Number of Timepoints: 36622
% Time Interval       : 22.56.43.250-22.56.45.546
% Time Differential   : 18.00.00.000
% EOF

%% this is example of file format revision 4

% St. Jude Medical Waveform  data export; file format revision 4
% Exported from study study_esi1179_2006:08:07:09:42:53
% Beginning sample number (time): 11:20:21:0000
% Ending sample number (time):    11:21:13:1199
%
% Number of waves (columns):      14
% Number of samples (rows):       63600
%
% *********************************************************************
% * Map Labels and ( x y z ) coordinates, surface name and surface ID *
% *********************************************************************
% 3DP APEX: ( x  y  z ) = (    -7.23   -16.16   -32.79 ) LV 0
% 3DP ROOT: ( x  y  z ) = (     9.40     7.55    46.07 ) LV 0
%     EA: ( x  y  z ) = (   -18.98    -8.30   -37.95 ) LV 0
%     BO: ( x  y  z ) = (     7.68     1.26   -34.30 ) LV 0
%
% *********************************************************************
% * Map Lesion and ( x y z ) coordinates, surface name and surface ID *
% *********************************************************************
%
% ****************************************
% * Wave Names and ( x y z ) coordinates *
% ****************************************
% Wave 0 = ECG I
% Wave 1 = ECG aVF
% Wave 2 = ECG V1
% Wave 3 = ABL D-2
% Wave 4 = ABL 3-4
% Wave 5 = ECG V2
% Wave 6 = ECG V6
% Wave 7 = Virtual: ( x  y  z ) = (     2.87    18.14     0.00 )
% Wave 8 = Virtual: ( x  y  z ) = (     4.84    20.54   -14.43 )
% Wave 9 = Virtual: ( x  y  z ) = (     5.81    21.92   -21.31 )
% Wave 10 = Virtual: ( x  y  z ) = (     5.09    20.58   -27.74 )
% Wave 11 = Virtual: ( x  y  z ) = (     2.95    16.87   -30.27 )
% Wave 12 = CS D-2
% Wave 13 = CS 3-4
%
% Begin data
%  Wave  0  Wave  1  Wave  2  Wave  3  Wave  4  Wave  5  Wave  6  Wave  7  Wave  8  Wave  9  Wave 10  Wave 11  Wave 12  Wave 13
%   -0.093    0.753    0.091    0.608    0.010    0.243    0.361    0.038    0.326    0.318    0.309    0.304    0.216   -0.215
%   -0.090    0.740    0.092    0.541    0.055    0.246    0.358   -0.047    0.098    0.081    0.062    0.052    0.241   -0.236
%   -0.087    0.723    0.092    0.296    0.097    0.247    0.352   -0.173   -0.268   -0.301   -0.338   -0.359    0.284   -0.260
%   ... lots of rows chopped out here
%   -0.085    0.700    0.090   -0.038    0.100    0.245    0.340   -0.271   -0.570   -0.618   -0.671   -0.700    0.305   -0.256
%   -0.083    0.666    0.087   -0.411    0.080    0.239    0.323   -0.348   -0.825   -0.885   -0.951   -0.989    0.310   -0.231
%EOF




%% this is example of file format revision 3

% St. Jude Medical Virtual  data export; file format revision 3
% Exported from study study_esi1264_2008:06:20:13:46:49
% Beginning sample number (time): 14:28:34:0822
% Ending sample number (time):    14:28:36:0103
%
% Number of signals (columns):    64
% Number of samples (rows):       1682
%
% **********************************************************
% * Map Labels and ( x y z ) coordinates and geometry name *
% **********************************************************
%     apex: ( x  y  z ) = (   -35.25   -25.61   -66.04 ) RV 0
%     Anterior: ( x  y  z ) = (    -1.35   -28.32     7.06 ) RV 0
%     Septal: ( x  y  z ) = (   -49.88   -16.86    14.13 ) RV 0
%     Outflow: ( x  y  z ) = (     5.00    23.80    15.77 ) RV 0
%
% **********************************************************
% * Map Lesion and ( x y z ) coordinates and geometry name *
% **********************************************************
%
% Begin data
%    7.33    3.33   -3.42   -8.19   -8.25   -3.30    2.99    7.00   21.28    8.94   -8.04  -23.06  -27.61   -9.81    7.62   16.46   23.58   10.38   -7.76  -28.79  -38.87  -11.71    9.18   19.13   20.24    9.50   -7.56  -26.14  -51.06  -12.63    7.15   14.75   17.01    8.57   -7.31  -23.42  -52.98  -13.73    6.73   11.42   12.03    6.82   -6.84  -18.96  -47.27  -11.15    5.30    9.51    8.38    4.16   -4.71  -14.00  -38.42   -8.22    4.24    7.80    3.72    1.62   -1.89   -5.90   -6.88   -2.46    1.82    3.80
%    3.04    8.04    8.25    3.39   -3.42   -7.96   -7.23   -2.90    8.81   21.59   19.41    9.55  -11.44  -23.68  -18.41   -6.82    9.77   25.05   18.73   11.93  -16.10  -28.27  -22.16   -7.92    8.39   22.94   18.24   10.83  -21.15  -30.50  -17.26   -6.11    7.05   20.69   17.64    9.70  -21.94  -33.15  -16.25   -4.73    4.98   16.46   16.52    7.86  -19.58  -26.91  -12.79   -3.94    3.47   10.05   11.36    5.80  -15.91  -19.83  -10.23   -3.23    1.54    3.91    4.56    2.44   -2.85   -5.94   -4.39   -1.57
%   39.90   43.76   44.90   44.59   44.89   43.34   39.32   38.08   34.47   34.98   31.44   37.35   44.73   38.36   29.82   26.66   17.06   18.12   13.55   20.82   28.11   20.45   16.03   13.83    4.36    4.94    3.93    5.63   10.99    6.57    3.72    3.17   -3.66   -4.45   -3.80   -5.04  -11.41   -7.14   -3.50   -2.46   -8.70  -11.91  -11.94  -13.72  -34.19  -19.47   -9.25   -6.88  -13.58  -16.29  -18.40  -22.68  -62.24  -32.13  -16.57  -12.64  -20.25  -21.30  -24.82  -32.10  -37.46  -32.32  -23.92  -20.67
%
%   0.005   0.014  -0.014  -0.062  -0.095  -0.093  -0.074  -0.037   0.037   0.078  -0.012  -0.134  -0.198  -0.191  -0.144  -0.071  -0.012   0.062  -0.102  -0.287  -0.357  -0.338  -0.246  -0.154  -0.136  -0.062  -0.300  -0.529  -0.579  -0.544  -0.385  -0.265  -0.347  -0.380  -0.604  -0.761  -0.780  -0.739  -0.578  -0.428  -0.562  -0.823  -0.900  -0.834  -0.844  -0.814  -0.663  -0.499  -0.709  -0.918  -0.855  -0.741  -0.778  -0.706  -0.196  -0.173  -0.438  -0.497  -0.519  -0.504  -0.467  -0.384  -0.140  -0.178
%   0.009   0.018  -0.009  -0.056  -0.088  -0.087  -0.068  -0.032   0.040   0.082  -0.005  -0.126  -0.191  -0.183  -0.136  -0.065  -0.007   0.066  -0.096  -0.280  -0.349  -0.330  -0.232  -0.140  -0.127  -0.062  -0.299  -0.524  -0.570  -0.536  -0.368  -0.242  -0.333  -0.385  -0.605  -0.751  -0.768  -0.729  -0.569  -0.410  -0.544  -0.816  -0.886  -0.814  -0.828  -0.797  -0.665  -0.496  -0.683  -0.880  -0.816  -0.710  -0.752  -0.679  -0.200  -0.181  -0.410  -0.456  -0.477  -0.466  -0.435  -0.355  -0.121  -0.166
%   0.014   0.024  -0.002  -0.049  -0.081  -0.080  -0.062  -0.026   0.045   0.088   0.003  -0.118  -0.183  -0.175  -0.127  -0.057   0.001   0.073  -0.089  -0.274  -0.342  -0.323  -0.220  -0.125  -0.114  -0.059  -0.299  -0.520  -0.562  -0.529  -0.355  -0.221  -0.318  -0.388  -0.605  -0.742  -0.758  -0.720  -0.563  -0.395  -0.526  -0.803  -0.869  -0.796  -0.812  -0.778  -0.664  -0.493  -0.660  -0.840  -0.778  -0.682  -0.726  -0.650  -0.200  -0.192  -0.386  -0.418  -0.438  -0.430  -0.404  -0.326  -0.104  -0.158
%   0.020   0.029   0.004  -0.043  -0.076  -0.074  -0.056  -0.020   0.051   0.093   0.009  -0.113  -0.178  -0.169  -0.119  -0.049   0.009   0.077  -0.086  -0.271  -0.337  -0.317  -0.210  -0.113  -0.103  -0.060  -0.301  -0.519  -0.556  -0.524  -0.348  -0.206  -0.307  -0.391  -0.606  -0.735  -0.751  -0.712  -0.561  -0.388  -0.511  -0.787  -0.854  -0.784  -0.799  -0.761  -0.660  -0.492  -0.643  -0.803  -0.747  -0.662  -0.704  -0.624  -0.199  -0.205  -0.371  -0.389  -0.410  -0.403  -0.380  -0.304  -0.092  -0.157
%   ... lots of rows chopped out here
%   0.021   0.031   0.005  -0.042  -0.074  -0.072  -0.053  -0.018   0.053   0.095   0.010  -0.112  -0.177  -0.167  -0.114  -0.044   0.015   0.078  -0.088  -0.274  -0.337  -0.316  -0.207  -0.105  -0.096  -0.063  -0.308  -0.523  -0.556  -0.524  -0.350  -0.201  -0.301  -0.394  -0.610  -0.734  -0.750  -0.708  -0.565  -0.390  -0.501  -0.770  -0.842  -0.778  -0.791  -0.748  -0.656  -0.494  -0.636  -0.772  -0.726  -0.652  -0.690  -0.603  -0.198  -0.222  -0.368  -0.372  -0.394  -0.387  -0.365  -0.290  -0.087  -0.165
%   0.023   0.032   0.006  -0.041  -0.072  -0.069  -0.050  -0.015   0.056   0.096   0.010  -0.113  -0.176  -0.164  -0.109  -0.038   0.019   0.077  -0.092  -0.277  -0.337  -0.316  -0.205  -0.100  -0.092  -0.068  -0.316  -0.526  -0.555  -0.524  -0.356  -0.201  -0.297  -0.395  -0.612  -0.733  -0.749  -0.703  -0.570  -0.397  -0.494  -0.749  -0.830  -0.775  -0.784  -0.736  -0.647  -0.496  -0.633  -0.746  -0.712  -0.649  -0.679  -0.586  -0.194  -0.237  -0.372  -0.364  -0.387  -0.379  -0.357  -0.282  -0.086  -0.177
%   0.029   0.036   0.010  -0.037  -0.068  -0.064  -0.043  -0.008   0.063   0.099   0.011  -0.111  -0.173  -0.159  -0.102  -0.030   0.026   0.079  -0.094  -0.278  -0.335  -0.313  -0.203  -0.094  -0.084  -0.068  -0.320  -0.527  -0.553  -0.521  -0.363  -0.203  -0.291  -0.388  -0.609  -0.730  -0.747  -0.697  -0.573  -0.406  -0.486  -0.721  -0.815  -0.773  -0.777  -0.722  -0.631  -0.493  -0.633  -0.722  -0.703  -0.649  -0.669  -0.570  -0.182  -0.247  -0.380  -0.361  -0.387  -0.377  -0.351  -0.276  -0.085  -0.189
%EOF

%% this is example of file format new velocity

%St. Jude Medical. File Revision : Velocity 2.0
%Export Data Element : Virtuals
%Export from Study : study_dwsG020599_2009_05_07_07_38_43
%%Export from Segment : 400 RVOT EDR
%Export User Comments : Anterior and Posterior Virtuals
%Export files Stored in Dir : /var/STJ/Clinical/systemStudy/study_dwsG020599_2009_05_07_07_38_43/export/Seg09_2010_07_08_14_27_52/
%Export Start Time (h.m.s): 03.45.10.753
%Export Start Time (secs usecs) : 1241685910 753038
%Export End Time(h.m.s) : 03.51.59.217
%Export End Time(secs usecs) : 1241686319 217038
%Export Duration(h.m.s) : 00.06.48.464
%
%t_dws: timepoint as seen in DWS in HH.MM.SS.mSmSmS format
%t_secs: DWS timepoint secs part
%t_usecs: DWS timepoint usecs part
%t_ref: timepoint relative to first timepoint in seconds
%V-En: enGuide virtual@active electrode
%V#-##: virtual bank# channel ##
%
%t_dws,t_secs,t_usecs,t_ref,V-En,V1-1,V1-2,V1-3,V1-4,V1-5,V1-6,V1-7,V1-8,V1-9,V1-10,V1-11,V1-12,V1-13,V1-14,V1-15,V1-16,V2-1,V2-2,V2-3,V2-4,V2-5,V2-6,V2-7,V2-8,V2-9,V2-10,V2-11,V2-12,V2-13,V2-14,V2-15,V2-16,
%03.45.10.753,1241685910,753032,0,0,-0.199679,-0.401143,-0.494229,-0.522749,0,0,0,0,0,0,0,0,0,0,0,0,-0.938491,-1.36356,-1.90761,-1.41113,0,0,0,0,0,0,0,0,0,0,0,0,0,
%03.45.10.753,1241685910,753523,0.000491,0,-0.201898,-0.40905,-0.507493,-0.53033,0,0,0,0,0,0,0,0,0,0,0,0,-0.941976,-1.36395,-1.90437,-1.40908,0,0,0,0,0,0,0,0,0,0,0,0,0,
%... lots of rows....
%03.51.59.216,1241686319,216551,408.032,0,-0.046215,-0.0508034,-0.0814276,-0.0643903,0,0,0,0,0,0,0,0,0,0,0,0,-0.0802791,-0.0917099,-0.122675,-0.128662,0,0,0,0,0,0,0,0,0,0,0,0,-0.133984,
%EOF