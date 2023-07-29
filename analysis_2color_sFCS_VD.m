

close all
fclose all;
clear reload
if exist('linetime', 'var')
 answer = questdlg('Do you want to re-use the parameters from the last image?', ...
                          '', ...
                          'Yes','No','No');
        % Handle response
        switch answer
         case 'No'  
            reload = 1;
         case 'Yes'
            reload = 0;  
        end

        
end

if exist('reload', 'var')
if reload==1
    clearvars -except reload
    clear global
end
else
    clear
    clear global
end

if exist('reload', 'var')
if reload==0
    clearvars -except reload czifilename czi_path metafile_czi2
    clear global
end

end


  

%%%%%% Define some parameters
SP=5;                     % <--------- structural parameter taken from the AF488 calibration
%SP=5.41;                     % <--------- structural parameter taken from the AF488 calibration

groupID='hIntegegrin-HApeptide';             % <--------- used for the post-hoc analysis (e.g. monomer, dimer, M1, etc.)
% lsp_Ch1=1;
% lsp_Ch2=0.5;
%%%%%% Define input & output path
input='H:\Anne\2022-10-01_Laurdan_NR12S_Di4ANEP\sFCS\Export\hIntegegrin-HApeptide';       
output='c:/';
czi_input='H:\Anne\2022-10-01_Laurdan_NR12S_Di4ANEP\sFCS\RAW\hIntegegrin-HApeptide';

sigbinning=1;  %a value of 2 makes you lose the first point of the AC/CC but decreases the noise considerably



%declare all necessary functions (located under ..\MATLAB\scanningFCS\sFCS functions)as global
global movingaveragewindowsize curveincl1_poly curveincl2_poly corfit1 corfit2 corfit1_poly corfit2_poly  ItraceIICh1_poly ItraceIICh2_poly blocksize spatialfilter depletioncorrection...
   polygrad intensityfilter intensitythreshold intensityfilterwindow membranewidth S waist1 waist Global_Metafile factorgreen factorred  correctionfitseg1aa correctionfitseg2aa

PIE=1;                  % 1 PIE mode, 0 normal mode                   

loadandpool=0;          % 0 load new raw file, 1 load previously computed ACFs

calculateACFonly=0;     % Calculate ACF only and return afterwards
ACFselection=1;         % 1 divide timetrace in ACFsegmentnumb segments and select or discard segment ACFs
intensityselection=1;
intensitycarpetselection=1;
greenredchannel=0;
backgroundcorrectionCh1=1; %1: background will be subtracted, 0: no subtraction
backgroundcorrectionCh2=1; %1: background will be subtracted, 0: no subtraction
polygrad=6; % polygrad bleaching and slow-fluctuation corrections, 8 is just an arbitrary number & does not need to be changed for the moment. You can/should test what happens though
Jonas=1;% polygrad bleaching correction from Jonas Ries, is important only if someone uses the "experimental" polynomial correction


if ACFselection==1
    global x0 y0 fixed  fitting_poly_GUI1  fitting_poly_GUI2 curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh1_poly correlationcurvesCh2_poly correlationcurvesChCC_poly  sigmas1curves_poly sigmas2curves_poly sigmascccurves_poly correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves corfit ItraceIICh1 ItraceIICh2 goon
end
if intensityselection==1
    global segmentsincl1 segmentsincl2 Ifull1 Ifull2 timeline1binned timeline2binned numberofsegments segmentlength correctionfitseg1 correctionfitseg2 ff1 ff2
end
if intensitycarpetselection==1
    global gooncarpets carpetsegments carpetCh1 carpetCh2
end

carpetsegments=50;
%ACFsegmentnumb=16;      % <----------------------------------------
                         % Input here!!!, everything between 10 and 20 is ok, but it should be even (e.g. 400,000 frames -> 4 per 100 000 scans)
%ACFsegmentnumb=20;       % 250,000 frames (8 per 100 000 scans) 
ACFsegmentnumb=10;      % 250,000 frames (4 per 100 000 scans)
numberofsegments=10;
loadsimulation=0;
kymoload=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
%if loadandpool==1              
%    S=6.84;
%    scantime=945.45e-06;
    %scantime=473e-06;
    %globalfit=1;
%     pixelsize=0.3; %in um
%     d=pixelsize;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correction of instabile measurements
fouriercorrection=0;
fourierwindow_env=20;
loadfourierenvelope=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%load a new raw file
if loadandpool==0
    if loadsimulation==0
        % Initialize framedata based on image aquisition parameters
        if kymoload==1
            path= uigetdir; 
            files=dir([path '\*.tif']);
            fra
            for i=1:size(files,1)
                [namedata,remain]=strtok(files(i).name,'.');
                inputfilename=files(i).name;
                framedatai=double(imread([path '/' inputfilename]));
                framedata(:,:,i)=framedatai;
            end
            line1array=framedata(:,:,1);
            line2array=framedata(:,:,2);
        else
            % Import of .tif files
            %[inputfilename, path]=uigetfile;
            [inputfilename, path]=uigetfile([input '\*.tif']);   % <----------------------------------------
            linedata=imread([path inputfilename]);
            linedata=blocksum(linedata, sigbinning);
            
            %Read in imagedatafile -> checkpoint for bitdepth
            imagemetadata=imfinfo([path inputfilename]);
            
            % checkpoint image criteria
%             bitdepth=imagemetadata.BitDepth/3;
%             pixelnumbX=imagemetadata.Width;
%             if bitdepth~=8
%                 error('image bitdepth is not 8-bit')
%             end
%             
%             if pixelnumbX~=256
%                 error('pixelnumber is not 256')
%             end
            
            % Import of *czi files
            if exist('reload', 'var')
            if reload==1
            [czifilename,czi_path]=uigetfile([czi_input '\*.czi']);
metafile_czi2=bfopen_justinfo([czi_path czifilename]);
            end
            else
                [czifilename,czi_path]=uigetfile([czi_input '\*.czi']);
metafile_czi2=bfopen_justinfo([czi_path czifilename]);
            end
            
imagemetadata=metafile_czi2{1,2}; 

allKeys = arrayfun(@char, imagemetadata.keySet.toArray, 'UniformOutput', false);
for i=1:size(allKeys,1)
allKeys{i,2}=imagemetadata.get(allKeys{i,1});
end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if greenredchannel==1
            line1array=double(linedata(:,:,2))+double(linedata(:,:,3)); % green channel ([r G b])
            else
            line1array=double(linedata(:,:,2));
            end
            line2array=double(linedata(:,:,1)); % red channel ([R g b])

            if length(line1array(:,1))==999999
                line1array=line1array(1:end-999,:);
                line2array=line2array(1:end-999,:);
            end
        end
    else
        % Code for loading simulations
        [inputfilename,path]=uigetfile;
        filepath=[path inputfilename];
        load(filepath);
        line1array=linearray(:,:,1);
        line2array=linearray(:,:,2);
    end


    movingaveragewindowsize=2000;
    blocksize=movingaveragewindowsize;
    if mod(length(line1array(:,1)),movingaveragewindowsize)>0 || mod(movingaveragewindowsize,2)==1;
        fprintf('Window Size not valid!\n')
        return
    end
    spatialfilter=2.5; % factor of stdev above and below mean membrane position
    depletioncorrection=1; % cell measurements always 1                 % <----------------------------------------
    intensityfilter=0;
    intensitythreshold=10^-5;
    intensityfilterwindow=100;
    binningwindow=100;
    binningwindow_GUI=50;
    %backgroundcorrectionCh1=0; %1: background will be subtracted, 0: no subtraction
    %backgroundcorrectionCh2=0; %1: background will be subtracted, 0: no subtraction

    %scantime_lines=1890e-06; %Time it takes to scan one frame (two lines)
    %scantime_lines=945.45e-06;         % <----------------------------------------
        scantime_lines=sigbinning*str2double(imagemetadata.get('Global Information|Image|Channel|LaserScanInfo|FrameTime #1'));
    %timebreak=5e-3; % Intervall in between subsequent scans (e.g. G-R_break_G-R_break...)
    timebreak=0;                       % <----------------------------------------
    scantime=scantime_lines+timebreak;
    %scantime=1890e-06+10e-3;
    %pixeltime=0.79e-06;             %two-color
        pixeltime=sigbinning*str2double(imagemetadata.get('Global Information|Image|Channel|LaserScanInfo|PixelTime #1'));
    %pixelsize=0.16; %in um
    %pixelsize=0.079;                     % <------ Zoom 10.6 correspond to a pixel size of 0.08 �m
    %pixelsize=imagemetadata.XResolution/1000;    % unit: Inch
        pixelsize=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'));  %unit m
    pixelsize=pixelsize*1e+06;              % unit �m
    %d=pixelsize;
    %pixels=length(linearray(1,:));
    %pixeltime=scantime/pixels;
    %linetime=pixeltime*length(linearray(1,:));
    
       laserpower=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|Laser|LaserPower #1'));
        LP_total=str2double(imagemetadata.get('Global Information|Instrument|LightSource|Power #1'));
        Attenuator_Transmission=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Attenuator|Transmission #1'));
        Attenuator_Detector=str2double(imagemetadata.get('Global Information|Image|Channel|Attenuation #1'));
        lsp_Ch1=((LP_total/laserpower)*Attenuator_Transmission)/10; %actual used laser power in %
    
    
     laserpower=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|Laser|LaserPower #2'));
        LP_total=str2double(imagemetadata.get('Global Information|Instrument|LightSource|Power #2'));
        Attenuator_Transmission=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Attenuator|Transmission #2'));
        Attenuator_Detector=str2double(imagemetadata.get('Global Information|Image|Channel|Attenuation #2'));
        lsp_Ch2=((LP_total/laserpower)*Attenuator_Transmission)/10; %actual used laser power in %
    
    
    
    
    if PIE==0
        linetime=scantime;
    else
       % linetime=scantime/2;
       linetime=0.5*scantime_lines;
    end
    S=SP;                             % <---------------------------------------- Structural parameter from the AF488 calibration
    x0=[100,0.05,S]; % Set in header or GUI manually later!
    y0=[100,0.05,S];
    fixed=[false false true];
    
    
    if depletioncorrection==1
    % Test of GUI for selection
    % Remove Signals close to the membrane (e.g. from vesicles)
    if intensitycarpetselection==1
        carpetCh1=line1array;
        carpetCh2=line2array;
        Intensitycarpets2ch
        gooncarpets=0;
        while gooncarpets==0
            pause(5)
        end
        line1array=carpetCh1;
        line2array=carpetCh2;
    end
    end
    
    % Figure: Fluorescence of all scanned lines
    scrsz=get(0,'ScreenSize');
    load('redmap.mat');
    load('greenmap.mat')
    
    blocksizedisp=1000;
    line1arraybinned=line1array;
    nbins=size(line1array,1)/blocksizedisp;
    for ii=1:nbins
        for kk=1:blocksizedisp
            line1arraybinned((ii-1)*blocksizedisp+kk,:)=mean(line1arraybinned((ii-1)*blocksizedisp+1:ii*blocksizedisp,:),1);
        end
    end
    
    line2arraybinned=line2array;
    nbins=size(line2array,1)/blocksizedisp;
    for ii=1:nbins
        for kk=1:blocksizedisp
            line2arraybinned((ii-1)*blocksizedisp+kk,:)=mean(line2arraybinned((ii-1)*blocksizedisp+1:ii*blocksizedisp,:),1);
        end
    end


    LineFl=figure('OuterPosition',[1 1 2*scrsz(3)/3 scrsz(4)],'Name','Line Fluorescence - Define ROI');
    subplot(2,2,3), plot(mean(line1array(1:size(line1array,1)/2,:),1),'b');
    hold on
    subplot(2,2,3), plot(mean(line1array(size(line1array,1)/2+1:end,:),1),'r');
    legend('1st half','2nd half')
    ax1=subplot(2,2,1);imagesc(line1arraybinned)
    colormap(ax1,greenmap)
    if isempty(factorgreen)
        factorgreen=1;
    end
    factorgreen=1;
    ax1.CLim(2)=ax1.CLim(2)/factorgreen
    xlabel('pixelposition (a.u.)')
    ylabel('time (linenumber)')
    mask1=roipoly; % Need two seperate ROIS? Probably yes for cell applications
%     line1arraymasked = mask.*line1array;
%     line2arraymasked = mask.*line2array;
    hold on
    % background correction
    if backgroundcorrectionCh1==1
        answer = questdlg('Where is the background with respect to the membrane?', ...
                          'Background selection', ...
                          'Left','Right','Left');
        % Handle response
        switch answer
         case 'Left'  
            backgleftCh1 = 1;
         case 'Right'
            backgleftCh1 = -1;  
        end
        %hold (handles.axes1,'on')
        mask_bgr1=roipoly;
        backgroundmasked1=mask_bgr1.*line1array;
        %hold (handles.axes1,'off')
        %linearraymasked=linearraymasked-0.5*mean(mean(backgroundmasked));
        % subtract from linearraymasked or linefluorescenceseries?
    end
    hold on
    subplot(2,2,4), plot(mean(line2array(1:size(line2array,1)/2,:),1),'b');
    hold on
    subplot(2,2,4), plot(mean(line2array(size(line2array,1)/2+1:end,:),1),'r');
    legend('1st half','2nd half')
    ax2=subplot(2,2,2);imagesc(line2arraybinned)
    colormap(ax2,redmap)
      if isempty(factorred)
        factorred=2;
      end
      factorred=2;
      ax2.CLim(2)=ax2.CLim(2)/factorred

    xlabel('pixelposition (a.u.)')
    ylabel('time (linenumber)')
    mask2=roipoly; % Need two seperate ROIS? Probably yes for cell applications
    hold on
    if backgroundcorrectionCh2==1
%         answer = questdlg('Where is the background with respect to the membrane?', ...
%                           'Background selection', ...
%                           'Left','Right','Left');
        % Handle response
        switch answer
         case 'Left'  
            backgleftCh2 = 1;
         case 'Right'
            backgleftCh2 = -1;  
        end
        %hold (handles.axes1,'on')
        mask_bgr2=roipoly;
        backgroundmasked2=mask_bgr2.*line2array;
        %hold (handles.axes1,'on')
        %linearraymasked=linearraymasked-0.5*mean(mean(backgroundmasked));
        % subtract from linearraymasked or linefluorescenceseries?
    end
    path2=uigetdir(output);        % <-----------------------------------------       
    savefig(LineFl,[path2 '\' inputfilename(1:end-4) ' LineFl ' int2str(i) '.fig'])
    
%     % background correction
%     if backgroundcorrectionCh1==1
%         hold (handles.axes1,'on')
%         mask_bgr1=roipoly;
%         backgroundmasked1=mask_bgr1.*line1array;
%         hold (handles.axes1,'off')
%         %linearraymasked=linearraymasked-0.5*mean(mean(backgroundmasked));
%         % subtract from linearraymasked or linefluorescenceseries?
%     end
%     if backgroundcorrectionCh2==1
%         hold (handles.axes1,'on')
%         mask_bgr2=roipoly;
%         backgroundmasked2=mask_bgr2.*line2array;
%         hold (handles.axes1,'on')
%         %linearraymasked=linearraymasked-0.5*mean(mean(backgroundmasked));
%         % subtract from linearraymasked or linefluorescenceseries?
%     end
    


    % Analysis

    % Polygonal Selection
    % mask=roipoly; % Need two seperate ROIS? Probably yes for cell applications
    line1arraymasked = mask1.*line1array;
    line2arraymasked = mask2.*line2array;
    
    line1arraymaskedbinned=line1arraymasked;
    nbins=size(line1array,1)/blocksizedisp;
    for ii=1:nbins
        for kk=1:blocksizedisp
            line1arraymaskedbinned((ii-1)*blocksizedisp+kk,:)=mean(line1arraymaskedbinned((ii-1)*blocksizedisp+1:ii*blocksizedisp,:),1);
        end
    end
    
    line2arraymaskedbinned=line2arraymasked;
    nbins=size(line2array,1)/blocksizedisp;
    for ii=1:nbins
        for kk=1:blocksizedisp
            line2arraymaskedbinned((ii-1)*blocksizedisp+kk,:)=mean(line2arraymaskedbinned((ii-1)*blocksizedisp+1:ii*blocksizedisp,:),1);
        end
    end
 
    figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Selected Lines')
    selected1=subplot(1,2,1);imagesc(line1arraymaskedbinned)
    colormap(selected1,greenmap)
    xlabel('pixelposition (a.u.)')
    ylabel('time (linenumber)')
    hold on 
    selected2=subplot(1,2,2);imagesc(line2arraymaskedbinned)
    colormap(selected2,redmap)
    xlabel('pixelposition (a.u.)')
    ylabel('time (linenumber)')

    % Alignment and Calculation of time series
    
    if backgroundcorrectionCh1==1
        [line1fluorescenceseries,line1arrayalignedba]=linealignfuncbgr(mask1,line1arraymasked,backgroundmasked1,backgleftCh1);
        membranewidth_Ch1=membranewidth;
        waist1_Ch1=waist1;
        waist_Ch1=waist1_Ch1*pixelsize;
    else
        [line1fluorescenceseries,line1arrayalignedba]=linealignfunc(line1arraymasked);
        membranewidth_Ch1=membranewidth;
        waist1_Ch1=waist1;
        waist_Ch1=waist1_Ch1*pixelsize;
    end
    
    if backgroundcorrectionCh2==1
        [line2fluorescenceseries,line2arrayalignedba]=linealignfuncbgr(mask2,line2arraymasked,backgroundmasked2,backgleftCh2);
        membranewidth_Ch2=membranewidth;
        waist1_Ch2=waist1;
        waist_Ch2=waist1_Ch2*pixelsize;
    else
        [line2fluorescenceseries,line2arrayalignedba]=linealignfunc(line2arraymasked);
        membranewidth_Ch2=membranewidth;
        waist1_Ch2=waist1;
        waist_Ch2=waist1_Ch2*pixelsize;
    end
    
    if (backgroundcorrectionCh1+backgroundcorrectionCh2)==2
        line1arrayalignedmasum=sum(line1arrayalignedba,1);
        line1arrayalignedmasum(isnan(line1arrayalignedmasum))=0;
        membfit1=fit([1:1:length(line1arrayalignedmasum)]',line1arrayalignedmasum'./sum(line1arrayalignedmasum),'gauss1');
        mu1=membfit1.b1;
        line2arrayalignedmasum=sum(line2arrayalignedba,1);
        line2arrayalignedmasum(isnan(line2arrayalignedmasum)==1)=0;
        membfit2=fit([1:1:length(line2arrayalignedmasum)]',line2arrayalignedmasum'./sum(line2arrayalignedmasum),'gauss1');
        mu2=membfit2.b1;
           % scrsz=   get(0,'ScreenSize');
        MemPeak=figure('OuterPosition',[scrsz(3)/3 scrsz(4)/4 2*scrsz(3)/3 scrsz(4)/4],'Name','Total Membrane fluorescence peak');
        plot(1:length(line1arrayalignedmasum),line1arrayalignedmasum./sum(line1arrayalignedmasum),'-g')
        hold on
        plot(1:length(line2arrayalignedmasum),line2arrayalignedmasum./sum(line2arrayalignedmasum),'-m')
        plot(membfit1,'--g')
        hold on
        plot(membfit2,'--m')
        xlabel('pixelposition')
        ylabel('intensity')
        legend('Line Ch1','Line Ch2','Fit Ch1', 'Fit Ch2')
        xlim([round(mu2-8) round(mu1+8)]);
        shiftfoci=abs(membfit1.b1-membfit2.b1)*pixelsize; %in um
        savefig(MemPeak,[path2 '\' inputfilename(1:end-4) ' _MemPeak ' int2str(i) '.fig'])
    end
    
    line1arrayalignedbabinned=line1arrayalignedba;
    nbins=size(line1array,1)/blocksizedisp;
    for ii=1:nbins
        for kk=1:blocksizedisp
            line1arrayalignedbabinned((ii-1)*blocksizedisp+kk,:)=mean(line1arrayalignedbabinned((ii-1)*blocksizedisp+1:ii*blocksizedisp,:),1);
        end
    end
    
    line2arrayalignedbabinned=line2arrayalignedba;
    nbins=size(line2array,1)/blocksizedisp;
    for ii=1:nbins
        for kk=1:blocksizedisp
            line2arrayalignedbabinned((ii-1)*blocksizedisp+kk,:)=mean(line2arrayalignedbabinned((ii-1)*blocksizedisp+1:ii*blocksizedisp,:),1);
        end
    end
    
    
    LineFl2=figure('OuterPosition',[scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Line Fluorescence (aligned)');
    aligned1=subplot(1,2,1);imagesc(line1arrayalignedbabinned);
    colormap(aligned1,greenmap)
    xlabel('pixelposition (a.u.)')
    ylabel('time (linenumber)') 
    hold on
    aligned2=subplot(1,2,2);imagesc(line2arrayalignedbabinned);
    colormap(aligned2,redmap)
    xlabel('pixelposition (a.u.)')
    ylabel('time (linenumber)') 
    savefig(LineFl2,[path2 '\' inputfilename(1:end-4) ' LineFl_aligned ' int2str(i) '.fig'])

    return
    
    % Figure: Fluorescence time series
%     figure('Name','Membrane Fluorescence time series')
    timeline=1:1:length(line1fluorescenceseries);
    if PIE==0
        timeline1=(timeline-1)'*scantime+scantime_lines;
        timeline2=(timeline-1)'*scantime+scantime_lines;
    else
        timeline1=(timeline-1)'*scantime+scantime_lines-linetime;
        timeline2=(timeline-1)'*scantime+scantime_lines;
    end
    f1=fit(timeline1,line1fluorescenceseries,'exp2');
    f2=fit(timeline2,line2fluorescenceseries,'exp2');
%     subplot(1,2,1),plot(f1,timeline1,line1fluorescenceseries);
%     xlabel('time (linenumber)')
%     ylabel('line fluorescence')
%     hold on
%     subplot(1,2,2),plot(f2,timeline2,line2fluorescenceseries);
%     xlabel('time (linenumber)')
%     ylabel('line fluorescence')
    bleachingfraction1=1-mean(line1fluorescenceseries(end-1000+1:end))/mean(line1fluorescenceseries(1:1000));
    bleachingfraction2=1-mean(line2fluorescenceseries(end-1000+1:end))/mean(line2fluorescenceseries(1:1000));
    
    % Figure: Binned fluorescence time series and Brightness
    lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
    line1fluorescenceseriesbinned=mean(reshape(line1fluorescenceseries(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries(1:end-lastbinlength))/binningwindow),1);
    line1fluorescenceseriesbinned=[line1fluorescenceseriesbinned mean(line1fluorescenceseries(end-lastbinlength+1:end))];
    line2fluorescenceseriesbinned=mean(reshape(line2fluorescenceseries(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries(1:end-lastbinlength))/binningwindow),1);
    line2fluorescenceseriesbinned=[line2fluorescenceseriesbinned mean(line2fluorescenceseries(end-lastbinlength+1:end))];
    correctionfit1=f1.a.*exp(f1.b*timeline1)+f1.c.*exp(f1.d*timeline1);
    correctionfit2=f2.a.*exp(f2.b*timeline2)+f2.c.*exp(f2.d*timeline2);
    %correctionfit=f.a1*sin(f.b1*timeline+f.c1)+f.a2*sin(f.b2*timeline+f.c2)+f.a3*sin(f.b3*timeline+f.c3)+f.a4*sin(f.b4*timeline+f.c4)+f.a5*sin(f.b5*timeline+f.c5)+f.a6*sin(f.b6*timeline+f.c6)+f.a7*sin(f.b7*timeline+f.c7)+f.a8*sin(f.b8*timeline+f.c8)+f.d1;
    
    correctionfit1binned=mean(reshape(correctionfit1(1:end-lastbinlength),binningwindow,length(correctionfit1(1:end-lastbinlength))/binningwindow),1);
    correctionfit1binned=[correctionfit1binned mean(correctionfit1(end-lastbinlength+1:end))];
    line1fluorescencebinnedresiduals=line1fluorescenceseriesbinned-correctionfit1binned;
    correctionfit2binned=mean(reshape(correctionfit2(1:end-lastbinlength),binningwindow,length(correctionfit2(1:end-lastbinlength))/binningwindow),1);
    correctionfit2binned=[correctionfit2binned mean(correctionfit2(end-lastbinlength+1:end))];
    line2fluorescencebinnedresiduals=line2fluorescenceseriesbinned-correctionfit2binned;
    
%     linebrightnessbinned=var(reshape(linefluorescenceseries(1:end-lastbinlength),binningwindow,length(linefluorescenceseries(1:end-lastbinlength))/binningwindow),1);
%     linebrightnessbinned=[linebrightnessbinned var(linefluorescenceseries(end-lastbinlength+1:end))];
%     meanbrightness=mean(linebrightnessbinned);
%     meanbrightnessfit=meanbrightness*ones(1,length(linebrightnessbinned));
%     linebrightnessresiduals=linebrightnessbinned-meanbrightness;
    hI=figure('OuterPosition',[scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Binned Membrane Fluorescence time series');
    timeline1_binned=1:1:length(line1fluorescenceseriesbinned);
    timeline2_binned=1:1:length(line2fluorescenceseriesbinned);
    if PIE==0
        timeline1_binned=(timeline1_binned-1)*binningwindow*scantime+scantime_lines*binningwindow;
        timeline2_binned=(timeline2_binned-1)*binningwindow*scantime+scantime_lines*binningwindow;
    else
        timeline1_binned=(timeline1_binned-1)*binningwindow*scantime+scantime_lines*binningwindow-binningwindow*linetime; %Correct?
        timeline2_binned=(timeline2_binned-1)*binningwindow*scantime+scantime_lines*binningwindow;
    end
    subplot(2,1,1),plot(timeline1_binned,line1fluorescenceseriesbinned)
    hold on
    subplot(2,1,1),plot(f1)
    xlabel('time')
    ylabel('line fluorescence Ch1')
    subplot(2,1,2),plot(timeline2_binned,line2fluorescenceseriesbinned)
    hold on
    subplot(2,1,2),plot (f2)
%     hold on
%     plot(timeline_binned,0*linefluorescencebinnedresiduals)
    xlabel('time')
    ylabel('line fluorescence Ch2')
    
    
    if intensityselection==1
        Ifull1=line1fluorescenceseriesbinned;
        Ifull2=line2fluorescenceseriesbinned;
        timeline1binned=timeline1_binned;
        timeline2binned=timeline2_binned;
        size(Ifull1);
        size(timeline1_binned);
        segmentlength=length(Ifull1)/numberofsegments;
        segmentsincl1=ones(1,numberofsegments);
        segmentsincl2=ones(1,numberofsegments);
        
        % GUI for selection of segment intensity traces
        Intselection2colors_new %Fehler in GUI (Klick in einem Kanal beeinflusst Fit in anderem Kanal!)
        goon=0;
        while goon==0
            pause(5)
        end
        
        answer2 = questdlg('Do you want to use the exponential (red curve) or polynomial (green) correction?', ...
	                            'Correction selection', ...
	                            'Exponential','Polynomial','Exponential');
            % Handle response
            switch answer2
                case 'Exponential'      
                    expocorr = 1;
                case 'Polynomial'
                    expocorr = 0;   
            end
            
          if expocorr==1   
        line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,ff1);
        line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,ff2);
          else
             line1fluorescenceseries_corrected=depletioncorrectionfunc_poly(timeline1,line1fluorescenceseries,correctionfitseg1aa);
        line2fluorescenceseries_corrected=depletioncorrectionfunc_poly(timeline2,line2fluorescenceseries,correctionfitseg2aa);
         
          end
        
        
        
        lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
        line1fluorescenceseries_corrected_binned=mean(reshape(line1fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
        line1fluorescenceseries_corrected_binned=[line1fluorescenceseries_corrected_binned mean(line1fluorescenceseries_corrected(end-lastbinlength+1:end))];
        line2fluorescenceseries_corrected_binned=mean(reshape(line2fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
        line2fluorescenceseries_corrected_binned=[line2fluorescenceseries_corrected_binned mean(line2fluorescenceseries_corrected(end-lastbinlength+1:end))];
        
        figure('OuterPosition',[scrsz(3)/3 scrsz(4)/4 2*scrsz(3)/3 scrsz(4)/4],'Name','Corrected Membrane Fluorescence time series')
        subplot(2,1,1),plot(timeline1_binned,line1fluorescenceseries_corrected_binned); % Diesen Plot auch gebinnt darstellen!
        xlabel('time (linenumber)')
        ylabel('line fluorescence Ch1')
        subplot(2,1,2),plot(timeline2_binned,line2fluorescenceseries_corrected_binned);
        xlabel('time (linenumber)')
        ylabel('line fluorescence Ch2')
        
    else
        if depletioncorrection==1
        %linefluorescenceseries_corrected=depletioncorrectionfunc(timeline,linefluorescenceseries,fsin);
        line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,f1); %Funktion umschreiben, sodass sie auf Fit zur�ckgreift!!
        line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,f2);
        lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
        line1fluorescenceseries_corrected_binned=mean(reshape(line1fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
        line1fluorescenceseries_corrected_binned=[line1fluorescenceseries_corrected_binned mean(line1fluorescenceseries_corrected(end-lastbinlength+1:end))];
        line2fluorescenceseries_corrected_binned=mean(reshape(line2fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
        line2fluorescenceseries_corrected_binned=[line2fluorescenceseries_corrected_binned mean(line2fluorescenceseries_corrected(end-lastbinlength+1:end))];
        
        figure('Name','Corrected Membrane Fluorescence time series')
        subplot(1,2,1),plot(timeline1_binned,line1fluorescenceseries_corrected_binned);
        xlabel('time (linenumber)')
        ylabel('line fluorescence Ch1')
        hold on
        subplot(1,2,2),plot(timeline2_binned,line2fluorescenceseries_corrected_binned);
        xlabel('time (linenumber)')
        ylabel('line fluorescence Ch2')
        end
    end
    
    
    
    % Fourier Spectrum of Fluorescence Time Series
    if fouriercorrection==1
        if depletioncorrection==1
            line1fluorescenceseries=line1fluorescenceseries_corrected;
            line2fluorescenceseries=line2fluorescenceseries_corrected;
        end
        
        for channel=1:2
            if channel==1
                linefluorescenceseries=line1fluorescenceseries;
            else
                linefluorescenceseries=line2fluorescenceseries;
            end
        linefluorescencespec=fft(linefluorescenceseries);
        Y=linefluorescencespec;
        L=length(linefluorescenceseries);
        Fs=1/(scantime);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        freqs = Fs*(0:(L/2))/L;

        freqs2=Fs*(1:L)/L;
        
        envelope=zeros(1,length(freqs)-fourierwindow_env+1);
        freqs_env=zeros(1,length(freqs)-fourierwindow_env+1);
        for i=1:length(freqs)-fourierwindow_env+1
            envelope(i)=mean(P1(i:i+fourierwindow_env-1))+3*std(P1(i:i+fourierwindow_env-1));
            freqs_env(i)=mean(freqs(i:i+fourierwindow_env-1));
        end
        
        
        %envfitfunc=@(x,t) x(1).*exp(-x(2).*t)+x(3).*exp(-x(4).*t)+x(5).*exp(-x(6).*t);
        %x0=[1 1 1 1 10 100];
        %fixed=[false false false false false false];
        if loadfourierenvelope==1
            [inputfilenameenv, path]=uigetfile('*.txt');
            envparameters=load([path '/' inputfilenameenv]);
            envscalefitfunc=@(x,t) x(1).*(envparameters(1).*exp(-envparameters(2).*t)+envparameters(3).*exp(-envparameters(4).*t)+envparameters(5).*exp(-envparameters(6).*t));
            x0FFT=0.5;
            fixedFFT=false;
            [envscalefactor,residuals,J,COVB,MSE] = nlinfitsome(fixedFFT,freqs_env(2:end),envelope(2:end),envscalefitfunc,x0FFT);
            envfit=envscalefitfunc(envscalefactor,freqs);
            
            Ycropped=Y;
            for kkk=2:length(freqs)
                cropfactor=P1(kkk)/envfit(kkk);
                if cropfactor>1
%                     Ycropped(kkk)=Y(kkk)/sqrt(cropfactor);
%                     Ycropped(end-kkk+2)=Y(end-kkk+2)/sqrt(cropfactor);
                    Ycropped(kkk)=Y(kkk)/cropfactor;
                    Ycropped(end-kkk+2)=Y(end-kkk+2)/cropfactor;
                end
            end
            % Cropped Fourier Spectrum
            P2cropped = abs(Ycropped/L);
            P1cropped = P2cropped(1:L/2+1);
            P1cropped(2:end-1) = 2*P1cropped(2:end-1);

            figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Filtered Fourier Spectrum')
%             plot(freqs,P2cropped(1:L/2+1),'.b')
%             hold on
%             plot(freqs,P2(1:L/2+1),'.m')
            plot(freqs,P1cropped,'.b')
            hold on
            plot(freqs,P1,'.m')
            hold on
            plot(freqs,envfit,'-r')
            xlim([0 50]);
            ylim([0 0.3]);
            %plot(linefluorescencespec)
            xlabel('Frequency [Hz]')
            ylabel('Fourier Amplitude')
            legend('Cropped Spectrum','Spectrum','Envelope')
            
            linefluorescenceseries_filtered=ifft(Ycropped);
            
            linefluorescenceseries_corrected=linefluorescenceseries_filtered;
            if channel==1
                line1fluorescenceseries_corrected=linefluorescenceseries_corrected;
            else
                line2fluorescenceseries_corrected=linefluorescenceseries_corrected;
            end 
            lastbinlength=mod(length(linefluorescenceseries),binningwindow)+binningwindow;
            linefluorescenceseries_corrected_binned=mean(reshape(linefluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
            linefluorescenceseries_corrected_binned=[linefluorescenceseries_corrected_binned mean(linefluorescenceseries_corrected(end-lastbinlength+1:end))];
            
            hIFT=figure('OuterPosition',[scrsz(3)/3 scrsz(4)/4 2*scrsz(3)/3 scrsz(4)/4],'Name','Corrected Membrane Fluorescence time series after FT cropping');
            subplot(2,1,channel),plot(timeline1_binned,linefluorescenceseries_corrected_binned); % Diesen Plot auch gebinnt darstellen!
            xlabel('time (linenumber)')
            ylabel('line fluorescence')
            hold on
            if depletioncorrection==1
                if channel==1
                    line1fluorescenceseries=line1fluorescenceseries_corrected;
                else
                    line2fluorescenceseries=line2fluorescenceseries_corrected;
                end
            end
        
            
        else
            envfitfunc=@(x,t) x(1).*exp(-x(2).*t)+x(3).*exp(-x(4).*t)+x(5).*exp(-x(6).*t);
            %x0=[0.1 10^-4 0.01 2 0.01 0.1];
            x0FFT=[10 10^-3 10 2 10 0.1];
            fixedFFT=[false false false false false false];
            [envfitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixedFFT,freqs_env(2:end),envelope(2:end),envfitfunc,x0FFT);

            envfit=envfitfunc(envfitparameters,freqs_env);
            
            %path2=uigetdir('E:\shared\Anne\microscope data\2021-05-21_sFCCS\');        % <------------------------------------------
            fid11=fopen([path2 '\' inputfilename(1:end-4) '_fourierenvelopefit' sprintf('%i',channel) '.txt'],'a'); % adjust path if necessary!
            fid12=fopen([path2 '\' inputfilename(1:end-4) '_fourierenvelope' sprintf('%i',channel) '.txt'],'a'); % adjust path if necessary!
            outputparameters=envfitparameters';
            output=zeros(size(envelope,2),2);
            output(:,1)=freqs_env;
            output(:,2)=envelope;
            fprintf(fid11,'%e\n',outputparameters');
            fprintf(fid12,'%e\t %e\n',output');
        end
                
        
%         return
% 
%         
%         
% 
        hFT=figure('OuterPosition',[scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Full Fourier Spectrum');
        plot(freqs,P1,'.b')
        hold on
        plot(freqs_env,envelope,'.r')
        hold on
        if loadfourierenvelope==1
            plot(freqs,envfit,'-m','LineWidth',2)
        else
            plot(freqs_env,envfit,'-r')
        end
        %plot(linefluorescencespec)
        xlabel('Frequency [Hz]')
        ylabel('Fourier Amplitude')
        xlim([0 5]);
        ylim([0 0.3]);
        
    
        if loadfourierenvelope==0
            savefig(hFT,[path2 '\' inputfilename(1:end-4) ' FT Spectrum' sprintf('%i',channel) '.fig'])
            if channel==2
            return
            end
        end
        end
    end
    
    
    
    

%     % Depletion correction of time trace;
%     if depletioncorrection==1
%         figure('Name','Membrane Fluorescence time series, depletion corrected')
%         % Add Poissonian threshold intensities?
%         line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,f1);
%         line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,f2);
%         subplot(1,2,1),plot(timeline1,line1fluorescenceseries_corrected);
%         xlabel('time (linenumber)')
%         ylabel('line fluorescence')
%         hold on
%         subplot(1,2,2),plot(timeline2,line2fluorescenceseries_corrected);
%         xlabel('time (linenumber)')
%         ylabel('line fluorescence')
%     end
    % Removal of single bright events
    if intensityfilter==1
        if depletioncorrection==1
            line1fluorescenceseries_corrected=brightintensityfilter(line1fluorescenceseries_corrected);
            line2fluorescenceseries_corrected=brightintensityfilter(line2fluorescenceseries_corrected);
            figure('Name','Line fluorescence distribution after removal of bright events')
            % maybe better not as histogram, but time series with threshold values of cut-off intensities!
            subplot(1,2,1),hist(line1fluorescenceseries_corrected,20)
            hold on
            subplot(1,2,2),hist(line2fluorescenceseries_corrected,20)
        else
            line1fluorescenceseries=brightintensityfilter(line1fluorescenceseries);
            line2fluorescenceseries=brightintensityfilter(line2fluorescenceseries);
            figure('Name','Line fluorescence distribution after removal of bright events')
            subplot(1,2,1),hist(line1fluorescenceseries,20)
            hold on
            subplot(1,2,2),hist(line2fluorescenceseries,20)
        end
    end

    % Calculation of autocorrelation function
    % fprintf('Calculating correlation...\n');
    % fcorr=autocorrFCS(linefluorescenceseries);
    % tcorr=(1:1:length(linefluorescenceseries)-1)*linetime;
    fprintf('Calculating correlation using multiple tau...\n');
    if depletioncorrection==1
        [tcorr1,fautocorr1,sigmas1]=autocorrFCSmultipletau(line1fluorescenceseries_corrected);
        [tcorr2,fautocorr2,sigmas2]=autocorrFCSmultipletau(line2fluorescenceseries_corrected);
        if PIE==0
            [tcorrcc,fcrosscorr,sigmascc]=crosscorrFCSmultipletau(line1fluorescenceseries_corrected,line2fluorescenceseries_corrected);
        else
            [tcorrcc,fcrosscorr,sigmascc]=crosscorrFCSmultipletauPIE(line1fluorescenceseries_corrected,line2fluorescenceseries_corrected);
        end
        % Compute cross-correlation function!!!!!!!!!!!!
    else
        [tcorr1,fautocorr1,sigmas1]=autocorrFCSmultipletau(line1fluorescenceseries);
        [tcorr2,fautocorr2,sigmas2]=autocorrFCSmultipletau(line2fluorescenceseries);
        if PIE==0
            [tcorrcc,fcrosscorr,sigmascc]=crosscorrFCSmultipletau(line1fluorescenceseries,line2fluorescenceseries);
        else
             [tcorrcc,fcrosscorr,sigmascc]=crosscorrFCSmultipletauPIE(line1fluorescenceseries,line2fluorescenceseries);
        end
    end
    if PIE==0
        tcorr1=tcorr1*scantime;
        tcorr2=tcorr2*scantime;
        tcorrcc=tcorrcc*scantime;
    else
        tcorr1=tcorr1*scantime;
        tcorr2=tcorr2*scantime;
        tcorrcc=tcorrcc*scantime+linetime;
    end
    
    weights1=fautocorr1./sigmas1;
    weights2=fautocorr2./sigmas2;
    weightscc=fcrosscorr./sigmascc;
    

    % Figure: Plot of autocorrelation function
    figure('OuterPosition',[scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Correlation functions')
    % semilogx(tcorr,fcorr,'b');
    xlabel('time')
    ylabel('autocorrelation')
    % hold on
    semilogx(tcorr1,fautocorr1,'g.')
    hold on
    semilogx(tcorr2,fautocorr2,'r.')
    hold on
    semilogx(tcorrcc,fcrosscorr,'b.')

    % figure('Name','Weights')
    % semilogx(tcorr2,1./sigmas,'b.')
    % Need to plot weights

    %path2=uigetdir;
    fid1=fopen([path2 '\' inputfilename(1:end-4) '.txt'],'a'); % adjust path if necessary!
    %fid3=fopen([path2 '\' inputfilename(1:end-4) '_files.out'],'a');
    
    output=zeros(length(tcorr2),9);
  

    output(:,1)=tcorr1;
    output(:,2)=fautocorr1;
    output(:,3)=weights1;
    output(:,4)=tcorr2;
    output(:,5)=fautocorr2;
    output(:,6)=weights2;
    if PIE==0
        output(:,7)=tcorrcc;
        output(:,8)=fcrosscorr;
        output(:,9)=weightscc;
    else
        output(:,7)=tcorrcc(1:end-1);
        output(:,8)=fcrosscorr(1:end-1);
        output(:,9)=weightscc(1:end-1);
    end
    fprintf(fid1,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output');
    
    
else
    % Modify this part!!!!!!!!!!!!
    % Calculate pooled ACF
    path= uigetdir; 
    files=dir([path '/*.txt']);
    namefile1=files(1).name;
    datafile1=load([path '/' namefile1]);
    ACFlength=size(datafile1,1);
    
    for j=2:3:8
        ACFvalues=zeros(ACFlength,size(files,1));
        timepoints=zeros(ACFlength,size(files,1));
        timepointscc=timepoints;
        weightsdata=zeros(ACFlength,size(files,1));
        for i=1:size(files,1)
            [namedata,remain]=strtok(files(i).name,'.');
            namefile=files(i).name;
            ACFdatai=load([path '/' namefile]);
            timepoints(:,i)=ACFdatai(:,1);
            timepointscc(:,i)=ACFdatai(:,7);
            ACFvalues(:,i)=ACFdatai(:,j);
            weightsdata(:,i)=ACFdatai(:,j+1);
        end
    
        % Calculation of average ACF
        fcorrj=zeros(ACFlength,1);
        tcorrj=timepoints(:,1);
        tcorrjcc=timepointscc(:,1);
        weightsj=zeros(ACFlength,1);
        sigmasj=weightsj;
        for ii=1:size(files,1)
            fcorrj=fcorrj+weightsdata(:,ii).*ACFvalues(:,ii);
            sigmasj=sigmasj+1./weightsdata(:,ii).^2;
        end
        fcorrj=fcorrj./sum(weightsdata,2);
        weightsj=1./sqrt(sigmasj);
        if j==2
            tcorr1=tcorrj';
            fautocorr1=fcorrj';
            weights1=weightsj';
        else
            if j==5
                tcorr2=tcorrj';
                fautocorr2=fcorrj';
                weights2=weightsj';
            else
                tcorrcc=tcorrjcc';
                fcrosscorr=fcorrj';
                weightscc=weightsj';
            end
        end
    end
end

if calculateACFonly==1
    return
end


% Fit of correlation function with diffusion model
fprintf('Fitting...\n');
lb=scantime;
%lb=linetime;
lbcc=lb;
%ub=500;
ub=max(tcorr1);
ubcc=ub;
tfit=tcorr1; % So far nothing done with lb and ub. Adjust later, take care of different lb in case of cross corr func (time axis shifted by linetime)
%tfit=lb:linetime:ub;

% Weights
weights1=abs(weights1);
weights2=abs(weights2);
weightscc=abs(weightscc);
% weights1=1./sigmas1;
% weights2=1./sigmas2;
% weightscc=1./sigmascc;
% weights=zeros(size(tfit));
% for i=1:length(tfit)
%     weights(i)=1/i;
% end    
tcorrfit1=tcorr1(1:length(tfit));
tcorrfit2=tcorr2(1:length(tfit));
tcorrfitcc=tcorrcc;
%tcorrfitcc=tcorrcc(1:length(tfit));

fcorrforfit1=fautocorr1(1:length(tfit));
fcorrforfit2=fautocorr2(1:length(tfit));
fcrosscorrforfit=fcrosscorr;
%fcrosscorrforfit=fcrosscorr(1:length(tfit));
lsautofitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
%x0=[2,0.04 S]; % Set in header or GUI manually later!
%fixed=[false false true];
flscrossfitfunc=lsautofitfunc;
%flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
%y0=[10,0.10,S]; %[N0,tau0,S], Set in header or GUI manually later!

% % For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
% [N1,taud1,Sfit1,CI1,fitcurve1,residuals1] = autocorrfit2Ddiff(tcorrfit1,fcorrforfit1,lb,ub,weights1,lsautofitfunc,x0,fixed);
% [N2,taud2,Sfit2,CI2,fitcurve2,residuals2] = autocorrfit2Ddiff(tcorrfit2,fcorrforfit2,lb,ub,weights2,lsautofitfunc,x0,fixed);
% [Ncc,taucc,Sfitcc,CI2f,fitcurve2f,residuals2f] =autocorrfit2Ddiff(tcorrfitcc,fcrosscorrforfit,lbcc,ubcc,weightscc,flscrossfitfunc,y0,fixed);
% %same for cross-correlation function and corresponding vectors!
% 
% % Figure: Plot of data and fit
% %N=10;
% %taud=0.002;
% figure('OuterPosition',[scrsz(3) scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Fit curves')
% subplot(2,1,1),semilogx(tcorrfit1,fitcurve1,'-g')
% hold on
% subplot(2,1,1),semilogx(tcorrfit2,fitcurve2,'-r')
% hold on
% subplot(2,1,1),semilogx(tcorrfitcc,fitcurve2f,'-b')
% hold on
% subplot(2,1,1),semilogx(tcorrfit1,fcorrforfit1,'gs')
% hold on
% subplot(2,1,1),semilogx(tcorrfit2,fcorrforfit2,'rd')
% hold on
% subplot(2,1,1),semilogx(tcorrfitcc,fcrosscorrforfit,'bx')
% xlabel('time')
% ylabel('autocorrelation')
% legend('auto Ch1','auto Ch2','2color cross')
% hold on
% subplot(2,1,2),semilogx(tcorrfit1,residuals1,'gs')
% hold on
% subplot(2,1,2),semilogx(tcorrfit2,residuals2,'rd')
% hold on
% subplot(2,1,2),semilogx(tcorrfitcc,residuals2f,'bx')
% xlabel('time')
% ylabel('residuals')
% hold on
% subplot(2,1,2),semilogx(tcorrfit1,zeros(size(tcorrfit1)),'-k');
% legend('auto Ch1','auto Ch2','2color cross')

%%% From Here: Keep modifying for 2 Channels!!!

if depletioncorrection==1
    line1fluorescenceseries=line1fluorescenceseries_corrected;
    line2fluorescenceseries=line2fluorescenceseries_corrected;
end

lastbinlength=mod(length(line1fluorescenceseries),binningwindow_GUI)+binningwindow_GUI;
line1fluorescenceseriesbinned_GUI=mean(reshape(line1fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line1fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line1fluorescenceseriesbinned_GUI=[line1fluorescenceseriesbinned_GUI mean(line1fluorescenceseries(end-lastbinlength+1:end))];
line2fluorescenceseriesbinned_GUI=mean(reshape(line2fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line2fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line2fluorescenceseriesbinned_GUI=[line2fluorescenceseriesbinned_GUI mean(line2fluorescenceseries(end-lastbinlength+1:end))];

if ACFselection==1
    segmentlength=size(line1array,1)/ACFsegmentnumb;
    Isegmentlength_GUI=segmentlength/binningwindow_GUI;
    Itrace1=zeros(segmentlength,ACFsegmentnumb);
    Itrace2=Itrace1;
    Itrace1_GUI=zeros(Isegmentlength_GUI,ACFsegmentnumb);
    Itrace2_GUI=Itrace1_GUI;
     Itrace1_poly=Itrace1;
    Itrace1_poly_GUI=Itrace1_GUI;
     Itrace2_poly=Itrace2;
    Itrace2_poly_GUI=Itrace2_GUI;
    
    % Fitting average curve in GUI
%    ?? correlationcurves=zeros(1,ACFsegmentnumb);
%    ?? sigmascurves=correlationcurves;
%    ?? corfit=zeros(size(correlationcurves,1),size(correlationcurves,2));
    
    figure('OuterPosition',[scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','ACFs segments - Press SPACE')
    %x0=[N,taud,Sfit];
    for i=1:size(Itrace1,2)
        Itrace1(:,i)=line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)';
        Itrace2(:,i)=line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)';
        Itrace1_GUI(:,i)=line1fluorescenceseriesbinned_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI)';
        Itrace2_GUI(:,i)=line2fluorescenceseriesbinned_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI)';
        assex=0:scantime:size(Itrace1,1)*scantime-scantime;
        assex_GUI=mean(reshape(assex, binningwindow_GUI, size(assex,2)/binningwindow_GUI),1); 
        
        
     fitting_poly_GUI1(1:polygrad+1,i) = polyfit(assex_GUI',Itrace1_GUI(:,i),polygrad);
        fitting_poly_GUI2(1:polygrad+1,i) = polyfit(assex_GUI',Itrace2_GUI(:,i),polygrad);

        %%% older version
        if Jonas==0
            Itrace1_poly(:,i)=Itrace1(:,i)-(polyval(fitting_poly_GUI1(:,i),assex))'+mean(polyval(fitting_poly_GUI1(:,i),assex));
        Itrace1_poly_GUI(:,i)=Itrace1_GUI(:,i)-(polyval(fitting_poly_GUI1(:,i),assex_GUI))'+mean(polyval(fitting_poly_GUI1(:,i),assex_GUI));
        Itrace2_poly(:,i)=Itrace2(:,i)-(polyval(fitting_poly_GUI2(:,i),assex))'+mean(polyval(fitting_poly_GUI2(:,i),assex));
        Itrace2_poly_GUI(:,i)=Itrace2_GUI(:,i)-(polyval(fitting_poly_GUI2(:,i),assex_GUI))'+mean(polyval(fitting_poly_GUI2(:,i),assex_GUI));
        end
        
        %%%%%
        %    Jonas version of the correction
        if Jonas==1
            effediti1=(polyval(fitting_poly_GUI1(:,i),assex))';
            effediti1_GUI=(polyval(fitting_poly_GUI1(:,i),assex_GUI))';
            Itrace1_poly(:,i)=Itrace1(:,i)./(effediti1./effediti1(1,1)).^0.5+effediti1(1,1).*(1-(effediti1./effediti1(1,1)).^0.5)   ;
            Itrace1_poly_GUI(:,i)=Itrace1_GUI(:,i)./(effediti1_GUI./effediti1_GUI(1,1)).^0.5+effediti1_GUI(1,1).*(1-(effediti1_GUI./effediti1_GUI(1,1)).^0.5)   ;
          effediti2=(polyval(fitting_poly_GUI2(:,i),assex))';
            effediti2_GUI=(polyval(fitting_poly_GUI2(:,i),assex_GUI))';
            Itrace2_poly(:,i)=Itrace2(:,i)./(effediti2./effediti2(1,1)).^0.5+effediti2(1,1).*(1-(effediti2./effediti2(1,1)).^0.5)   ;
            Itrace2_poly_GUI(:,i)=Itrace2_GUI(:,i)./(effediti2_GUI./effediti2_GUI(1,1)).^0.5+effediti2_GUI(1,1).*(1-(effediti2_GUI./effediti2_GUI(1,1)).^0.5)   ;
        
        end
        %%%%%%%%%%%%%%%%%%
        
        
        
        
        [tcorr1i,fcorr1i,sigmas1i]=autocorrFCSmultipletau(Itrace1(:,i));
        [tcorr2i,fcorr2i,sigmas2i]=autocorrFCSmultipletau(Itrace2(:,i));
          [tcorr1i_poly,fcorr1i_poly,sigmas1i_poly]=autocorrFCSmultipletau(Itrace1_poly(:,i));
        [tcorr2i_poly,fcorr2i_poly,sigmas2i_poly]=autocorrFCSmultipletau(Itrace2_poly(:,i));
        if PIE==0
            [tcorrcci,fcrosscorri,sigmascci]=crosscorrFCSmultipletau(Itrace1(:,i),Itrace2(:,i));
        [tcorrcci_poly,fcrosscorri_poly,sigmascci_poly]=crosscorrFCSmultipletau(Itrace1_poly(:,i),Itrace2_poly(:,i));

        else
            [tcorrcci_poly,fcrosscorri_poly,sigmascci_poly]=crosscorrFCSmultipletauPIE(Itrace1_poly(:,i),Itrace2_poly(:,i));
 [tcorrcci,fcrosscorri,sigmascci]=crosscorrFCSmultipletauPIE(Itrace1(:,i),Itrace2(:,i));

        end
        if PIE==0
            tcorr1i=tcorr1i*scantime;
            tcorr2i=tcorr2i*scantime;
            tcorrcci=tcorrcci*scantime;
              tcorr1i_poly=tcorr1i_poly*scantime;
            tcorr2i_poly=tcorr2i_poly*scantime;
            tcorrcci_poly=tcorrcci_poly*scantime;
        else
            tcorr1i_poly=tcorr1i_poly*scantime;
            tcorr2i_poly=tcorr2i_poly*scantime;
            tcorrcci_poly=tcorrcci_poly*scantime+linetime;
                  tcorr1i=tcorr1i*scantime;
            tcorr2i=tcorr2i*scantime;
            tcorrcci=tcorrcci*scantime+linetime;
        end
        correlationcurvesCh1(:,i)=fcorr1i';
        sigmas1curves(:,i)=sigmas1i';
        correlationcurvesCh2(:,i)=fcorr2i';
        sigmas2curves(:,i)=sigmas2i';
        correlationcurvesChCC(:,i)=fcrosscorri';
        sigmascccurves(:,i)=sigmascci';
        lbi=min(tcorr1i);
        ubi=max(tcorr1i);
        
              correlationcurvesCh1_poly(:,i)=fcorr1i_poly';
        sigmas1curves_poly(:,i)=sigmas1i_poly';
        correlationcurvesCh2_poly(:,i)=fcorr2i_poly';
        sigmas2curves_poly(:,i)=sigmas2i_poly';
        correlationcurvesChCC_poly(:,i)=fcrosscorri_poly';
        sigmascccurves_poly(:,i)=sigmascci_poly';
        lbi_poly=min(tcorr1i_poly);
        ubi_poly=max(tcorr1i_poly);
% Fit ACFs of segments
         %[Ni,taudi,Sfiti,CIi,fitcurvei,residualsi] = autocorrfit2Ddiff(tcorri*linetime,fcorri,lbi,ubi,1./sigmasi,lsfitfunc,x0);
         %x0=[Ni,taudi,Sfiti];
         %corfit(:,i)=fitcurvei';
         
         %[N1i,taud1i,Sfit1i,CI1i,fitcurve1i,residuals1i] = autocorrfit2Ddiff(tcorr1i*linetime,fcorr1i,lbi,ubi,1./sigmas1i,lsfitfunc,x0);
         %x0=[N1i,taud1i,Sfit1i];
         %corfit(:,i)=fitcurve1i';

%         subplot(2,2,1),semilogx(tcorr1i,fcorr1i,'.g')
%         hold on
%         subplot(2,2,2),semilogx(tcorr2i,fcorr2i,'.r')
%         hold on
%         subplot(2,2,3),semilogx(tcorrcci,fcrosscorri,'.b')
        
        %p1=subplot(2,2,3),semilogx(tcorr1i,fcorr1i,'Color',[44/255 123/255 182/255],'LineStyle','none','Marker','.','MarkerFaceColor',[44/255 123/255 182/255]);
        p1=subplot(2,2,3),title('ACF_1'),xlabel('time'),ylabel('autocorrelation'),semilogx(tcorr1i,fcorr1i,'Color',[44/255 123/255 182/255],'LineStyle','none','Marker','.','MarkerFaceColor',[44/255 123/255 182/255]);
        %p1(1).Marker = '.';
        hold on
   p1_poly=subplot(2,2,3),title('ACF_1'),xlabel('time'),ylabel('autocorrelation'),semilogx(tcorr1i_poly,fcorr1i_poly,'.r');

        %p2=subplot(2,2,4),semilogx(tcorr2i,fcorr2i,'Color',[44/255 123/255 182/255],'LineStyle','none','Marker','.','MarkerFaceColor',[44/255 123/255 182/255]);
        p2=subplot(2,2,4),title('ACF_2'),xlabel('time'),ylabel('autocorrelation'),semilogx(tcorr2i,fcorr2i,'Color',[44/255 123/255 182/255],'LineStyle','none','Marker','.','MarkerFaceColor',[44/255 123/255 182/255]);
        hold on
           p2_poly=subplot(2,2,4),title('ACF_2'),xlabel('time'),ylabel('autocorrelation'),semilogx(tcorr2i_poly,fcorr2i_poly,'.r');
%p2(1).Marker = '.';
        %p3=subplot(2,2,2),semilogx(tcorrcci,fcrosscorri,'Color', [44/255 123/255 182/255],'LineStyle','none','Marker','.','MarkerFaceColor',[44/255 123/255 182/255]);
        p3=subplot(2,2,2),title('CCF'),xlabel('time'),ylabel('crosscorrelation'),semilogx(tcorrcci,fcrosscorri,'Color', [44/255 123/255 182/255],'LineStyle','none','Marker','.','MarkerFaceColor',[44/255 123/255 182/255]);
        hold on
           p3_poly=subplot(2,2,2),title('CCF'),xlabel('time'),ylabel('crosscorrelation'),semilogx(tcorrcci_poly,fcrosscorri_poly,'.r');
%p3(1).Marker = '.';
        %semilogx(tcorri,fitcurvei,'-r')
        pause
    end
    corfit1=zeros(length(tcorr1i),ACFsegmentnumb);
    corfit2=corfit1;
        corfit1_poly=corfit1;
    corfit2_poly=corfit1;

    if PIE==0
        corfitcc=corfit1;
        corfitcc_poly=corfit1;
    else
        corfitcc=zeros(size(corfit1,1)+1,size(corfit1,2));
                corfitcc_poly=zeros(size(corfit1,1)+1,size(corfit1,2));

    end
    
%     % Binning for better display
%     Binfactor=100;
%     ItraceII=zeros(size(Itrace,1)/Binfactor,size(Itrace,2));
%     for i=1:size(Itrace,2)
%         ItraceII(:,i)=nanmean(reshape(Itrace(:,i),Binfactor,length(Itrace(:,1))/Binfactor),1);
%     end

ItraceIICh1=Itrace1_GUI;
ItraceIICh2=Itrace2_GUI;
%    segmenttime=linetime:linetime:linetime*segmentlength;
%    ItraceII=[segmenttime' ItraceII];



%     ItraceIICh1=Itrace1;
%     ItraceIICh2=Itrace2;
    if PIE==0
%         segmenttime1=scantime:scantime:scantime*segmentlength;
%         segmenttime2=segmenttime1;
        segmenttime1_GUI=scantime:scantime*binningwindow_GUI:scantime*segmentlength;
        segmenttime2_GUI=segmenttime1_GUI;
        %segmenttimecc=segmenttime1;
    else
%         segmenttime1=scantime:scantime:scantime*segmentlength;
%         segmenttime2=scantime:scantime:scantime*segmentlength;
        segmenttime1_GUI=scantime:scantime*binningwindow_GUI:scantime*segmentlength;
        segmenttime2_GUI=scantime:scantime*binningwindow_GUI:scantime*segmentlength; % Is this correct?
        %segmenttimecc=scantime-linetime:scantime:scantime*segmentlength-linetime;
    end
%     ItraceIICh1=[segmenttime1' ItraceIICh1];
%     ItraceIICh2=[segmenttime2' ItraceIICh2];
    
    ItraceIICh1=[segmenttime1_GUI' ItraceIICh1];
    ItraceIICh2=[segmenttime2_GUI' ItraceIICh2];
     ItraceIICh1_poly= ItraceIICh1;
        ItraceIICh2_poly= ItraceIICh2;

    
%     segmenttime_binned=linetime:Binfactor*linetime:linetime*segmentlength;
%     ItraceII=[segmenttime_binned' ItraceII];
    correlationcurvesCh1=[tcorr1i' correlationcurvesCh1];

    correlationcurvesCh2=[tcorr2i' correlationcurvesCh2];
    correlationcurvesChCC=[tcorrcci' correlationcurvesChCC];
    corfit1=[tcorr1i' corfit1];
    corfit2=[tcorr2i' corfit2];
    corfitcc=[tcorrcci' corfitcc];
    curveincl1=ones(size(corfit1,2)-1,1);
    curveincl2=curveincl1;
    
        correlationcurvesCh1_poly=[tcorr1i' correlationcurvesCh1_poly];

     correlationcurvesCh2_poly=[tcorr2i' correlationcurvesCh2_poly];
    correlationcurvesChCC_poly=[tcorrcci' correlationcurvesChCC_poly];
    corfit1_poly=[tcorr1i' corfit1_poly];
    corfit2_poly=[tcorr2i' corfit2_poly];
    corfitcc_poly=[tcorrcci' corfitcc_poly];
    curveincl1_poly=zeros(size(corfit1_poly,2)-1,1);
    curveincl2_poly=curveincl1_poly;

    % GUI for selection of segment ACFs
    % =====================================================================
    % ============== Modify axis of windows!!!!!!!!============
    % =====================================================================
    
    %Corrselection2Channelsindividuell
    Corrselection2Channelsindividuellpreview_new
    goon=0;
    while goon==0
        pause(5)
    end
    fprintf('Done!')
    for i=1:length(curveincl1)
     
         if curveincl1_poly(i)==1
            line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)-(polyval(fitting_poly_GUI1(:,i),assex))'+mean(polyval(fitting_poly_GUI1(:,i),assex));
            %Jonas version
            if Jonas==1
                effediti=(polyval(fitting_poly_GUI1(:,i),assex))';
                line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)./(effediti./effediti(1,1)).^0.5+effediti(1,1).*(1-(effediti./effediti(1,1)).^0.5);
            end
        end
           if curveincl1(i)==0
            line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=NaN;
           end
        
           
        if curveincl2_poly(i)==1
            line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)-(polyval(fitting_poly_GUI2(:,i),assex))'+mean(polyval(fitting_poly_GUI2(:,i),assex));
            %Jonas version
            if Jonas==1
                effediti=(polyval(fitting_poly_GUI1(:,i),assex))';
                line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)./(effediti./effediti(1,1)).^0.5+effediti(1,1).*(1-(effediti./effediti(1,1)).^0.5);
            end
        end
        
        
        
        if curveincl2(i)==0
            line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=NaN;
        end
    end
    
    % Calculation of final ACF
    %need to check whether multipletau scripts (also for CC) work correctly
    %for normalization when segments are removed. Also, why do they give
    %NaN sometime?
    [tcorr1final,fcorr1final,sigmas1final]=autocorrFCSmultipletau(line1fluorescenceseries);
    [tcorr2final,fcorr2final,sigmas2final]=autocorrFCSmultipletau(line2fluorescenceseries);
    if PIE==0
        [tcorrccfinal,fcrosscorrfinal,sigmasccfinal]=crosscorrFCSmultipletau(line1fluorescenceseries,line2fluorescenceseries);
    else
        [tcorrccfinal,fcrosscorrfinal,sigmasccfinal]=crosscorrFCSmultipletauPIE(line1fluorescenceseries,line2fluorescenceseries);
    end
    
    if PIE==0
        tcorr1final=tcorr1final*scantime;
        tcorr2final=tcorr2final*scantime;
        tcorrccfinal=tcorrccfinal*scantime;
    else
        tcorr1final=tcorr1final*scantime;
        tcorr2final=tcorr2final*scantime;
        tcorrccfinal=tcorrccfinal*scantime+linetime;
    end
    
     weights1final=abs(fcorr1final./sigmas1final);
    weights2final=abs(fcorr2final./sigmas2final);
    weightsccfinal=abs(fcrosscorrfinal./sigmasccfinal);
    
     % ------> average of segments
    fcorr1segfinal=zeros(length(tcorr1i),1);
    sigmas1segfinal=fcorr1segfinal;
    weight1seg=fcorr1segfinal;
    fcorr2segfinal=zeros(length(tcorr2i),1);
    sigmas2segfinal=fcorr2segfinal;
    weight2seg=fcorr2segfinal;
    
    fcrosscorrsegfinal=zeros(length(tcorrcci),1);
    sigmasccsegfinal=fcrosscorrsegfinal;
    weightsccseg=fcrosscorrsegfinal;
    for i=1:length(curveincl1)
        if curveincl1(i)==1
            if curveincl1_poly(i)==0
             fcorr1segfinal=fcorr1segfinal+correlationcurvesCh1(:,1+i);
             
            %fcorrsegfinal=fcorrsegfinal+correlationcurves(:,i+1)./sigmascurves(:,i);
             sigmas1segfinal=sigmas1segfinal+sigmas1curves(:,i);
            else
                 fcorr1segfinal=fcorr1segfinal+correlationcurvesCh1_poly(:,1+i);
             
            %fcorrsegfinal=fcorrsegfinal+correlationcurves(:,i+1)./sigmascurves(:,i);
             sigmas1segfinal=sigmas1segfinal+sigmas1curves_poly(:,i);
            %sigmassegfinal=sigmassegfinal+sigmascurves(:,i);
            %weightsseg=weightseg+1./sigmascurves(:,i);
            end
        end
        if curveincl2(i)==1
            if curveincl2_poly(i)==0
            fcorr2segfinal=fcorr2segfinal+correlationcurvesCh2(:,1+i);
            sigmas2segfinal=sigmas2segfinal+sigmas2curves(:,i);
            else
                fcorr2segfinal=fcorr2segfinal+correlationcurvesCh2_poly(:,1+i);
            sigmas2segfinal=sigmas2segfinal+sigmas2curves_poly(:,i);
            end
        end
        
        
        if curveincl1(i) && curveincl2(i)
            if curveincl1_poly(i) || curveincl2_poly(i)
        fcrosscorrsegfinal=fcrosscorrsegfinal+correlationcurvesChCC_poly(:,1+i);
        sigmasccsegfinal=sigmasccsegfinal+sigmascccurves_poly(:,i);
            else
                fcrosscorrsegfinal=fcrosscorrsegfinal+correlationcurvesChCC(:,1+i);
        sigmasccsegfinal=sigmasccsegfinal+sigmascccurves(:,i);
            end
            
        end
    end
    fcorr1segfinal=fcorr1segfinal/sum(curveincl1);
    fcorr2segfinal=fcorr2segfinal/sum(curveincl2);
    fcrosscorrsegfinal=fcrosscorrsegfinal/sum((curveincl1+curveincl2)==2);
    %fcorrsegfinal=fcorrsegfinal./weightsseg;
    sigmas1segfinal=sigmas1segfinal/sum(curveincl1);
    sigmas2segfinal=sigmas2segfinal/sum(curveincl2);
    sigmasccsegfinal=sigmasccsegfinal/sum((curveincl1+curveincl2)==2);
    
        
    %fcorrsegfinal=nansum(correlationcurves(:,2:end)./sigmascurves,2);
    %fcorrsegfinal=fcorrsegfinal./nansum(sigmascurves,2);
    %sigmassegfinal=nanmean(sigmascurves,2);
    weights1segfinal=abs(fcorr1segfinal./sigmas1segfinal);
    tcorr1segfinal=tcorr1i;
    weights2segfinal=abs(fcorr2segfinal./sigmas2segfinal);
    tcorr2segfinal=tcorr2i;
    weightsccsegfinal=abs(fcrosscorrsegfinal./sigmasccsegfinal);
    tcorrccsegfinal=tcorrcci;
    
    
    % Fit of final correlation function with diffusion model
    fprintf('Final Fitting...\n');
    lb=scantime;
    lbcc=lb;
    lbseg=lb;
    lbccseg=lb;
    %ub=1;
    ub=max(tcorr1final);
    ubcc=ub;
    tfit=tcorr1final;
    ubseg=max(tcorr1segfinal);
    ubccseg=ubseg;
    %tfit=lb:linetime:ub;
    %tfit=tcorr2(1:40);

    % Weights
    %weightsfinal=abs(fcorrfinal./sigmasfinal);
    % weights=zeros(size(tfit));
    % for i=1:length(tfit)
    %     weights(i)=1/i;
    % end    
    tcorr1fit=tcorr1final(1:length(tfit));
    fcorr1fit=fcorr1final(1:length(tfit));
    weights1fit=weights1final(1:length(tfit));
    weights1fit(weights1fit==0)=10^-4;
    
    tcorr2fit=tcorr2final(1:length(tfit));
    fcorr2fit=fcorr2final(1:length(tfit));
    weights2fit=weights2final(1:length(tfit));
    weights2fit(weights2fit==0)=10^-4;
    
    tcorrccfit=tcorrccfinal(1:length(tfit));
    fcrosscorrfit=fcrosscorrfinal(1:length(tfit));
    weightsccfit=weightsccfinal(1:length(tfit));
    weightsccfit(weightsccfit==0)=10^-4;
    tcorrccsegfinalfit=tcorrccsegfinal(1:length(tcorr1segfinal));
    fcrosscorrsegfinalfit=fcrosscorrsegfinal(1:length(tcorr1segfinal));
    weightsccsegfinalfit=weightsccsegfinal(1:length(tcorr1segfinal));
    
    lsautofitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
    flscrossfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);
    
    
    fixedcc=[false false true false];
    %flscrossfitfunc=lsautofitfunc;
    %flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
    %y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
    
     weights1fit(isnan(fcorr1fit))=10^-10;
    weights2fit(isnan(fcorr2fit))=10^-10;
    weightsccfit(isnan(fcrosscorrfit))=10^-10;
    weights1fit(isnan(weights1fit))=10^-10;
        weights2fit(isnan(weights2fit))=10^-10;
            weightsccfit(isnan(weightsccfit))=10^-10;
     weights1fit(isinf(weights1fit))=10^-10;
        weights2fit(isinf(weights2fit))=10^-10;
            weightsccfit(isinf(weightsccfit))=10^-10;
    
    fcorr1fit(isnan(fcorr1fit))=0;
fcorr2fit(isnan(fcorr2fit))=0;
fcorrccfit(isnan(fcrosscorrfit))=0;
        
   try
    y0=[1000,0.1,S]; %[N0,tau0,S], Set in header or GUI manually later!
    

    % For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
    %fcorr1fit & co should not have NaN. I am checking the
    %multipletauscript. In the meantime, trying to neglect NaN
    
       

    
    %%%%%%
    [N1final,taud1final,Sfit1final,CI1final,fitcurve1final,residuals1final] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfunc,x0,fixed);
    [N2final,taud2final,Sfit2final,CI2final,fitcurve2final,residuals2final] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfunc,x0,fixed);
    %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
    [Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfunc,y0,fixed);
    
    [N1segfinal,taud1segfinal,Sfit1segfinal,CI1segfinal,fitcurve1segfinal,residuals1segfinal] = autocorrfit2Ddiff(tcorr1segfinal,fcorr1segfinal',lbseg,ubseg,weights1segfinal',lsautofitfunc,x0,fixed);
    [N2segfinal,taud2segfinal,Sfit2segfinal,CI2segfinal,fitcurve2segfinal,residuals2segfinal] = autocorrfit2Ddiff(tcorr2segfinal,fcorr2segfinal',lbseg,ubseg,weights2segfinal',lsautofitfunc,x0,fixed);
    %[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);
    [Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiff(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',lsautofitfunc,y0,fixed);
    errorfull=0
    catch
       
     
        errorfull=1

%     % For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
%     [N1final,taud1final,Sfit1final,CI1final,fitcurve1final,residuals1final] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfunc,x0,fixed);
%     [N2final,taud2final,Sfit2final,CI2final,fitcurve2final,residuals2final] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfunc,x0,fixed);
%     %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
%     [Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfunc,y0,fixed);
%     
try
    [N1segfinal,taud1segfinal,Sfit1segfinal,CI1segfinal,fitcurve1segfinal,residuals1segfinal] = autocorrfit2Ddiff(tcorr1segfinal,fcorr1segfinal',lbseg,ubseg,weights1segfinal',lsautofitfunc,x0,fixed);
catch
     x0=[1,0.01,S]; %[N0,tau0,S], Set in header or GUI manually later!
        [N1segfinal,taud1segfinal,Sfit1segfinal,CI1segfinal,fitcurve1segfinal,residuals1segfinal] = autocorrfit2Ddiff(tcorr1segfinal,fcorr1segfinal',lbseg,ubseg,weights1segfinal',lsautofitfunc,x0,fixed);
end

try
        x0=[1000,0.1,S]; %[N0,tau0,S], Set in header or GUI manually later!

    [N2segfinal,taud2segfinal,Sfit2segfinal,CI2segfinal,fitcurve2segfinal,residuals2segfinal] = autocorrfit2Ddiff(tcorr2segfinal,fcorr2segfinal',lbseg,ubseg,weights2segfinal',lsautofitfunc,x0,fixed);
catch
     x0=[5000,10,S]; %[N0,tau0,S], Set in header or GUI manually later!
    [N2segfinal,taud2segfinal,Sfit2segfinal,CI2segfinal,fitcurve2segfinal,residuals2segfinal] = autocorrfit2Ddiff(tcorr2segfinal,fcorr2segfinal',lbseg,ubseg,weights2segfinal',lsautofitfunc,x0,fixed);
end

try
        y0=[1000,0.1,S]; %[N0,tau0,S], Set in header or GUI manually later!

 %[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);
    [Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiff(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',lsautofitfunc,y0,fixed);
   catch
     y0=[5000,10,S]; %[N0,tau0,S], Set in header or GUI manually later!
    [Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiff(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',lsautofitfunc,y0,fixed);
end

    end
%CI1segfinal=CI2segfinal*0;    
    % Relative cross-correlation amplitude
    if errorfull==0
        relCC=(N1final+N2final)/(2*Nccfinal);
        us1=0.5*(CI1final(1:3,2)-CI1final(1:3,1));
        us2=0.5*(CI2final(1:3,2)-CI2final(1:3,1));
        uscc=0.5*(CIccfinal(1:3,2)-CIccfinal(1:3,1));
    end
    relCCseg=(N1segfinal+N2segfinal)/(2*Nccsegfinal);
    relCCsegnew=max([N1segfinal/Nccsegfinal N2segfinal/Nccsegfinal]);
    
    us1seg=0.5*(CI1segfinal(1:3,2)-CI1segfinal(1:3,1));
    us2seg=0.5*(CI2segfinal(1:3,2)-CI2segfinal(1:3,1));
    usccseg=0.5*(CIccsegfinal(1:3,2)-CIccsegfinal(1:3,1));
    
    % Figure: Plot of data and fit
    % ================ Beschriftung in Graphen!!!!========================
    % Und Daten und Graphen abspeichern!
    %N=10;
    %taud=0.002;
    if errorfull==0
        h=figure('OuterPosition',[1 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
        positionvector1=[0.1 0.35 0.8 0.55];
        positionvector2=[0.1 0.1 0.8 0.15];
        subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1final,'-g',tcorr2fit,fitcurve2final,'-r',tcorrccfit,fitcurveccfinal,'-b')
        hold on
        subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
        xlabel('time')
        ylabel('autocorrelation')
        subplot('Position',positionvector2),semilogx(tcorr1fit,residuals1final,'-gs',tcorr2fit,residuals2final,'-rd',tcorrccfit,residualsccfinal,'-bx')
        xlabel('time')
        
        ylabel('residuals')
        hold on
        subplot('Position',positionvector2),semilogx(tcorr1fit,zeros(size(tcorr1fit)),'-k');
        textN1=['N_1 = %.1f ' char(177) ' %.1f'];
        strtextN1=sprintf(textN1,N1final,us1(1));
        textN2=['N_2 = %.1f ' char(177) ' %.1f'];
        strtextN2=sprintf(textN2,N2final,us2(1));
        textrelcc='rel.cc. = %.2f ';
        strtextrelcc=sprintf(textrelcc,relCC);
         textemprelcc='emp. rel.cc. = %.2f ';
        strtextemprelcc=sprintf(textemprelcc,max([N1segfinal*mean(fcrosscorrfit(1:5)) N2segfinal*mean(fcrosscorrfit(1:5))]) );
        texttaud1=[' = %.3f ' char(177) ' %.3f'];
        strtexttaud1=['\tau_1' sprintf(texttaud1,taud1final,us1(2)) ' s'] ;
        texttaud2=[' = %.3f ' char(177) ' %.3f'];
        strtexttaud2=['\tau_2' sprintf(texttaud2,taud2final,us2(2)) ' s'] ;
        texttaudcc=[' = %.3f ' char(177) ' %.3f'];
        strtexttaudcc=['\tau_{cc}' sprintf(texttaudcc,tauccfinal,uscc(2)) ' s'] ;
        %textSseg=['S = %.1f ' char(177) ' %.1f'];
        %strtextSseg=sprintf(textSseg,Sfit1segfinal,us1seg(3));
    %     textw0=['w0 = %.2f ' char(177) ' %.2f'];
    %     strtextw0=[sprintf(textw0,w0final,us(3)) ' \mu m'];
        dim = [0.68 0.58 0.3 0.3];
    %     str = {strtextC,strtextD,strtextw0};
        str = {strtextN1,strtextN2,strtextrelcc,strtextemprelcc, strtexttaud1,strtexttaud2,strtexttaudcc};
        t=annotation('textbox',dim,'String',str,'FitBoxToText','on');
        set(t,'FontSize',9)
    end
      % Figure: Plot of data and fit
    %N=10;
    %taud=0.002;
    hh=figure('OuterPosition',[800 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
    positionvector1=[0.1 0.35 0.8 0.55];
    positionvector2=[0.1 0.1 0.8 0.15];
    subplot('Position',positionvector1),semilogx(tcorr1segfinal,fitcurve1segfinal,'-g',tcorr2segfinal,fitcurve2segfinal,'-r',tcorrccsegfinalfit,fitcurveccsegfinal,'-b')
    hold on
    subplot('Position',positionvector1),semilogx(tcorr1segfinal,fcorr1segfinal,'gs',tcorr2segfinal,fcorr2segfinal,'rd',tcorrccsegfinalfit,fcrosscorrsegfinalfit,'bx')
    xlabel('time')
    ylabel('autocorrelation')
    subplot('Position',positionvector2),semilogx(tcorr1segfinal,residuals1segfinal,'-gs',tcorr2segfinal,residuals2segfinal,'-rd',tcorrccsegfinalfit,residualsccsegfinal,'-bx')
    xlabel('time')
    ylabel('residuals')
    hold on
    subplot('Position',positionvector2),semilogx(tcorr1segfinal,zeros(size(tcorr1segfinal)),'-k');
    textN1seg=['N_1 = %.1f ' char(177) ' %.1f'];
    strtextN1seg=sprintf(textN1seg,N1segfinal,us1seg(1));
    textN2seg=['N_2 = %.1f ' char(177) ' %.1f'];
    strtextN2seg=sprintf(textN2seg,N2segfinal,us2seg(1));
    textrelccseg='rel.cc. = %.2f ';
    strtextrelccseg=sprintf(textrelccseg,relCCsegnew);
      textrelcc='emp. rel.cc. = %.2f ';
        strtextemprelcc=sprintf(textrelcc,max([N1segfinal*mean(fcrosscorrsegfinalfit(1:5)) N2segfinal*mean(fcrosscorrsegfinalfit(1:5))]) );
       
    texttaud1seg=[' = %.3f ' char(177) ' %.3f'];
    strtexttaud1seg=['\tau_1' sprintf(texttaud1seg,taud1segfinal,us1seg(2)) ' s'] ;
    texttaud2seg=[' = %.3f ' char(177) ' %.3f'];
    strtexttaud2seg=['\tau_2' sprintf(texttaud2seg,taud2segfinal,us2seg(2)) ' s'] ;
    texttaudccseg=[' = %.3f ' char(177) ' %.3f'];
    strtexttaudccseg=['\tau_{cc}' sprintf(texttaudccseg,tauccsegfinal,usccseg(2)) ' s'] ;
    %textSseg=['S = %.1f ' char(177) ' %.1f'];
    %strtextSseg=sprintf(textSseg,Sfit1segfinal,us1seg(3));
%     textw0=['w0 = %.2f ' char(177) ' %.2f'];
%     strtextw0=[sprintf(textw0,w0final,us(3)) ' \mu m'];
    dim = [0.68 0.58 0.3 0.3];
%     str = {strtextC,strtextD,strtextw0};
    strseg = {strtextN1seg,strtextN2seg,strtextrelccseg, strtextemprelcc, strtexttaud1seg,strtexttaud2seg,strtexttaudccseg};
    t=annotation('textbox',dim,'String',strseg,'FitBoxToText','on');
    set(t,'FontSize',9)
    
   
    membranetime_Ch1=membranewidth_Ch1*pixeltime;
    membranetime_Ch2=membranewidth_Ch2*pixeltime;
    
    % path2=uigetdir;
    fid1=fopen([path2 '\' inputfilename(1:end-4) '_final_ACF.txt'],'a'); % adjust path if necessary!
    fid2=fopen([path2 '\' inputfilename(1:end-4) '_final_fitparameters.txt'],'a'); % adjust path if necessary!
    fid3=fopen([path2 '\' inputfilename(1:end-4) '_final_ACFseg.txt'],'a'); % adjust path if necessary!
    fid4=fopen([path2 '\' inputfilename(1:end-4) '_final_fitparameters_seg.txt'],'a'); % adjust path if necessary!
    fid5=fopen([path2 '\' inputfilename(1:end-4) '_final_curveincls.txt'],'a'); % adjust path if necessary!
    %fid3=fopen([path2 '\' inputfilename(1:end-4) '_files.out'],'a');#
    fid6=fopen([path2 '\' inputfilename(1:end-4) '_groupID.txt'],'a'); % adjust path if necessary!
    saveas(hh,[path2 '\' inputfilename(1:end-4) ' CCFs_seg.fig'])
    if errorfull==0
        saveas(h,[path2 '\' inputfilename(1:end-4) ' CCFs_full.fig'])
    end
    saveas(hI,[path2 '\' inputfilename(1:end-4) ' Is.fig'])
    writecell(allKeys,[path2 '\' inputfilename(1:end-4) ' _Czi_Metafile.csv'],'delimiter',',')
    
    outputseg=zeros(length(tcorr1segfinal),9);    
    outputseg(:,1)=tcorr1segfinal;
    outputseg(:,2)=fcorr1segfinal;
    outputseg(:,3)=weights1segfinal;
    outputseg(:,4)=tcorr2segfinal;
    outputseg(:,5)=fcorr2segfinal;
    outputseg(:,6)=weights2segfinal;
    if PIE==0
        outputseg(:,7)=tcorrccsegfinal;
        outputseg(:,8)=fcrosscorrsegfinal;
        outputseg(:,9)=weightsccsegfinal;
    else
%         output(:,7)=tcorrccfinal;
%         output(:,8)=fcrosscorrfinal;
%         output(:,9)=weightsccfinal;
        outputseg(:,7)=tcorrccsegfinal(1:end-1);
        outputseg(:,8)=fcrosscorrsegfinal(1:end-1);
        outputseg(:,9)=weightsccsegfinal(1:end-1);
    end
    
    if errorfull==0
        output=zeros(length(tcorr2),9);    
        output(:,1)=tcorr1final;
        output(:,2)=fcorr1final;
        output(:,3)=weights1final;
        output(:,4)=tcorr2final;
        output(:,5)=fcorr2final;
        output(:,6)=weights2final;
        if PIE==0
            output(:,7)=tcorrccfinal;
            output(:,8)=fcrosscorrfinal;
            output(:,9)=weightsccfinal;
        else
    %         output(:,7)=tcorrccfinal;
    %         output(:,8)=fcrosscorrfinal;
    %         output(:,9)=weightsccfinal;
            output(:,7)=tcorrccfinal(1:end-1);
            output(:,8)=fcrosscorrfinal(1:end-1);
            output(:,9)=weightsccfinal(1:end-1);
        end
    end
    
    if errorfull==0
        outputparameters=zeros(28,3);
        outputparameters(1:3,1)=[N1final;taud1final;Sfit1final];
        outputparameters(1:3,2:end)=CI1final;
        outputparameters(4:6,1)=[N2final;taud2final;Sfit2final];
        outputparameters(4:6,2:end)=CI2final;
        outputparameters(7:9,1)=[Nccfinal;tauccfinal;Sfitccfinal];
        outputparameters(7:9,2:end)=CIccfinal(1:3,:);
        outputparameters(10,1)=bleachingfraction1;
        outputparameters(11,1)=bleachingfraction2;
        outputparameters(12,1)=nanmean(line1fluorescenceseries);
        outputparameters(13,1)=nanmean(line2fluorescenceseries);
        outputparameters(14,1)=lsp_Ch1;
        outputparameters(15,1)=lsp_Ch2;
        outputparameters(16,1)=waist_Ch1;
        outputparameters(17,1)=waist_Ch2;
        outputparameters(18,1)=waist1_Ch1;
        outputparameters(19,1)=waist1_Ch2;
        outputparameters(20,1)=membranewidth_Ch1;
        outputparameters(21,1)=membranewidth_Ch2;
        outputparameters(22,1)=pixelsize;
        outputparameters(23,1)=pixeltime;
        outputparameters(24,1)=membranetime_Ch1;
        outputparameters(25,1)=membranetime_Ch2;
        outputparameters(26,1)=scantime_lines;
        outputparameters(27,1)=backgroundcorrectionCh1;
        outputparameters(28,1)=backgroundcorrectionCh2;
        
    end
    
    outputparametersseg=zeros(28,3);
    outputparametersseg(1:3,1)=[N1segfinal;taud1segfinal;Sfit1segfinal];
    outputparametersseg(1:3,2:end)=CI1segfinal;
    outputparametersseg(4:6,1)=[N2segfinal;taud2segfinal;Sfit2segfinal];
    outputparametersseg(4:6,2:end)=CI2segfinal;
    outputparametersseg(7:9,1)=[Nccsegfinal;tauccsegfinal;Sfitccsegfinal];
    outputparametersseg(7:9,2:end)=CIccsegfinal(1:3,:);
    outputparametersseg(10,1)=bleachingfraction1;
    outputparametersseg(11,1)=bleachingfraction2;
    outputparametersseg(12,1)=nanmean(line1fluorescenceseries);
    outputparametersseg(13,1)=nanmean(line2fluorescenceseries);
    outputparametersseg(14,1)=lsp_Ch1;
    outputparametersseg(15,1)=lsp_Ch2;
    outputparametersseg(16,1)=waist_Ch1;
    outputparametersseg(17,1)=waist_Ch2;
    outputparametersseg(18,1)=waist1_Ch1;
    outputparametersseg(19,1)=waist1_Ch2;
    outputparametersseg(20,1)=membranewidth_Ch1;
    outputparametersseg(21,1)=membranewidth_Ch2;
    outputparametersseg(22,1)=pixelsize;
    outputparametersseg(23,1)=pixeltime;
    outputparametersseg(24,1)=membranetime_Ch1;
    outputparametersseg(25,1)=membranetime_Ch2;
    outputparametersseg(26,1)=scantime_lines;
    outputparametersseg(27,1)=backgroundcorrectionCh1;
    outputparametersseg(28,1)=backgroundcorrectionCh2;
        
    if errorfull==0
        fprintf(fid1,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output');
        fprintf(fid2,'%e\t %e\t %e\n',outputparameters');
    end
    fprintf(fid3,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',outputseg');
    fprintf(fid4,'%e\t %e\t %e\n',outputparametersseg');
    
    curveincls=[curveincl1 curveincl2];
    fprintf(fid5,'%e\t %e\n',curveincls');
    fprintf(fid6,'%s\n',groupID);
    
    fid7=fopen([path2 '\' inputfilename(1:end-4) '_rel_cc_seg.txt'],'a'); % adjust path if necessary!
    outputcc=[relCCsegnew;...
        1/max([mean(fcorr1segfinal(1:5))/mean(fcrosscorrsegfinal(1:5)) mean(fcorr2segfinal(1:5))/mean(fcrosscorrsegfinal(1:5))]);...
        1/max([mean(fcorr1segfinal(1:10))/mean(fcrosscorrsegfinal(1:10)) mean(fcorr2segfinal(1:10))/mean(fcrosscorrsegfinal(1:10))])];
    fprintf(fid7,'%e\t %e\t %e\n',outputcc');
end
 
% path2=uigetdir;
%fid1=fopen([path2 '\' inputfilename(1:end-4) '_final_ACF.txt'],'a'); % adjust path if necessary!
% fid2=fopen([path2 '\' inputfilename(1:end-4) '_final_fitparameters.txt'],'a'); % adjust path if necessary!
%fid3=fopen([path2 '\' inputfilename(1:end-4) '_files.out'],'a');
% output=zeros(length(tcorr2),3);    
% 
% output(:,1)=tcorr1final;
% output(:,2)=fcorr1final;
% output(:,3)=weights1final;
% output(:,4)=fcorr2final;
% output(:,5)=weight2final;
% output(:,6)=fcorsscorrfinal;
% output(:,7)=weightsccfinal;

% outputparameters=zeros(4,3);
% outputparameters(:,1)=[Nfinal;taudfinal;Sfitfinal;bleachingfraction];
% outputparameters(1:3,2:end)=CIfinal;


%fprintf(fid1,'%e\t %e\t %e\n',output');
% fprintf(fid2,'%e\t %e\t %e\n',outputparameters');
