function [tcorr,fcorr,sigmas] = crosscorrFCSmultipletauPIE(timeseries1,timeseries2)
% Calculation of autocorrelation function of measured fluorescence
% intensity using multiple tau algorithm
% Evt.: Vern�nftiges maxtau �bergeben, bzw. so lassen und GUI mit Slider in
% Analyse einf�hren
meansignal1=nanmean(timeseries1);
meansignal2=nanmean(timeseries2);
deltatimeseries1=timeseries1-meansignal1;
deltatimeseries2=timeseries2-meansignal2;
maxtau=length(timeseries1)-1;
fcorr=zeros(1,length(maxtau)+1); % noch falsch, aber funktioniert vorerst...
tcorr=zeros(1,length(maxtau)+1);
sigmas=zeros(1,length(maxtau)+1);
m=16;
for i=0:m-1 
%for i=1:m-1    %changed for 2c2f
    fcorri=0;
    sigmaisq=0;
    counter=0;
        for ii=1:maxtau-i+1
            if isnan(deltatimeseries1(ii))==0 && isnan(deltatimeseries2(ii+i))==0
                fcorri=fcorri+deltatimeseries1(ii)*deltatimeseries2(ii+i);
                sigmaisq=sigmaisq+deltatimeseries1(ii)^2*deltatimeseries2(ii+i)^2;
                counter=counter+1;
            end
        end
        sigmai=1/(sqrt(counter)*meansignal1*meansignal2)*sqrt(1/(counter)*(sigmaisq-1/(counter)*fcorri^2));
        sigmas(i+1)=sigmai; 
        fcorr(i+1)=fcorri/(counter);
        tcorr(i+1)=i;
%         sigmas(i)=sigmai;   %changed for 2c2f
%         fcorr(i)=fcorri/(counter);
%         tcorr(i)=i;
end
tlag=1;
i=m+1; 
%i=m;   %changed for 2c2f
iterations=1;
while length(timeseries1)>m
    if mod(length(timeseries1),2)==1
        timeseries1=timeseries1(1:end-1);
        timeseries2=timeseries2(1:end-1);
    end
    timeseries1=nanmean(reshape(timeseries1,2,length(timeseries1)/2),1); % bin every 2 intensities
    timeseries2=nanmean(reshape(timeseries2,2,length(timeseries2)/2),1); % bin every 2 intensities
    deltatimeseries1=timeseries1-meansignal1;
    deltatimeseries2=timeseries2-meansignal2;
    maxtau=length(timeseries1)-1;
    tlag=tlag*2;
    iterations=iterations+1;
    
    if maxtau<m
        maxindex=maxtau;
    else
        maxindex=m-1;
    end
    for j=m/2:1:maxindex
        fcorri=0;
        sigmaisq=0;
        counter=0;
        for ii=1:maxtau-j+1
            if isnan(deltatimeseries1(ii))==0 && isnan(deltatimeseries2(ii+j))==0
                fcorri=fcorri+deltatimeseries1(ii)*deltatimeseries2(ii+j);
                sigmaisq=sigmaisq+deltatimeseries1(ii)^2*deltatimeseries2(ii+j)^2;
                counter=counter+1;
            end
        end
        sigmai=1/(sqrt(counter)*meansignal1*meansignal2)*sqrt(1/(counter)*(sigmaisq-1/(counter)*fcorri^2));
        sigmas(i)=sigmai;
        fcorr(i)=fcorri/(counter);
        tcorr(i)=j*tlag;
        i=i+1;
    end
end
fcorr=fcorr./(meansignal1*meansignal2);
%sigmas(end)=sigmas(end-1);
sigmas(end-m/2+1:end)=sigmas(end-m/2);


% tcorr=tcorr(1:end-1);
% fcorr=fcorr(1:end-1);
% sigmas=sigmas(1:end-1);
%sigmas(end-m/2+1:end)=sigmas(end-m/2); % last data point: only one summand, i.e. variance zero!