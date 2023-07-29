function [tcorr,fcorr,sigmas] = autocorrFCSmultipletau(timeseries)
% Calculation of autocorrelation function of measured fluorescence
% intensity using multiple tau algorithm
% Evt.: Vern�nftiges maxtau �bergeben, bzw. so lassen und GUI mit Slider in
% Analyse einf�hren
meansignal=nanmean(timeseries);
deltatimeseries=timeseries-meansignal;
maxtau=length(timeseries)-1;
fcorr=zeros(1,length(maxtau));
tcorr=zeros(1,length(maxtau));
sigmas=zeros(1,length(maxtau));
m=16;
for i=1:m-1
    fcorri=0;
    sigmaisq=0;
    counter=0;
        for ii=1:maxtau-i+1
            if (isnan(deltatimeseries(ii)+deltatimeseries(ii+i)))==0
                fcorri=fcorri+deltatimeseries(ii)*deltatimeseries(ii+i);
                sigmaisq=sigmaisq+deltatimeseries(ii)^2*deltatimeseries(ii+i)^2;
                counter=counter+1;
            end
        end
        sigmai=1/(sqrt(counter)*meansignal^2)*sqrt(1/(counter)*(sigmaisq-1/(counter)*fcorri^2));
        sigmas(i)=sigmai;
        fcorr(i)=fcorri/(counter);
        tcorr(i)=i;
end



tlag=1;
i=m;    
iterations=1;
while length(timeseries)>m
    if mod(length(timeseries),2)==1
        timeseries=timeseries(1:end-1);
    end
    timeseries=nanmean(reshape(timeseries,2,length(timeseries)/2),1); % bin every 2 intensities  &&&ok nanmean or nansum. Compare with Jonas
    deltatimeseries=timeseries-meansignal;
    maxtau=length(timeseries)-1;
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
            if (isnan(deltatimeseries(ii)+deltatimeseries(ii+j)))==0
                fcorri=fcorri+deltatimeseries(ii)*deltatimeseries(ii+j);
                sigmaisq=sigmaisq+deltatimeseries(ii)^2*deltatimeseries(ii+j)^2;
                counter=counter+1;
            end
        end
        sigmai=1/(sqrt(counter)*meansignal^2)*sqrt(1/(counter)*(sigmaisq-1/(counter)*fcorri^2));
        sigmas(i)=sigmai;
        fcorr(i)=fcorri/(counter);
        tcorr(i)=j*tlag;
        i=i+1;
    end
end
%i
%maxtau-j+1
fcorr=fcorr./meansignal^2;
%sigmas(end)=sigmas(end-1);
sigmas(end-m/2+1:end)=sigmas(end-m/2); % last data point: only one summand, i.e. variance zero!