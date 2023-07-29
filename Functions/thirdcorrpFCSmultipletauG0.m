function [tcorr1,tcorr2,fcorr,sigmas] = thirdcorrpFCSmultipletauG0(timeseries)
% Calculation of third order correlation function of measured fluorescence
% intensity using multiple tau algorithm
% Evt.: Vernünftiges maxtau übergeben, bzw. so lassen und GUI mit Slider in
% Analyse einführen
meansignal=nanmean(timeseries);
deltatimeseries=timeseries-meansignal;
maxtau=length(timeseries)-1;
fcorr=zeros(length(maxtau),length(maxtau));
tcorr1=zeros(1,length(maxtau));
tcorr2=zeros(1,length(maxtau));
sigmas=zeros(1,length(maxtau));
m=16;
%for i=1:m-1
for i=0:m-1
    for k=0:m-1
        fcorri=0;
        %sigmaisq=0;
            for ii=1:maxtau-i-k+1
                %if isnan(deltatimeseries(ii)+deltatimeseries(ii+i))==0
                    fcorri=fcorri+deltatimeseries(ii)*deltatimeseries(ii+i)*deltatimeseries(ii+k);
                    %sigmaisq=sigmaisq+deltatimeseries(ii)^2*deltatimeseries(ii+i)^2;
                %end
            end
            %sigmai=1/(sqrt(maxtau+1-i)*meansignal^2)*sqrt(1/(maxtau-i+1)*(sigmaisq-1/(maxtau-i+1)*fcorri^2));
            %sigmas(i+1)=sigmai;
            fcorr(i+1,k+1)=fcorri/(maxtau-i-k+1);
            tcorr2(k+1)=k;
    end
    tcorr1(i+1)=i;
end



tlag=1;
%i=m;
i=m+1;
k=m+1;
iterations=1;
while length(timeseries)>m
    if mod(length(timeseries),2)==1
        timeseries=timeseries(1:end-1);
    end
    timeseries=nanmean(reshape(timeseries,2,length(timeseries)/2),1); % bin every 2 intensities
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
        kstep=0;
        for l=m/2:1:maxindex
        fcorri=0;
        %sigmaisq=0;
        for ii=1:maxtau-j-l+1
            %if isnan(deltatimeseries(ii)+deltatimeseries(ii+j))==0
                fcorri=fcorri+deltatimeseries(ii)*deltatimeseries(ii+j)*deltatimeseries(ii+l);
                %sigmaisq=sigmaisq+deltatimeseries(ii)^2*deltatimeseries(ii+j)^2;
            %end
        end
        %sigmai=1/(sqrt(maxtau+1-j)*meansignal^2)*sqrt(1/(maxtau-j+1)*(sigmaisq-1/(maxtau-j+1)*fcorri^2));
        %sigmas(i)=sigmai;
        fcorr(i,k+kstep)=fcorri/(maxtau-j-l+1);
        tcorr2(k+kstep)=l*tlag;
        kstep=kstep+1;
        end
        tcorr1(i)=j*tlag;
        i=i+1;
    end
    k=i;
end
%i
%maxtau-j+1
fcorr=fcorr./meansignal^3;
%sigmas(end)=sigmas(end-1);
%sigmas(end-m/2+1:end)=sigmas(end-m/2); % last data point: only one summand, i.e. variance zero!