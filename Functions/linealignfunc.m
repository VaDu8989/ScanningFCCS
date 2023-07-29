function [linefluorescenceseries,linearrayalignedba] = linealignfunc(linearraymasked)
global  blocksize spatialfilter membranewidth waist1  waist
% movingaveragewindowsize spatialfilter depletioncorrection...
% ...intensityfilter intensitythreshold

% % Moving average (to determine membrane positions)
% fprintf('Moving average...\n');
% membranepositions=zeros(length(linearraymasked(:,1))-movingaveragewindowsize+1,1);
% for i=1:length(linearraymasked(:,1))-movingaveragewindowsize+1
%     linearraypart=sum(linearray(i:i+movingaveragewindowsize-1,:),1);
%     [maxpart,membranepositionpart]=max(linearraypart);
%     membranepositions(i)=membranepositionpart;
% end
% 
% % Alignment of scan lines (maximum fluorescence)
% fprintf('Aligning lines...\n');
% maxdriftma=max(membranepositions)-min(membranepositions);
% linearrayalignedma=zeros(length(linearraymasked(:,1)),length(linearraymasked(1,:))+maxdriftma);
% 
% firstlineindex=membranepositions(1);
% maxdeviation=max(membranepositions)-firstlineindex;
% for j=movingaveragewindowsize/2:length(linearraymasked(:,1))-movingaveragewindowsize/2
%     linedrift=membranepositions(j-movingaveragewindowsize/2+1)-firstlineindex;
%     shiftindex=maxdeviation-linedrift;
%     if j==1
%         linearrayalignedma(1:movingaveragewindowsize/2,1+shiftindex:shiftindex+length(linearraymasked(1,:)))=linearraymasked(1:movingaveragewindowsize/2,:);
%     else
%         linearrayalignedma(j,1+shiftindex:shiftindex+length(linearraymasked(1,:)))=linearraymasked(j,:);
%     end
%     if j==length(linearraymasked(:,1))-movingaveragewindowsize/2
%         linearrayalignedma(j:end,1+shiftindex:shiftindex+length(linearraymasked(1,:)))=linearraymasked(j:end,:);
%     end
% end

% Block average
fprintf('Block average...\n');
membranepositions=zeros(length(linearraymasked(:,1))/blocksize,1);
waist=zeros(length(membranepositions),1);
for i=1:length(membranepositions)
    linearrayblock=sum(linearraymasked((i-1)*blocksize+1:i*blocksize,:),1);
    [~,membranepositionblock]=max(linearrayblock);
    membranepositions(i)=membranepositionblock;
    membfittemp=fit([1:1:length(linearrayblock)]',linearrayblock'./sum(linearrayblock),'gauss1');
waist(i)=membfittemp.c1*sqrt(2);
end
figure
try
    histfit(waist, 60,'burr')
output=fitdist(waist,'burr')
pause
waist1=output.alpha*(((output.c-1)/(output.c*output.k+1))^(1/output.c))
end

% Alignment of scan lines (maximum fluorescence)
fprintf('Aligning lines...\n');
maxdriftba=max(membranepositions)-min(membranepositions);
linearrayalignedba=zeros(length(linearraymasked(:,1)),length(linearraymasked(1,:))+maxdriftba);
firstblockindex=membranepositions(1);
maxdeviation=max(membranepositions)-firstblockindex;
% blockdrifts=zeros(size(membranepositions));
for j=1:length(membranepositions)
    blockdrift=membranepositions(j)-firstblockindex;
    shiftindex=maxdeviation-blockdrift;
    linearrayalignedba((j-1)*blocksize+1:j*blocksize,1+shiftindex:shiftindex+length(linearraymasked(1,:)))=linearraymasked((j-1)*blocksize+1:j*blocksize,:);
end


% Sum of all lines and fit
linearrayalignedmasum=sum(linearrayalignedba,1);
size(linearrayalignedmasum)
% Try out: Gauss function plus heavyside function
% ustep = @(x,range) range>=x;
% g=fittype(@(a1,a2,a3,b1,c1,x) a1*exp(-((x-b1)./c1).^2)+a2*ustep(b1,x)+a3);
% coeffnames(g);
%membfit=fit([1:1:length(linearrayalignedmasum)]',linearrayalignedmasum'./sum(linearrayalignedmasum),g,'StartPoint',[1 0.1 0 64 1]);
membfit=fit([1:1:length(linearrayalignedmasum)]',linearrayalignedmasum'./sum(linearrayalignedmasum),'gauss1');
%  membfit.a1
%  membfit.a2
%  membfit.a3
mu=membfit.b1;
sigma=membfit.c1/sqrt(2);
waist1=2*sigma
indexmean=round(membfit.b1);
indexcutlow=round(mu-spatialfilter*sigma);
indexcutup=round(mu+spatialfilter*sigma);
membranewidth=indexcutup-indexcutlow-1;
%membranefitfunc=normpdf(1:length(linearrayalignedmasum),mu,sigma)
scrsz=   get(0,'ScreenSize');
figure('OuterPosition',[2*scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Total Membrane fluorescence peak')
plot(1:length(linearrayalignedmasum),linearrayalignedmasum./sum(linearrayalignedmasum),'-b')
hold on
plot(membfit,'-r')
xlabel('Pixel Number')
ylabel('Value')
sum(linearrayalignedmasum)

fitprofile=membfit.a1*exp(-(([1:1:length(linearrayalignedmasum)]-membfit.b1)/membfit.c1).^2);
profile=linearrayalignedmasum./sum(linearrayalignedmasum);
save('fitprofile.mat','fitprofile');
save('profile.mat','profile');

% filtering of membrane region
linearrayalignedba(:,1:indexcutlow)=NaN;
linearrayalignedba(:,indexcutup:end)=NaN;

% Calculation of time series
linefluorescenceseries=nansum(linearrayalignedba,2);