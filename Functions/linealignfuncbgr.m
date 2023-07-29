function [linefluorescenceseries,linearrayalignedba] = linealignfuncbgr(mask,linearraymasked,backgroundmasked,backgleft)
global  blocksize spatialfilter membranewidth waist waist1
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
backgroundblocks=membranepositions;
gsigmoid=fittype(@(a1,b1,c1,a2,c2,d1,x) a1*exp(-((x-b1)./c1).^2)+c2./(1+exp(-a2*(x-b1)))+d1);
nursigmoid=@(b1,a2,c2,d1,x) c2./(1+exp(-a2*(x-b1)))+d1;
% backgroundmasked=backgroundmasked(~isnan(backgroundmasked));
coeff_b1=zeros(length(membranepositions),1);
coeff_a2=zeros(length(membranepositions),1);
coeff_c2=zeros(length(membranepositions),1);
coeff_d1=zeros(length(membranepositions),1);
waist=zeros(length(membranepositions),1);

for i=1:length(membranepositions)
    linearrayblock=mean(linearraymasked((i-1)*blocksize+1:i*blocksize,:),1);
    
    backgroundblock=nanmean(backgroundmasked((i-1)*blocksize+1:i*blocksize,:));
    backgroundblock(backgroundblock==0)=NaN;
    backgroundblock=backgroundblock(~isnan(backgroundblock));
    backgroundblock=nanmean(backgroundblock);
    if isnan(backgroundblock)
        backgroundblock=0;
    end
    [~,membranepositionblock]=max(linearrayblock);
%     membranepositionblock
    membranepositions(i)=membranepositionblock;
    linearrayblock(linearrayblock==0)=NaN;
    linearrayblock=linearrayblock(~isnan(linearrayblock));
    [membmax,membranepositionblock2]=max(linearrayblock);
   % membfittemp=fit((1:1:length(linearrayblock))',linearrayblock',gsigmoid, 'StartPoint',[membmax-0.5*backgroundblock membranepositionblock2 2 backgleft*-2 backgroundblock 0.01 ]);
membfittemp=fit((1:1:length(linearrayblock))',linearrayblock',gsigmoid, 'StartPoint',[membmax-0.5*backgroundblock membranepositionblock2 2 backgleft*-2 ...
        backgroundblock 0.01 ], 'Lower', [0 0 0 0 -10^3 0 ] );
  if i==13
    figure(13)
    plot(linearrayblock, 'b.')
    hold on
    plot(membfittemp,'k-')
    plot(nursigmoid(membfittemp.b1, membfittemp.a2,membfittemp.c2,membfittemp.d1,(1:1:length(linearrayblock))'),'r-')
    plot(gsigmoid(membfittemp.a1, membfittemp.b1,membfittemp.c1, 1, 0, 0,(1:1:length(linearrayblock))'),'b-')
    % plot(gsigmoid(membmax-0.5*backgroundblock, membranepositionblock2, 2, backgleft*-2, backgroundblock, 0.01 ,(1:1:length(linearrayblock))'),'g-')
    xlabel("pixelposition [a.u.]")
    ylabel("Intensity")
    ytickformat('%.2f')
    pause
    hold off
    close (13)
  end

  waist(i)=membfittemp.c1*sqrt(2);
  coeff_b1(i)=membfittemp.b1;
  coeff_a2(i)=membfittemp.a2;
  coeff_c2(i)=membfittemp.c2;
  coeff_d1(i)=membfittemp.d1;
end


figure('Name','Waist Histogram')
try
histfit(waist, 60,'burr')
output=fitdist(waist,'burr');
xlabel("waist")
ylabel("frequency")
xtickformat('%.1f')
pause
waist1=output.alpha*(((output.c-1)/(output.c*output.k+1))^(1/output.c));
end


% Alignment of scan lines (maximum fluorescence), including background
% correction
fprintf('Aligning lines...\n');
maxdriftba=max(membranepositions)-min(membranepositions);
linearrayalignedba=zeros(length(linearraymasked(:,1)),length(linearraymasked(1,:))+maxdriftba);
firstblockindex=membranepositions(1);
maxdeviation=max(membranepositions)-firstblockindex;

%figure(33)% blockdrifts=zeros(size(membranepositions));

for jj=1:length(membranepositions)
    blockdrift=membranepositions(jj)-firstblockindex;
    shiftindex=maxdeviation-blockdrift;
%   linearrayalignedba((j-1)*blocksize+1:j*blocksize,1+shiftindex:shiftindex+length(linearraymasked(1,:)))=linearraymasked((j-1)*blocksize+1:j*blocksize,:);   
%   plot(mean(linearraymasked((jj-1)*blocksize+1:jj*blocksize,:),1),'r-')
%   hold on
%   plot(mean(mask((jj-1)*blocksize+1:jj*blocksize,:),1).* ...
%   (nursigmoid(membranepositions(jj), coeff_a2(jj),coeff_c2(jj),coeff_d1(jj),(1:1:size(linearraymasked,2))'))', 'k-')
    linearrayalignedba((jj-1)*blocksize+1:jj*blocksize,1+shiftindex:shiftindex+length(linearraymasked(1,:)))= ...
    linearraymasked((jj-1)*blocksize+1:jj*blocksize,:)-mean(mask((jj-1)*blocksize+1:jj*blocksize,:),1).* ...
          (nursigmoid(membranepositions(jj), coeff_a2(jj),coeff_c2(jj),coeff_d1(jj),(1:1:size(linearraymasked,2))'))';
%   plot(mean(linearrayalignedba((jj-1)*blocksize+1:jj*blocksize,1+shiftindex:shiftindex+length(linearraymasked(1,:))),1),'g-')
%   hold off
%   pause
%   
end
% figure(100)
% plot(backgroundblocks)

% Sum of all lines and fit
linearrayalignedmasum=sum(linearrayalignedba,1);

% Try out: Gauss function plus heavyside function
% ustep = @(x,range) range>=x;
% g=fittype(@(a1,a2,a3,b1,c1,x) a1*exp(-((x-b1)./c1).^2)+a2*ustep(b1,x)+a3);

% coeffnames(gsigmoid);
% membfit=fit([1:1:length(linearrayalignedmasum)]',linearrayalignedmasum'./sum(linearrayalignedmasum),gsigmoid,'StartPoint',[1 100 1 2 0.05 0.01]);
% membfit=fit([1:1:length(linearrayalignedmasum)]',linearrayalignedmasum'./sum(linearrayalignedmasum),g,'StartPoint',[1 0.1 0 64 1]);
membfit=fit((1:1:length(linearrayalignedmasum))',linearrayalignedmasum'./sum(linearrayalignedmasum),'gauss1');
mu=membfit.b1;
sigma=membfit.c1/sqrt(2);
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
xlabel('pixelposition')
ylabel('intensity')
ytickformat('%.1f')
sum(linearrayalignedmasum)
fitprofile=membfit.a1 .*gaussmf(1:length(linearrayalignedmasum), [sigma membfit.b1 ]);
profile=linearrayalignedmasum./sum(linearrayalignedmasum);
%%waist1=2*sigma;
save('fitprofile.mat','fitprofile');
save('profile.mat','profile');

% filtering of membrane region
linearrayalignedba(:,1:indexcutlow)=NaN;
linearrayalignedba(:,indexcutup:end)=NaN;

% Calculation of time series
linefluorescenceseries=nansum(linearrayalignedba,2);