function [correctionfit1,correctionfit2] = expfit2(linefluorescenceseries1,timeline1,timelinefull1,linefluorescenceseries2,timeline2,timelinefull2)
global ff1 ff2
ff1=fit(timeline1,linefluorescenceseries1,'exp2' );



correctionfit1=ff1.a.*exp(ff1.b*timelinefull1)+ff1.c.*exp(ff1.d*timelinefull1);
ff2=fit(timeline2,linefluorescenceseries2,'exp2' );

correctionfit2=ff2.a.*exp(ff2.b*timelinefull1)+ff2.c.*exp(ff2.d*timelinefull2); %Korrekt mit timelinefull1 in erster Zeile?)
% figure(100)
% clf
% plot(ff1)
% hold on
% plot(ff2)
end
