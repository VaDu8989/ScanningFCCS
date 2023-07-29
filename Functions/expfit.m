function [correctionfit] = expfit(linefluorescenceseries,timeline,timelinefull)
global ff
ff=fit(timeline,linefluorescenceseries,'exp2');
correctionfit=ff.a.*exp(ff.b*timelinefull)+ff.c.*exp(ff.d*timelinefull);
end