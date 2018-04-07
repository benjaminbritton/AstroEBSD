function FilterMap=Plot_FilterGen(MapOut,Filtervals)

%create a blank map
FilterMap=ones(size(MapOut.Err));

%threshold on MAE
FilterMap(MapOut.Err> Filtervals.MAE_Thresh) = 0;

%threshol on IQ
FilterMap(MapOut.IQ < Filtervals.IQ_Thresh) = 0;

