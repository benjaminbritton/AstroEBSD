function [ error1,num_peak_out,error2] = PC_GAOpt( PC,Peak_Centre_ok,ScreenSize,Crystal_LUT,Crystal_UCell,phase_num)
%PC_GAOpt Hunt for the pattern centre - used for the genetic algorithm
nhat_gnom = EBSP_NormConv( Peak_Centre_ok,ScreenSize,PC);
rotdata = EBSP_Index(nhat_gnom,Crystal_LUT{phase_num},1/Crystal_UCell{phase_num}.efac,Crystal_UCell{phase_num});

num_peak_in=size(Peak_Centre_ok,1);
num_peak_out=rotdata.maxok;

peak_weight=num_peak_in-num_peak_out+1;
error1=peak_weight*rotdata.error;
error2=rotdata.error;
end
