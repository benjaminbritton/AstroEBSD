function [ScaleDetails]=Plot_SetScaleBar(MapOut)
%pick a scale bar size

fs_ok=[0.005 0.01 0.05 0.1 0.2 0.5 1 5 10 20 50 100 200]; %allowable scale bare lengths in um
x_range=range(MapOut.X_axis);

fs_check=x_range/5; %find a fifth of the side length
[~,ix]=min(abs(fs_check-fs_ok)); %find the scale bar that it closest to this

ScaleDetails.Len=fs_ok(ix);

ScaleDetails.X=0.9*max(MapOut.X_axis)-[0 ScaleDetails.Len];
ScaleDetails.Y=[0 0]+0.9*(max(MapOut.Y_axis));
ScaleDetails.Step=round(MapOut.XStep,3);

if ScaleDetails.Len >= 0.5
    ScaleDetails.Step=round(ScaleDetails.Step,3);
    ScaleDetails.Text=['Scalebar = ' num2str(ScaleDetails.Len) ' \mum; Step = ' num2str(ScaleDetails.Step) ' \mum'];
else %convert to nm
    ScaleDetails.Step=round(ScaleDetails.Step,5);
    ScaleDetails.Text=['Scalebar = ' num2str(ScaleDetails.Len*1000) ' nm; Step = ' num2str(ScaleDetails.Step*1000) ' nm'];
end
 