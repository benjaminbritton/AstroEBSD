function [fitresult,gof] = fitdata_nonorm(xdata,ydata)
[xData, yData] = prepareCurveData( xdata, ydata );

% % Set up fittype and options.
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% % opts.Lower = [-Inf 0];
% opts.Normalize = 'off';
% opts.Robust = 'LAR';
% % opts.Upper = [Inf 0];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );

% Set up fittype and options.
ft = fittype( 'p1*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = 0;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end

