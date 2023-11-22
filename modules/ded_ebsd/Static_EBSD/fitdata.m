function [fitresult,gof] = fitdata(xdata,ydata)
[xData, yData] = prepareCurveData( xdata, ydata );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf 0];
opts.Normalize = 'off';
opts.Robust = 'LAR';
opts.Upper = [Inf 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end

