function densityPlot(PSI, x, t, tIntr, xIntr, handle)
% densityPlot: used to plot flat density maps of RWs and ABs.
% Input: PSI: this is the matrix with the data points. Note that it can't
%             be complex, if you are interested in plotting intensity then 
%             input abs(PSI.^2) when calling the function.
%        x:   x-array
%        t:   y-axis array, usually time.
%        tIntr: how many points to sample from t-array.
%        xIntr: how many points to sample from x-array.
%        handle: handle to the axes where you want to plot it. If you are
%                using this function manually, you can use h=axes(); 
%                and pass h.

% Set current axes
axes(handle);

% Plot
surf(x(1:xIntr:end), t(1:tIntr:end), PSI(1:tIntr:end, 1:xIntr:end), ...
                                                'EdgeColor', 'none'); 

% Adjust misc plot components
colormap('Jet');
colorbar('eastoutside') 
ylabel('t'); xlabel('x'); zlabel('|\psi|^2'); 
view([0 0 90]) 

hXLabel = xlabel('x');
hYLabel = ylabel('t');

pk_max = max(max(PSI));
pk_min = min(min(PSI));
caxis([pk_min pk_max]);

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel         ], ...
    'FontName'   , 'Helvetica');
set([hXLabel, hYLabel       ]  , ...
    'FontSize'   , 10          );

set(gca, ...
  'Box'         , 'off'           , ...
  'TickDir'     , 'out'           , ...
  'TickLength'  , [.02 .02]       , ...
  'XMinorTick'  , 'off'           , ...
  'YMinorTick'  , 'off'           , ...
  'XGrid'       , 'off'           , ...
  'YGrid'       , 'off'           , ...
  'XColor'      , [.3 .3 .3]      , ...
  'YColor'      , [.3 .3 .3]      , ...
  'LineWidth'   , 1               , ...
  'XLim'        , [min(x) max(x)] ,...
  'YLim'        , [floor(min(t)) ceil(max(t))]);
end