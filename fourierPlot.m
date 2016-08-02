function [] = fourierPlot(PSI_k, t_k, num, mult, handle)
% FourierPlot: Plots fourier modes of an AB/RW
% INPUT:
%       PSI_k: Fourior transform of the spatiotemporal WF.
%       t_k: time array
%       num: number of modes to plot.
%       mult: when used with talbot carpets. Just keep at 1.

% Set axes and extract data
cla(handle)
axes(handle);
data = PSI_k(:, 1:mult:end);          % Our data

% Plot
CO = lines(num);
LSO = {'-', '-.', '--'};

hData = cell(1,num);
legendInfo = cell(1,num);
for i = 1:num
    hData{i} = plot(t_k, data(:, i), 'color', CO(i, :), 'LineStyle', LSO{floor(i/8)+1}); 
    hold all;
end

for i=1:num
    legendInfo{i} = sprintf('A_{%d}', mult*(i-1));
end
hLegend = legend(legendInfo);
set(hLegend, 'Location', 'BestOutside');

hXLabel = xlabel('t');
hYLabel = ylabel('ln(|A_k|)');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel         ], ...
    'FontName'   , 'Helvetica');

set([hXLabel, hYLabel       ]  , ...
    'FontSize'   , 10          );

set(gca, ...
  'Box'         , 'off'        , ...
  'TickDir'     , 'out'        , ...
  'TickLength'  , [.02 .02]    , ...
  'XMinorTick'  , 'off'         , ...
  'YMinorTick'  , 'off'         , ...
  'XGrid'       , 'on'         , ...
  'YGrid'       , 'on'         , ...
  'XColor'      , [.3 .3 .3]   , ...
  'YColor'      , [.3 .3 .3]   , ...
  'LineWidth'   , 1            , ...
  'YLim'        , [-10, 2]    ,...
  'XLim'        , [min(t_k), ceil(max(t_k))]);

end