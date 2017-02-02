function nlsePlot(PSI, x, t, tIntr, xIntr, handle, mode)
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
delete(handle.Children);

% Plot
switch lower(mode)
    case 'density'
        surf(x(1:xIntr:end), t(1:tIntr:end), PSI(1:tIntr:end, 1:xIntr:end), ...
                                                'EdgeColor', 'none'); 
        % Adjust misc plot components
        colormap('Jet'); colorbar('eastoutside'); view([0 0 90]);
    case '3d'
        h = surf(x(1:xIntr:end), t(1:tIntr:end), PSI(1:tIntr:end, 1:xIntr:end), 'EdgeColor', 'none'); 
        colorbar off;
        view(-62,42)
        shading interp
        lightangle(-40,50);
        h.FaceLighting = 'gouraud';
        h.AmbientStrength = 0.3;
        h.DiffuseStrength = 0.8;
        h.SpecularStrength = 0.5;
        h.SpecularExponent = 3;
        h.BackFaceLighting = 'reverselit';
        colormap(jet(256));
        grid off;
    otherwise
        error('Unknown plotting mode');
end

hXLabel = xlabel('x');
hYLabel = ylabel('t');
hZLabel = zlabel('|\psi|^2'); 

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel, hZLabel         ], ...
    'FontName'   , 'Helvetica');
set([hXLabel, hYLabel, hZLabel         ]  , ...
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
  'XLim'        , [-max(x) max(x)] ,...
  'YLim'        , [floor(min(t)) ceil(max(t))]);

hx = xlabel('$t$','FontSize',11,'Interpreter','latex');
zlabel('$|\psi|^2$','FontSize',11,'Interpreter','latex');
hy = ylabel('$x$','FontSize',11,'Interpreter','latex');

end