function stop = optimplotresidual(x,optimValues,state)
    stop = false;
    switch state
        case 'iter'
            ax =gca;
            shift = max(abs(optimValues.residual)) + abs(ax.YLim(1));
            plot(optimValues.residual - shift);
              % Make updates to plot or guis as needed
        case 'interrupt'
              % Probably no action here. Check conditions to see  
              % whether optimization should quit.
        case 'init'
              hObject = figure;
              plot(optimValues.residual);
              hold on
              shift = 0;
              % Setup for plots or guis
        case 'done'
              % Cleanup of plots, guis, or final plot
        otherwise
    end
    
end

