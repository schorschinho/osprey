function gui = addCellToGuiAsTable(gui, CellToAdd)
% Very involved procedure for making table display work for
% newer MATLAB versions. 

% Check which path we need to go
if isMATLABReleaseOlderThan("R2025a")
    % This is how we've integrated the Quantify tables in the past
    temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', CellToAdd,...
        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    set( temp, 'BackgroundColor',gui.colormap.Background);
    delete(gui.Plot.quant.Children)
    set( temp, 'Parent', gui.Plot.quant ); %Update table
else
    % For newer version, we have seen that the HTML doesn't render well
    % using uimulticollist (see
    % https://github.com/schorschinho/osprey/issues/838)

    % We can here create a workaround that is actually much nicer rendered
    % (using uitable), but it DOESN'T work in MATLAB 2025 and older (I
    % think).

    % Convert back to a table so we can use uitable in the GUI
    back2Table = array2table(CellToAdd);
    % Use first row as variable names, then delete first row
    back2Table.Properties.VariableNames = back2Table{1,:};
    back2Table(1,:) = [];
    delete(gui.Plot.quant.Children);
    temp = uitable(gui.Plot.quant, "Data", back2Table);
end
