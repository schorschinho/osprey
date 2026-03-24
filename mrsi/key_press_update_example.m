function key_press_update_example(spec)
    % Create a figure and set the key press function

    
    op_plotspec(spec);
    set(gcf,'KeyPressFcn', @keyPressCallback)

    spec.phase1 = 0;
    spec.phase0 = 0;
    % Store the data in the figure's app data memory
    setappdata(gcf, 'data', spec);

    % Callback function to handle key press
    function keyPressCallback(~, event)
        % Retrieve stored data
        data = getappdata(gcf, 'data');

        switch event.Key
            case 'uparrow'  % If 'up arrow' is pressed, increase the counter
                data=op_addphase(data ,0,0.00001);

            case 'downarrow' % If 'down arrow' is pressed, decrease the counter
                data=op_addphase(data ,0,0.00001);

            case 'leftarrow' % If space is pressed, reset the counter
                data.counter = 0;

            case 'rightarrow' % If 'Esc' is pressed, close the figure
                delete(fig);
                return;

            otherwise
                % Do nothing on other key presses
                return;
        end

        % Update the displayed counter value
        set(data.textHandle, 'String', sprintf('Counter: %d', data.counter));

        % Save updated data back to the figure
        setappdata(fig, 'data', data);
    end
end