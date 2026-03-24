function str_out = osp_capitalize_first(str_in)
    % Capitalize first letter of string
    if isempty(str_in)
        str_out = str_in;
    else
        str_out = [upper(str_in(1)), str_in(2:end)];
    end
end