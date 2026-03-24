function out_str = ConvertRelativePath(in_str)
if strcmp(in_str(1:5),'which')                                          % Full path is not given we have to eval
    strdef = append('out_str = ', in_str, ';');
    pattern = "''"; % Match find double quotes
    replacement = "'"; % Replace the matched string with just a quote
    strdef = regexprep(strdef, pattern, replacement);
    eval(strdef);
else
    out_str = in_str;
end


if ~isfile(out_str) && ~isempty(out_str)                                % Intercept if file doesn't exist
    error('basis file %s does not exist.', out_str);
end

end