function str_out = osp_capitalize_first(str_in)
%% str_out = osp_capitalize_first(str_in)
%   This function capitalizes the first letter of an input string while
%   leaving all remaining characters unchanged.
%
%   The function handles empty strings gracefully by returning them
%   unchanged. This is useful for formatting display labels, metabolite
%   names, or other text elements in the Osprey GUI.
%
%   USAGE:
%       str_out = osp_capitalize_first(str_in);
%
%   INPUTS:
%       str_in      = Input string or character array.
%
%   OUTPUTS:
%       str_out     = String with first letter capitalized.
%
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%
    if isempty(str_in)
        str_out = str_in;
    else
        str_out = [upper(str_in(1)), str_in(2:end)];
    end
end