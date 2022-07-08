
function [out] = op_combine_water_subspecs(in,average)
% Some formats end up having subspectra in their reference scans
% (e.g. Philips), as well as empty lines. Intercept these cases
% here.
if nargin < 2
    average = 1;
end
    out = in;
    if out.subspecs == 2
        raw_A               = op_takesubspec(out,1);
        [raw_A]             = op_rmempty(raw_A);            % Remove empty lines
        raw_B               = op_takesubspec(out,2);
        [raw_B]             = op_rmempty(raw_B);            % Remove empty lines
        out                 = op_concatAverages(raw_A,raw_B);
    end

    if out.subspecs == 4
        raw_A               = op_takesubspec(out,1);
        [raw_A]             = op_rmempty(raw_A);            % Remove empty lines
        raw_B               = op_takesubspec(out,2);
        [raw_B]             = op_rmempty(raw_B);            % Remove empty lines
        raw_C               = op_takesubspec(out,3);
        [raw_C]             = op_rmempty(raw_C);            % Remove empty lines
        raw_D               = op_takesubspec(out,4);
        [raw_D]             = op_rmempty(raw_D);            % Remove empty lines
        out                 = op_concatAverages(op_concatAverages(raw_A,raw_B),op_concatAverages(raw_C,raw_D));
    end
    if out.subspecs == 1
        [out]             = op_rmempty(out);            % Remove empty lines
    end
    if average
        % Align and verage the refernce data
        if out.averages > 1 && ~out.flags.averaged
            [out]             = op_rmempty(out);
            [out,~,~]               = op_alignAverages(out, 1, 'n');
            out                     = op_averaging(out);            % Average
        else
            out.flags.averaged  = 1;
            out.dims.averages   = 0;
        end
    end
    out.names = {'A'};
end