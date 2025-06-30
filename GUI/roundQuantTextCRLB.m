function QuantTextOut = roundQuantTextCRLB(QuantTextIn)
% This small function takes an OspreyQuantify table that is ready to be printed
% It then rounds the CRLB values and adds a [%] unit tag to the column
% name.
% (but only if a CRLB column exists)
%
% Georg Oeltzschner, Johns Hopkins University 2025

QuantTextOut = QuantTextIn;

% Check CRLB column
findCRLBCol = strcmp(QuantTextIn, 'CRLB');
if any(findCRLBCol(:))
    % Determine the column number
    CRLBCol = find(sum(findCRLBCol));
    % Change column title to include [%]
    QuantTextOut{1,CRLBCol} = 'CRLB [%]';
    for m = 2:size(QuantTextIn,1)
        QuantTextOut{m,4} = round(QuantTextIn{m,4}, 1);
    end

else
    % do nothing

end

end