%  A UIMULTICOLLIST must be created before it can be operated on
%     The initial values are arg pairs as per a UICONTROL of style listbox.
%
%  Calling with no arguments provides a dynamic list of the methods available:
%    UIMULTICOLLIST
%   
%
%    h = UIMULTICOLLIST ( figHandle );
%    h = UIMULTICOLLIST ( 'figHandle, units', 'normalized', 'position', [0 0 1 1] );
%  If no fig handle supplier - gcf is used
%    h = UIMULTICOLLIST ( 'units', 'normalized', 'position', [0 0 1 1] );
% 
%  After a UIMULTICOLLIST is created the following methods can be applied:
% 
%  Add a new column to the listbox
%     UIMULTICOLLIST ( h, 'addCol', colItems, colIndex,  colColour )   
%     UIMULTICOLLIST ( h, 'addCol', colItems, colColour, colIndex )   
%                           colItems : { 'A' ... 'B' }  (cell)   | cell array of items 
%             optional: 
%                The order of these 2 optional arguments is not important.
%                           colIndex : nCols+1          (int)    | index where column to be inserted
%                           colColour: 'BLACK'          (char)   | colour of column (HTML format colour allowed)
%
%  Add a new tow to the listbox
%     UIMULTICOLLIST ( h, 'addRow', rowItems , rowIndex )
%                           rowItems : { 'A' ... 'B' }  (cell)   | cell array of items 
%             optional:
%                           rolIndex : nRows+1          (int)    | index where row to be inserted
%
%  Change a specific item in the listbox
%     UIMULTICOLLIST ( h, 'changeItem', newValue, rowIndex, colIndex )
%                           newValue : the new list value        | char, int or double of new item in table
%                           rowIndex : int                       | row index, positive int
%                           colIndex : int                       | column index, positive int
%
%  Change an entire row in the table
%     UIMULTICOLLIST ( h, 'changeRow', newItems, rowIndex )
%                           newValue : the new list value        | char, int or double of new item in table
%                           rowIndex : int                       | row index, positive int
%
%  Change the colour of a specific column
%     UIMULTICOLLIST ( h, 'columnColour', colIndex, 'Colour' )
%                           colIndex : int                       | column index, positive int
%                           colColour: 'BLACK'          (char)   | colour of column (HTML format colour allowed)
%
%  Delete an entire column from the table
%     UIMULTICOLLIST ( h, 'delCol', colIndex )
%                           colIndex : int                       | column index, int
%                                                                | 0 -> remove the last column
%                                                                | negative number is the same as end-INT
%
%  Delete an entire row from the table
%     UIMULTICOLLIST ( h, 'delRow', rowIndex )
%                           rowIndex : int                       | row index, positive int
%                                                                | 0 -> remove the last row
%                                                                | negative number is the same as end-INT
%
%  Number of cols in the table
%     UIMULTICOLLIST ( h, 'nCols' )
%
%  Number of cols in the table
%     UIMULTICOLLIST ( h, 'nRows' )
% 
%  Return the current selection in numbers (if char, NaN is returned in its place)
%     UIMULTICOLLIST ( h, 'selectionInNumbers, colIndex )
%             optional:
%                           colIndex : 1          (int)    | only return the selected columns (default = all)
%
%  Return the selected row(s) - but only the stated col
%     UIMULTICOLLIST ( h, 'selectedStrCol', colIndex )
%                           colIndex : 1          (int)    | only return the selected columns (default = all)
%
%  Return the selected strings (use this instead of get ( h, 'string' ) as this is a HTML formatted string
%     UIMULTICOLLIST ( h, 'selectedString' )
%
%  Change the column separator
%    UIMULTICOLLIST ( h, 'separator', newColSeparator )
%                           newColSeparator : ' '  (char)    | replace the default ' ' as the column separator
%
%  Set the 1st row to be a header - this influences the filtering below
%     UIMULTICOLLIST ( h, 'setRow1Header', flag )
%                           flag : boolean flag (logical)    | turn on/off internal flag that the 1st row is a header
%
%  Interact with the string cell array
%    Return the full string
%     str = UIMULTICOLLIST ( h, 'string' )
%    update the full string
%     UIMULTICOLLIST ( h, 'string', newCellArray )
%                           newCellArray : square cell array    | see examples
%
%  Interact with the table value  
%    return the selected rows
%      value = UIMULTICOLLIST ( h, 'value' )
%    update the selected rows
%      value = UIMULTICOLLIST ( h, 'value', selectedRows )
%                           selectedRows : int scalar/array   | selected rows
%
%  FILTERING
%  ---------
%  Special case being trialed in this version of UIMULTICOLLIST to apply filters to rows
%  All of the above row interactions work on the underlying cell table.   
%    e.g. ADDROW, CHANGEROW, CHANGEITEM, DELROW, NROWS
%  They DO NOT act on the filtered table.
%  Any action on SELECTIONINNUMBERS, SELECTEDSTRCOL, SELECTEDSTR, VALUE 
%     will act on the filtered table data
%  
%  A command window method to apply filters:
%     UIMULTICOLLIST ( h, 'addFilter', colIndex, filterValue )
%                           colIndex    : 1   (int)      | column to be filtered
%                           filterValue :  filterItem(s) | cell array of filter to be applied to this col
%
%  A UICONTEXTMENU can be added to allow for interactive filtering
%     UIMULTICOLLIST ( h, 'applyUIFilter', useRow1AsLabels )
%             optional:
%                           useRow1AsLabels : boolean flag  | whether to use the 1st row as a flag (ignored in filtering)
%     If a listbox already has a UICONTEXTMENU - the method will add a menu item to the users context menu.
%
%  Remove the interactive filter
%     UIMULTICOLLIST ( h, 'removeUIFilter' )
%
%  Reset all the filters
%     UIMULTICOLLIST ( h, 'resetFilter' )
%
%  Ask the code what filters are currently applied 
%     output = UIMULTICOLLIST ( h, 'filtersApplied' )
%
%  The user can use get and set on the handle as per normal - but it is advised 
%   to use the above methods for interacting with the string, as the users string
%   is converted to html for formatting reasons.
%
%   see also uicontrol, uicontextmenu, uimenu, get, set
%
% %% Example code:
%     clear str
%     d = dir;
%     for i=1:length(d)
%       str{i,1} = d(i).name;
%       str{i,2} = d(i).date;
%       str{i,3} = d(i).bytes;
%       str{i,4} = d(i).isdir;
%       [dummy,dummy,ext] = fileparts ( d(i).name );
%       str{i,5} = strcmp ( ext, '.m' );  
%     end
%     h=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', str, 'columnColour', { 'RED' 'RED' 'BLUE' 'BLACK' 'GREEN' } );
% 
%     %% now add a header
%     header = { 'FileName' 'Date' 'Bytes' 'isDir' 'isMFile' };
%     uimulticollist ( h, 'addRow', header )
%     %% Say we wanted to add it to the beginning - remove and add again
%     uimulticollist ( h, 'delRow', 0 ); % 0 deletes the last row
%     %uimulticollist ( h, 'delRow', -1 ); % -N deletes the Nth from the end
%     %uimulticollist ( h, 'delRow', [1:6] ); % delete multiple rows
%     uimulticollist ( h, 'addRow', header, 1 )
%     %% add an extra column
%     nItems = uimulticollist ( h, 'nRows' );
%     uimulticollist ( h, 'addCol', [0:nItems-1] )
%     %% delete the end column and add a new one at the start
%     uimulticollist ( h, 'delCol', 0 ); % 0 deletes the last col
%     uimulticollist ( h, 'addCol', [0:nItems-1], 1 )
%     %% can also specify the colour when adding
%     uimulticollist ( h, 'delCol', 1 ); % 0 deletes the last col
%     uimulticollist ( h, 'addCol', [0:nItems-1], 1, 'PURPLE' )
%     %% The order of the index and colour doesn't matter
%     uimulticollist ( h, 'delCol', 1 ); % 0 deletes the last col
%     uimulticollist ( h, 'addCol', [0:nItems-1], 'PURPLE', 1 )
%     %% change an item in the list, index (row, col)
%     uimulticollist ( h, 'changeItem', 'ORDER', 1, 1 )
%     %% change an entire row
%     uimulticollist ( h, 'changeRow', { NaN 'NAME', 'Data' 'bytes' 'isDir' 'isMFile' }, 1 )
%     %% change an column col
%     uimulticollist ( h, 'columnColour', 1, 'YELLOW' )
%     %% retrieve selected items by column
%     set ( h, 'value', 1 )
%     uimulticollist ( h, 'selectedStrCol', [2] )
%     set ( h, 'Max', 2, 'value', 1:3 );
%     uimulticollist ( h, 'selectedStrCol', [2 3] )
%     %% return the full rows of the data seleted
%     uimulticollist ( h, 'selectedString' )
%     %% change the default separator
%     uimulticollist ( h, 'separator', ' | ' )
%     %% you can get the whole string - do soem work and put it back
%     str = uimulticollist ( h, 'string' );
%     str{end,end} = 'XYZ';
%     uimulticollist ( h, 'string', str );
%     %% you can get the value & set the value
%     value = uimulticollist ( h, 'value' );
%     uimulticollist ( h, 'value', [2 5]);
%     %% convert selection to numbers (if char -> NaN)
%     data = uimulticollist ( h, 'selectionInNumbers');
%     data = uimulticollist ( h, 'selectionInNumbers', 3);
%     %% add a context menu to filter
%     uimulticollist ( h, 'applyUIFilter');
%     %% remove a filter
%     uimulticollist ( h, 'removeUIFilter');
%     %% user already had a context menu
%     uic = uicontextmenu;
%     uimenu ( uic, 'label', 'Users Items' );
%     set ( h, 'uicontextmenu', uic );
%     uimulticollist ( h, 'applyUIFilter', 1); % 3rd arg uses row 1 as labels
%     %% reset a filter form commandline (if things go wrong)
%     uimulticollist ( h, 'resetFilter');
%     uimulticollist ( h, 'removeUIFilter');
%     %% filtering can be applied from commandline
%     uimulticollist ( h, 'setRow1Header', 1);
%     uimulticollist ( h, 'addFilter', 1, [10 11]); % 3rd arg uses row 1 as labels
%     %% what filtes are applied
%     filters = uimulticollist ( h, 'filtersApplied' )
% 
%   TODO:
%   -----
% A number of methods are planned for the future:
% TODO: If user passed in numbers - return numbers via the selectedStrCol method.
% TODO: Sorting lists
% TODO: If user using header - keep the haeder at the top of the list when scrolling
%
%  Copyright Robert Cumming
%            rcumming@matpi.com
%            robertjcumming@yahoo.co.uk
%            www.matpi.com
%
%  The following websites were used in the development of the code
%   http://www.undocumentedmatlab.com
%   http://blogs.mathworks.com/loren/2010/10/14/how-many-digits-to-write/
%
%
function varargout = uimulticollist ( varargin )
  newList = 1;                                   
  pairStart = 1;
  parentProvided = 0;
  if nargin == 0
    ParseFileForCase ( which ( mfilename ) );
    return
  end
  % is arg 1 a handle
  if ishandle ( varargin{1} )
    % is the handle a lisbox or a figure/dialog/uipanel
    if strcmp ( get ( varargin{1}, 'type' ), 'uicontrol' ) && strcmp ( get ( varargin{1}, 'style' ), 'listbox' )
      h = varargin{1};
      parent = ancestor ( varargin{1}, 'figure' );
      newList = 0;
    else
      parent = varargin{1};
      pairStart = 2;
      parentProvided = 1;
    end
  end
  % is this preparting a new list for the first time
  if newList
    if parentProvided == 0; parent = gcf; end % if its a new list with no parent -> use gcf
    if mod(length ( varargin ) - pairStart + 1 ,2) % check for pairs of inputs
      error ( 'uimulticollist:IncorrectInput', 'User did not provide valid pairs on input' );
    end

    % has the user provided a any special pairs as input?
    [varargin cString]    = CheckForSpecialInputs ( 'string',       varargin(pairStart:end), ''  ); 
    [varargin colours]    = CheckForSpecialInputs ( 'columnColour', varargin(pairStart:end), NaN ); % 3rd input is the detault.
    [varargin separator]  = CheckForSpecialInputs ( 'separator',    varargin(pairStart:end), ' ' );
    [varargin background] = CheckForSpecialInputs ( 'separator',    varargin(pairStart:end), 'white' );
    
    % create the basic uilistbox
    h = uicontrol ( 'style', 'listbox', 'parent', parent, varargin{:} );
    % if user hasn't specified a background colour -> set the default to white
%     if max ( userBackground ) == 0  % override the default (if the user didn't provide it)
      set ( h, 'background', background );
%     end
    s = get ( h, 'ApplicationData' );
    % make sure the strings are cells of cells
    if ~iscell ( cString ) % user didn't provide, or only one item,
      s.string{1,1} = cString;
    else
      s.string = EnsureAllStringInputCorrect ( cString );
    end
    % set up the other items
    s.colours = colours;
    s.separator = separator;
    set ( h, 'ApplicationData', s )
    varargout{1} = h;
  else
    % user interacting with an existing list
    iseven = mod(length ( varargin ),2)==0;
    if nargin > 1
      s = get ( h, 'ApplicationData' );
      [nRows,nCols]= size(s.string);
      switch varargin{2}
        case 'addCol'             % ( h, 'addCol', colItems , ColIndex )
          newData = CheckInputIsCell ( varargin{3} );
          if length(newData)~=nRows
            error ( 'uimulticollist:addCol', 'A new col must have an item in every row, input %i items, required %i', length(newData), nRows )
          end
          newColour = 'BLACK';
          userSpecifedLocation=0;
          if nargin >= 4
            if ischar ( varargin{4} ); newColour = varargin{4}; else userSpecifedLocation = 1; insertIndex = varargin{4}; end
          end
          if nargin == 5
            if ischar ( varargin{5} ); newColour = varargin{5}; else userSpecifedLocation = 1; insertIndex = varargin{5}; end
          end
          s = CheckColourVariable ( s, nCols );
          if userSpecifedLocation % user specified a location
            if insertIndex < 1 || insertIndex > nCols
              error ( 'uimulticollist:addCol', 'user specified a column which is outwith the bounds of the data, input %i, data %i', insertIndex, nCols );
            end
            if insertIndex > 1
              beforeIndex = s.string(:,1:insertIndex-1);
              s.string = [beforeIndex newData s.string(:,insertIndex:end)];
              s.colours = [s.colours(1:insertIndex-1) newColour s.colours(insertIndex:end)];
              if isfield ( s, 'filters' )
                s.filters = [s.filters(:,1:insertIndex-1) ones(nRows,1) s.filters(:,insertIndex:end)];
                s.filterDetails = [s.filterDetails(1:insertIndex) {[]}, s.filterDetails(insertIndex:end) ];
              end
            else
              s.string = [newData s.string];
              s.colours = [newColour s.colours];
              if isfield ( s, 'filters' )
                s.filters = [ones(nRows,1) s.filters];
                s.filterDetails = [ {[]} s.filterDetails];
              end
            end
          else
            s.string = [s.string newData];
            s.colours{end+1} = newColour;
            if isfield ( s, 'filters' )
              s.filters(:,end+1) = 1;
              s.filterDetails{end+1} = {};
            end
          end
          s.string = EnsureAllStringInputCorrect ( s.string );
        case 'addRow'             % ( h, 'addRow', rowItems , rowIndwex )
          newData = varargin{3};
          if ~iscell ( newData )
            error ( 'uimulticollist:addRow', 'A new row must be cell array of %i items', nCols )
          end
          if length(newData)~=nCols
            error ( 'uimulticollist:addRow', 'A new row must have an item in every col, input %i items, required %i', length(newData), nCols )
          end
          if nargin == 4 % user specified a location
            insertIndex = varargin{4};
            if insertIndex > 1
              beforeIndex = s.string(1:insertIndex-1,:);
              s.string = [beforeIndex; newData; s.string(insertIndex:end,:)];
            else
              s.string = [newData; s.string];
            end
          else
            s.string = [s.string; newData];
          end
          s.string = EnsureAllStringInputCorrect ( s.string );
        case 'addFilter'          % ( h, 'addFilter', colIndex, filterValue )
          if nargin ~= 4 
            error ( 'uimulticollist:addFilter', 'addFilter requires 4 input arguments' );
          end
          if varargin{3} > nCols
            error ( 'uimulticollist:addFilter', 'User provided col index, %i, which is greater than number cols, %i', varargin{3}, nCols );
          end
          FilterTableData ( varargin{4}, NaN, h, varargin{3} )
          return
        case 'applyUIFilter'      % ( h, 'applyUIFilter', useRow1AsLabels )
          s.uic = get ( h, 'uicontextmenu' );
          if isempty ( s.uic ) % if no menu -> create one
            s.uic = uicontextmenu ( 'parent', parent' );
          end
          names = get ( get ( s.uic, 'children' ), 'tag' );
          if ~isfield ( s, 'filterMenuName' ) % check to see if it already has a filter.
            s.filterMenuName = 'uimulticollist_uic';
            if max ( strcmp ( names, s.filterMenuName ) ); % does the name already exist in the menu
              count = 0;                                   % loop until a unique name is found (needed for deleting later)
              while max(max ( strcmp ( names, s.filterMenuName ) ))
                s.filterMenuName = sprintf ( '%s%i',s.filterMenuName, count );
                count = count + 1;
              end
            end
            uimenu ( s.uic, 'label', 'Filter', 'tag', s.filterMenuName );
            set ( h, 'uicontextmenu', s.uic );
            s.uicRow1Labels = nargin ==3;
          end
        case 'changeItem'         % ( h, 'changeItem, 'NewValue', rowIndex, colIndex )
          if nargin ~= 5
            error ( 'uimulticollist:changeItem', 'Incorrect number of input arguments - required 5' )
          end
          newItem = varargin{3};
          row = varargin{4};
          col = varargin{5};
          if row < 0 || row > nRows
            error ( 'uimulticollist:changeItem', 'Row index not valid, input %i items, required %i', row, length(s.string) )
          end
          if col < 0 || row > nCols
            error ( 'uimulticollist:changeItem', 'Row index not valid, input %i items, required %i', col, length(s.string{1}) )
          end
          s.string{row,col} = newItem;
          s.string = EnsureAllStringInputCorrect ( s.string );
        case 'changeRow'          % ( h, 'changeRow', 'NewItems', rowIndex )
          if nargin ~= 4 
            error ( 'uimulticollist:changeRow', 'Change row requires 4 input arguments' );
          end
          rowIndex = varargin{4};
          if rowIndex < 1 || rowIndex > nRows
            error ( 'uimulticollist:changeRow', 'User requested to change a row < or > data in list' );
          end
          if length(varargin{3}) == nCols
            for ii=1:nCols
              s.string{rowIndex,ii} = varargin{3}{ii};
            end;
          else
            error ( 'uimulticollist:changeRow', 'New row does not contain the correct number of items, input %i, required %i', length(varargin{2}), nCols );
          end
          s.string = EnsureAllStringInputCorrect ( s.string );
        case 'columnColour'       % ( h, 'columnColor', colIndex, 'Colour' )
          if nargin ~= 4 
            error ( 'uimulticollist:columnColour', 'Change colour requires 4 input arguments' );
          end
          s = CheckColourVariable ( s, nCols );
          s.colours{varargin{3}} = varargin{4};
        case 'delCol'             % ( h, 'delCol', colIndex )
          if nargin ~= 3 
            error ( 'uimulticollist:delCol', 'Delete column requires 3 input arguments' );
          end
          if varargin{3} <= 0
            s.string(:,end+varargin{3}) = [];
            s.colours(end+varargin{3}) = [];
            if isfield ( s, 'filters' );
            end
          else
            s.string(:,varargin{3}) = [];
            s.colours(varargin{3}) = [];
          end
        case 'delRow'             % ( h, 'delRow', rowIndex )
          if nargin ~= 3 
            error ( 'uimulticollist:delRow', 'Delete row requires 3 input arguments' );
          end
          if varargin{3} <= 0
            s.string(end+varargin{3},:) = [];
          else
            s.string(varargin{3},:) = [];
          end
        case 'filtersApplied'     % ( h, 'filtersApplied )
          varargout{1} = s.filterDetails;
          return
        case 'nCols'              % ( h, 'nCols' )
          varargout{1} = nCols;
        case 'nRows'              % ( h, 'nRows' )
          varargout{1} = nRows;
        case 'removeUIFilter'     % ( h, 'removeUIFilter' )
          fMenu = get ( s.uic, 'children' );
          fMenuIndex = strcmp ( get ( fMenu, 'tag' ), s.filterMenuName );
          delete ( get ( fMenu(fMenuIndex), 'children' ));
          delete ( fMenu(fMenuIndex) );
          if length(fMenu) == 1 
            delete ( s.uic );
          end
          s = rmfield ( s, { 'uic' 'filterMenuName' 'uicRow1Labels'} );
        case 'resetFilter'        % ( h, 'resetFilter' )
          if isfield ( s, 'filters' );        s = rmfield ( s, 'filters' ); end
          if isfield ( s, 'filterDetails' ) ; s = rmfield ( s, 'filterDetails' ); end
        case 'selectionInNumbers' % ( h, 'seletionInNumbers', colIndex )
          val = get ( h, 'value' );
          if nargin==2
            varargout{1} = cellfun ( @str2double, s.string(val,:));
          else
            varargout{1} = cellfun ( @str2double, s.string(val,varargin{3}));
          end 
        case 'selectedStrCol'     % ( h, 'selectedStrCol', colIndex )
          if nargin ~= 3 
            error ( 'uimulticollist:selectedStrCol', 'selectedStrCol requires 3 input arguments' );
          end
          val = get ( h, 'value' );
          colIndex = varargin{3};
          if length(val) > 1
            for jj=1:length(val)
              for kk=1:length(colIndex)
                varargout{1}{jj,kk} = s.string{val(jj),colIndex(kk)};
              end
            end
          else
%             varargout{1} = {s.string{val,varargin{3}}};
            varargout{1} = s.string(val,varargin{3});
          end
        case 'selectedString'     % ( h, 'selectedString' )
          val = get ( h, 'value' );
          if length(val) > 1
            varargout{1} = s.string(val,:);
          else
%             varargout{1} = {s.string{val,:}};
            varargout{1} = s.string(val,:);
          end
        case 'separator'          % ( h, 'separator', 'newColSeparator' )
          if nargin ~= 3 
            error ( 'uimulticollist:separator', 'separator requires 3 input arguments' );
          end
          s.separator = varargin{3};
        case 'setRow1Header'      % ( h, 'setRow1Header', OnOffFlag )
          if nargin ~= 3 
            error ( 'uimulticollist:setRow1Header', 'setRow1Header requires 3 input arguments' );
          end
          s.uicRow1Labels = varargin{3};
        case 'string'             % ( h, 'string', cellArrayStrings )
          if iseven
            varargout{1} = s.string;
          else
            s.string = varargin{3};
          end
        case 'value'              % ( h, 'value', selectionValue )
          if iseven
            varargout{1} = get ( h, 'value' );
          else
            set ( h, 'value', varargin{3} );
          end
        otherwise
          text = ParseFileForCase ( which ( mfilename ) );       
          error ( 'uimulticollist:Unsupported', 'Unsupported request "%s"\n Suitable methods are:\n%s', varargin{2}, text )
      end
      set ( h, 'ApplicationData', s );
    else
      error ( 'uimulticollist:InputError', 'User must provide extra inputs to interact with a live listbox' );
    end
  end
  % has filters been applied?
  if isfield ( s, 'filters' );
    s = get ( h, 'ApplicationData' );
    s.filters = ones(nRows,nCols);
    for i=1:length(s.filterDetails)
      for j=1:length(s.filterDetails{i})
        s.filters(~strcmp ( s.string(:,i), s.filterDetails{i}{j} ),i) = 0; 
      end
    end
    set ( h, 'ApplicationData', s );
  end
  % Update the html string in the listbox
  UpdateStringHTML(h);
  % if a uicontextmenu -> update it.
  if isfield ( s, 'uic' )
    UpdateFilterMenu ( h );
  end
end
%% Update the filter menu
function UpdateFilterMenu ( h )
  s = get ( h, 'ApplicationData' );
  [nRows,nCols]= size(s.string);
  if ~isfield ( s, 'uic' ); return; end
  fMenu = get ( s.uic, 'children' );
  fMenuIndex = strcmp ( get ( fMenu, 'tag' ), s.filterMenuName );
  delete ( get ( fMenu(fMenuIndex), 'children' ));
  for i=1:nCols
    if s.uicRow1Labels
      label = s.string{1,i};
    else
      label = sprintf ( 'Col %i', i );
    end
    uimenu ( fMenu(fMenuIndex), 'label', label, 'Callback', {@FilterTableData h i} );
  end
  if isfield ( s, 'filters' )
    uimenu ( fMenu(fMenuIndex), 'label', 'Reset Filters', 'Callback', {@FilterTableData h 0}, 'separator', 'on' );
  end    
end
%% Filter the string data
function FilterTableData ( obj, eventdata, h, colIndex )
  s = get ( h, 'ApplicationData' );
  [nRows,nCols]= size(s.string);
  if colIndex == 0 % removing the filters
    s = rmfield ( s, { 'filters' 'filterDetails' } );
  else
    if ~isfield ( s, 'filters' )
      s.filters = ones(nRows,nCols);
    end
    if ~isfield ( s, 'uicRow1Labels' ); s.uicRow1Labels = 0; end
    if s.uicRow1Labels
      items = s.string(2:end,colIndex);
    else
      items = s.string(:,colIndex);
    end
    items = ['None'; unique(items)];
    try
      [dummy,initValue]=intersect ( items, s.filterDetails{colIndex} ); %#ok<*ASGLU>
    catch %#ok<CTCH>  Ok to crash - if no matches found or no filter applied.
      initValue = 1;
    end
    if ~isnan ( eventdata )
      [a,b] = listdlg('PromptString','Filter Items','SelectionMode','multi','ListString',items, 'initialValue', initValue);
      if b == 0; return; end
    else
      obj = CheckInputIsCell ( obj );
      filterValue= EnsureAllStringInputCorrect (obj);
      
      [dummy,a] = intersect ( items, filterValue );
    end
    if length(a)==1 && a == 1 % resetting the filter.
      s.filterDetails{colIndex} = [];
      s.filters(:,colIndex) = 1;
    else
      s.filterDetails{colIndex} = items(a);
      s.filters(:,colIndex) = 0;
      for ii=1:length(a)
        s.filters(strcmp ( s.string(:,colIndex), items(a(ii)) ),colIndex) = 1; 
      end
%       s.filters(:,colIndex)=s.filters(:,colIndex)==2;
    end
    if s.uicRow1Labels
      s.filters(1,colIndex) = 1;
    end
    % check the filtering hasn't make the listbox disappear
    val = get ( h, 'value' );
    val(val > length(find(sum(s.filters,2)==6))) = [];
    if isempty ( val );
      val = 1;
    end
     set ( h, 'value', val );
  end
  
  set ( h, 'ApplicationData', s );
  UpdateStringHTML( h );
  UpdateFilterMenu ( h );
end
%% Update the text to HTML format.
function UpdateStringHTML (h)
  s = get ( h, 'ApplicationData' );
  str = s.string;
  sz=size(str);
  nRow = sz(1);
  % work out the width of each column
  nCol = sz(2);
  colWidths = zeros(nCol,1);
  colDefined = iscell( s.colours );
  if colDefined
    if length(s.colours) ~= nCol
      error ( 'uimulticollist:IncorrectInput', 'length of column colours (%i) does not match the number of columns (%i)', length(s.colours), nCol )
    end
  end
  for ii=1:nCol
    for jj=1:nRow
      colWidths(ii) = max ( colWidths(ii), length(str{jj,ii}) );
    end
  end
  %
%   if s.headedList
%     s.header = ConvertStringToHTML ( s.header, ' ', colWidths );
%   end
  details = cell(nRow,1);
  for i=1:nRow
    details{i} = '<HTML>';
    if colDefined
      details{i} = ConvertStringToHTML ( str(i,:), s.separator, colWidths, s.colours );
    else
      details{i} = ConvertStringToHTML ( str(i,:), s.separator, colWidths );
    end
  end
  if isfield ( s, 'filters' )
    details(~logical(sum(s.filters,2)==nCol)) = [];
    % TODO: If filter wipes out seleted values?
  end
  s.listboxString = details;
  set ( h, 'string', details, 'fontname', 'monospaced', 'ApplicationData', s );
end
%% HTML string conversion
function output = ConvertStringToHTML ( str, sep, colWidths, colours )
  nCol = length(str);
  colDefined = nargin == 4;
  output = '<HTML>';
  for j=1:nCol
    text = str{j};
    text = LeftFormatString(text,colWidths(j));
    if colDefined
      text = sprintf ( '<FONT color="%s">%s</FONT>', colours{j}, text );
    end
    if j==1
      output = sprintf ( '%s%s', output, text );
    else
      output = sprintf ( '%s%s%s', output, sep, text );
    end
    
  end
  output = sprintf ( '%s</HTML>', output );
end
%% format the string so its left formatted
function output = LeftFormatString(input, nLen)
  formatString = sprintf ( '%%%is', nLen );
  output = fliplr(sprintf ( formatString, fliplr(input) ));
  spaces = length(input)+1;
  output = sprintf ( '%s%s', output(1:spaces-1), regexprep ( output(spaces:end), ' ', '&nbsp' ) );
end
%% convert the string to a double.
function str = double2CompactStr ( inputStr )
  if iscell ( inputStr )
    inputStr = cell2mat( inputStr );
  end
  str = '';
  for i=1:length(inputStr)
    t_str = sprintf ( '%f', inputStr(i) );
    n = regexp(t_str,'\.0*$');
    if ~isempty(n)
      t_str(n:end) = [];
    else
      m = regexp(t_str,'0*$');
      if ~isempty(m)
        t_str(m:end) = [];
      end
    end
    str = sprintf ( '%s %s', str, t_str );
  end
  str = strtrim ( str );
end
%% check the input strings.
function cString = EnsureAllStringInputCorrect ( cString )
%   nCols = length(cString{1});
  % check all sub items are cell arrays
%   for ii=1:length(cString)
%     if ~iscell ( cString{ii} )
%       error ( 'uimulticollist:IncorrectInput', 'Each row needs to be a cell array' )
%     end
%   end
  % check the input data is correct
%   for ii=1:length(cString)
%     if length(cString{ii})~=nCols
%       error ( 'uimulticollist:IncorrectInput', 'All Input data does not contain the same number of rows' )
%     end
%   end
  % check all items in the cells are strings
  s = size(cString);
  for ii=1:s(1)%length(cString)
    for jj=1:s(2)%length(cString{ii})
      if isnumeric ( cString{ii,jj} )
        cString{ii,jj} = double2CompactStr ( cString{ii,jj} );
      elseif ischar ( cString{ii,jj} )
      elseif islogical ( cString{ii,jj} )
        cString{ii,jj} = double2CompactStr ( cString{ii,jj} );
      else
        error ( 'uimulticollist:IncorrectInput', 'Each item in a row can be numeric or char only' )
      end
    end
  end  
end
%% check the input for any special pairs
function [output data] = CheckForSpecialInputs ( name, output, data )
  colIndex = strcmp ( output(1:2:end), name );
  if max ( colIndex ) == 1
    vIndex = 2*find(colIndex==1);
    data = output{(vIndex)};
    output(vIndex-1) = [];
    output(vIndex-1) = []; % remove it from the varargin passed to the uicontrol
  end
end
%% parse file for help
function text = ParseFileForCase ( filename )
  fid = fopen ( filename, 'r' );
  text = '';
  if fid ~= -1
    while ( true )
      line = fgetl ( fid );
      if line == -1; break; end
      if length(line)>4
        if ~isempty ( strfind ( line, sprintf ( 'ca%s', 'se' ) ) )
          trimLine = strtrim ( line );
          if ~strcmp ( trimLine(1), '%' )
            if nargout == 0
              disp ( line );
            else
              text = sprintf ( '%s\n%s', text, line );
            end
          end
        end        
      end
    end
  end
end
%% check the input is a cell
function newData = CheckInputIsCell ( newData )
  if ~iscell(newData)
    for jj=1:length(newData)
      temp{jj,1} = newData(jj); %#ok<AGROW>
    end
    newData = temp;
  end
end
function s = CheckColourVariable ( s, nCols )
  if ~iscell ( s.colours ) % if it isn't set - make them all black;
    s = rmfield ( s, 'colours' );
    for jj=1:nCols;
      s.colours{jj} = 'BLACK';
    end
  end
end
%%  %% some initial coding for the headed listbox....
%     s.headedList = headedList;
%     s.lastIndex = 1;
%     if ~iscell ( header )
%       s.header = {header};
%     else
%       s.header = header;
%     end
%     if s.headedList
%       
%       if exist ( 'findjobj', 'file' ) ~= 2
%         error ( 'uimulticollist:headedList', 'Required function findjobj not found - can be downloaded from matlab FEX' );
%       end
%       jScrollPane = findjobj(h,'propert',{'tag',tag});
%       if length ( jScrollPane ) ~= 1
%         error ( 'uimulticollist:findjobj', 'findjobj found to many objects to work on' );
%       end
%       if isempty ( jScrollPane ) % this happens if the figure is invisible.
%         figHandle = ancestor ( parent, 'figure' );
%         s = get ( 0, 'screensize' );
%         p = get ( figHandle, 'position' );
%         set ( figHandle, 'position', [0-p(3) 0-s(4), 0.5*p(3), 0.5*p(4)], 'visible', 'on' );
%         jScrollPane = findjobj(tmp);
%         set ( figHandle, 'visible', 'off', 'position', p ); % put it back and turn it off again.
%       end
%       jScroll=handle ( jScrollPane, 'CallbackProperties' );
%       set ( jScroll, 'AdjustmentValueChangedCallback', {@UpdateHeadedListCallback h} );
%     end      
