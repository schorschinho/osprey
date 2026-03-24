%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright:
% Jun Tan
% University of Texas Southwestern Medical Center
% Department of Radiation Oncology
% Last edited: 08/19/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = line_measurement(varargin)
% LINE_MEASUREMENT MATLAB code for line_measurement.fig
%      LINE_MEASUREMENT, by itself, creates a new LINE_MEASUREMENT or raises the existing
%      singleton*.
%
%      H = LINE_MEASUREMENT returns the handle to a new LINE_MEASUREMENT or the handle to
%      the existing singleton*.
%
%      LINE_MEASUREMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LINE_MEASUREMENT.M with the given input arguments.
%
%      LINE_MEASUREMENT('Property','Value',...) creates a new LINE_MEASUREMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before line_measurement_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to line_measurement_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help line_measurement

% Last Modified by GUIDE v2.5 16-Aug-2014 00:56:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @line_measurement_OpeningFcn, ...
    'gui_OutputFcn',  @line_measurement_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before line_measurement is made visible.
function line_measurement_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to line_measurement (see VARARGIN)

% Choose default command line output for line_measurement
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

parse_args(hObject, varargin);

setappdata(hObject, 'entry_update_data', @parse_args);

set(hObject, 'Visible', getappdata(hObject, 'initialVisible'));


% -------------------------------------------------------------------
function parse_args(hFig, args)

handles = guidata(hFig);

numArgs = length(args);

assert(1 == numArgs || 2 == numArgs || 3 == numArgs, 'Must have 1, 2, or 3 arguments');

setappdata(handles.figure_lp, 'xyLims', []);
setappdata(handles.figure_lp, 'mainFigHandle', []);
setappdata(handles.figure_lp, 'initialVisible', 'on');

if 1 == numArgs
    
    if is_handle(args{1})
        setappdata(handles.figure_lp, 'mainFigHandle', args{1});
        setappdata(handles.figure_lp, 'initialVisible', 'off');
    else
        parse_line_data(handles, args{1});
    end
    
elseif 2 == numArgs
    
    parse_line_data(handles, args{1})
    
    if is_handle(args{2})
        setappdata(handles.figure_lp, 'mainFigHandle', args{2});
    else
        parse_lim_data(handles, args{2});
    end
    
else % if 3 == numArgs
    
    parse_line_data(handles, args{1})
    
    parse_lim_data(handles, args{2});
    
    assert(is_handle(args{3}), 'The third argument must be a valid GUI handle.');
    setappdata(handles.figure_lp, 'mainFigHandle', args{3});
    
end

hSel = getappdata(handles.figure_lp, 'hSelectedPoint');
if ~isempty(hSel)
    delete(hSel);
    setappdata(handles.figure_lp, 'hSelectedPoint', []);
end

update_measurement(handles);


% -------------------------------------------------------------------
function parse_line_data(handles, lineData)

% lineData: Each row is 1 point. First number is value.
% lineType:
% 1: coordinate is row number.
% 2: coordinate is 1D (x).
% 3: coordinate is 2D (x, y).

assert(isnumeric(lineData) || islogical(lineData), 'Data type must be numeric or logical.');

pointValues = lineData(:, 1);
pointPos = lineData(:, 2:end);

numPoints = numel(pointValues);
assert(numPoints > 1, 'Must have at least 2 points.');

numPosDims = size(pointPos, 2);
assert(numPosDims <= 2, 'Dimenson ust be 1D or 2D.');

if 0 == numPosDims
    pointPos = (1 : numPoints)';
    lineType = 1;
elseif 1 == numPosDims
    lineType = 1;
else
    lineType = 2;
end

assert(2 == lineType || 1 == lineType && all(diff(pointPos) > 0), ...
    'Coordinates must be monotonically increasing.');

setappdata(handles.figure_lp, 'lineType', lineType);
setappdata(handles.figure_lp, 'pointValues', pointValues);
setappdata(handles.figure_lp, 'pointPos', pointPos);

get_data_format(handles, pointValues, pointPos);


% -------------------------------------------------------------------
function parse_lim_data(handles, xyLims)

assert(isnumeric(xyLims), 'Coordinate limits must be numeric.');
assert(4 == numel(xyLims) && numel(xyLims) == length(xyLims) ...
    && xyLims(2) > xyLims(1) && xyLims(4) > xyLims(3), ...
    'X and Y lims must be [xLower xUpper yLower yUpper] for .');

setappdata(handles.figure_lp, 'xyLims', xyLims);


% -------------------------------------------------------------------
function get_data_format(handles, pointValues, pointPos)

if all(0 == rem(pointValues, 1))
    valueFormat = '%.0f';
    valueMaxExp = 0;
else
    [valueFormat, valueMaxExp] = find_float_text_format(max(abs(pointValues(:))));
end

if all(0 == rem(pointPos(:), 1))
    posFormat = '%.0f';
    posMaxExp = 0;
else
    [posFormat, posMaxExp] = find_float_text_format(max(abs(pointPos(:))));
end

setappdata(handles.figure_lp, 'valueFormat', valueFormat);
setappdata(handles.figure_lp, 'valueMaxExp', valueMaxExp);
setappdata(handles.figure_lp, 'posFormat', posFormat);
setappdata(handles.figure_lp, 'posMaxExp', posMaxExp);


% -------------------------------------------------------------------
function [fmt, maxExp] = find_float_text_format(v)

v = abs(v); % Force v to be >= 0, though not always necessary.

if  v >= 1e3 || v <= 0
    fmt = '%.3e';
    maxExp = floor(log10(v));
elseif v >= 1e2
    fmt = '%.1f';
    maxExp = 0;
elseif v >= 1e1
    fmt = '%.2f';
    maxExp = 0;
else
    fmt = '%.3f';
    maxExp = 0;
end


% -------------------------------------------------------------------
function s = format_float(v, fmt, maxExp)

if fmt(end) == 'e'
    v = v / (10 ^ maxExp);
    s = [sprintf('%.3f', v) 'e+03'];
else
    s = sprintf(fmt, v);
end

if '-' ~= s(1)
    s = [' ' s];
end


% -------------------------------------------------------------------
function update_measurement(handles)

stats = cell(6, 2);
stats{1, 1} = ' #Points';
stats{2, 1} = ' Length';
stats{3, 1} = ' Max';
stats{4, 1} = ' Min';
stats{5, 1} = ' Mean';
stats{6, 1} = ' SD';

pointValues = getappdata(handles.figure_lp, 'pointValues');

stats{1, 2} = length(pointValues);

if stats{1, 2} > 0
    stats{3, 2} = max(pointValues);
    stats{4, 2} = min(pointValues);
    if ~isempty(pointValues)
        stats{5, 2} = mean(pointValues);
    end
    
    if ~isempty(pointValues)
        stats{6, 2} = std(pointValues);
    end
    
    pointPos = getappdata(handles.figure_lp, 'pointPos');
    
    if 1 == getappdata(handles.figure_lp, 'lineType')
        set(handles.checkbox_stretch, 'Visible', 'off');
        stats{2, 2} = range(pointPos);
        update_2d_profile(handles);
    else % if 2 == lineType
        set(handles.checkbox_stretch, 'Visible', 'on');
        numSegments = length(pointPos) - 1;
        segmentLength = zeros(1, numSegments);
        for i = 1 : numSegments
            segmentLength(i) = pdist2(pointPos(i, :), pointPos(i+1, :));
        end
        stats{2, 2} = sum(segmentLength);
        setappdata(handles.figure_lp, 'segmentLength', segmentLength);
        update_3d_profile(handles);
    end
    
    valueFormat = getappdata(handles.figure_lp, 'valueFormat');
    valueMaxExp = getappdata(handles.figure_lp, 'valueMaxExp');
    valueText = arrayfun(@(x)format_float(x, valueFormat, valueMaxExp), ...
        pointValues, 'UniformOutput', false);
    posFormat = getappdata(handles.figure_lp, 'posFormat');
    posMaxExp = getappdata(handles.figure_lp, 'posMaxExp');
    posText = arrayfun(@(x)format_float(x, posFormat, posMaxExp), ...
        pointPos, 'UniformOutput', false);
    set(handles.uitable_points, 'Data', [valueText posText]);
end

set(handles.uitable_stats, 'Data', stats);

check_grid(handles);
check_marker(handles);


% -------------------------------------------------------------------
function update_2d_profile(handles)

pointValues = getappdata(handles.figure_lp, 'pointValues');
pointPos = getappdata(handles.figure_lp, 'pointPos');

h = plot(handles.axes_lineProfile, pointPos, pointValues, 'linewidth', 2);
setappdata(handles.figure_lp, 'hProfile', h);
rotate3d(handles.axes_lineProfile, 'off');
box(handles.axes_lineProfile, 'off');
xlabel(handles.axes_lineProfile, 'x');
ylabel(handles.axes_lineProfile, 'Value');
setappdata(handles.figure_lp, 'segmentLength', []);

set(handles.uitable_points, 'ColumnName', {'Value', 'x'});
set(handles.uitable_points, 'ColumnWidth', {65, 40});


% -------------------------------------------------------------------
function update_3d_profile(handles)

pointValues = getappdata(handles.figure_lp, 'pointValues');
pointPos = getappdata(handles.figure_lp, 'pointPos');

if isempty(pointValues) || isempty(pointPos)
    return;
end

if 1 == get(handles.checkbox_stretch, 'Value')
    segmentLength = getappdata(handles.figure_lp, 'segmentLength');
    h = plot(handles.axes_lineProfile, [0 cumsum(segmentLength)], pointValues, 'linewidth', 1);
    setappdata(handles.figure_lp, 'hProfile', h);
    rotate3d(handles.axes_lineProfile, 'off');
    box(handles.axes_lineProfile, 'off');
    xlabel(handles.axes_lineProfile, 'Cumulative distance from first point');
    ylabel(handles.axes_lineProfile, 'Value');
else
    h = plot3(handles.axes_lineProfile, pointPos(:, 1), pointPos(:, 2), pointValues, 'linewidth', 1);
    setappdata(handles.figure_lp, 'hProfile', h);
    rotate3d(handles.axes_lineProfile, 'on');
    xyLims = getappdata(handles.figure_lp, 'xyLims');
    if isempty(xyLims)
        xlim(handles.axes_lineProfile, 'auto');
        ylim(handles.axes_lineProfile, 'auto');
    else
        xlim(handles.axes_lineProfile, xyLims(1:2));
        ylim(handles.axes_lineProfile, xyLims(3:4));
    end
    box(handles.axes_lineProfile, 'off');
    xlabel(handles.axes_lineProfile, 'x');
    ylabel(handles.axes_lineProfile, 'y');
    zlabel(handles.axes_lineProfile, 'Value');
end

set(handles.uitable_points, 'ColumnName', {'Value', 'x', 'y'});
set(handles.uitable_points, 'ColumnWidth', {50, 40, 40});


% -------------------------------------------------------------------
function check_grid(handles)

if 1 == get(handles.checkbox_grid, 'Value')
    grid(handles.axes_lineProfile, 'on');
else
    grid(handles.axes_lineProfile, 'off');
end


% -------------------------------------------------------------------
function check_marker(handles)

h = getappdata(handles.figure_lp, 'hProfile');

if 1 == get(handles.checkbox_marker, 'Value')
    set(h, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
else
    set(h, 'Marker', 'none');
end


% --- Outputs from this function are returned to the command line.
function varargout = line_measurement_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_stretch.
function checkbox_stretch_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_stretch

hSel = getappdata(handles.figure_lp, 'hSelectedPoint');
if is_handle(hSel)
    delete(hSel);
end
setappdata(handles.figure_lp, 'hSelectedPoint', []);

update_3d_profile(handles);
check_grid(handles);
check_marker(handles);


% --- Executes on button press in checkbox_grid.
function checkbox_grid_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_grid

check_grid(handles);


% --- Executes on button press in checkbox_marker.
function checkbox_marker_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_marker

check_marker(handles);


% --- Executes when selected cell(s) is changed in uitable_points.
function uitable_points_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable_points (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata.Indices)
    return;
end

r = eventdata.Indices(1);
hProfile = getappdata(handles.figure_lp, 'hProfile');

if isempty(hProfile)
    return;
end

xd = get(hProfile, 'XData');
yd = get(hProfile, 'YData');
zd = get(hProfile, 'ZData');

lineType = getappdata(handles.figure_lp, 'lineType');

hSel = getappdata(handles.figure_lp, 'hSelectedPoint');
if isempty(hSel)
    hold(handles.axes_lineProfile, 'on');
    if 1 == lineType || 1 == get(handles.checkbox_stretch, 'Value')
        hSel = plot(handles.axes_lineProfile, xd(r), yd(r));
    else
        hSel = plot3(handles.axes_lineProfile, xd(r), yd(r), zd(r));
    end
    set(hSel, 'Marker', 'o', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'g');
    hold(handles.axes_lineProfile, 'off');
    setappdata(handles.figure_lp, 'hSelectedPoint', hSel);
else
    set(hSel, 'XData', xd(r), 'YData', yd(r));
    if 2 == lineType && 0 == get(handles.checkbox_stretch, 'Value')
        set(hSel, 'ZData', zd(r));
    end
end


% --------------------------------------------------------------------
function uitable_points_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitable_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hSel = getappdata(handles.figure_lp, 'hSelectedPoint');

if ~isempty(hSel)
    delete(hSel);
    setappdata(handles.figure_lp, 'hSelectedPoint', []);
end


% --- Executes when user attempts to close figure_lp.
function figure_lp_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_lp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainFigHandle = getappdata(handles.figure_lp, 'mainFigHandle');
if is_handle(mainFigHandle)
    hFun = getappdata(mainFigHandle, 'hFunCallbackLineMeasurementClosed');
    if isa(hFun, 'function_handle')
        feval(hFun, mainFigHandle);
    end
else
    delete(hObject);
end

