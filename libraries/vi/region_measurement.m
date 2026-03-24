%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright:
% Jun Tan
% University of Texas Southwestern Medical Center
% Department of Radiation Oncology
% Last edited: 08/19/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = region_measurement(varargin)
% REGION_MEASUREMENT MATLAB code for region_measurement.fig
%      REGION_MEASUREMENT, by itself, creates a new REGION_MEASUREMENT or raises the existing
%      singleton*.
%
%      H = REGION_MEASUREMENT returns the handle to a new REGION_MEASUREMENT or the handle to
%      the existing singleton*.
%
%      REGION_MEASUREMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGION_MEASUREMENT.M with the given input arguments.
%
%      REGION_MEASUREMENT('Property','Value',...) creates a new REGION_MEASUREMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before region_measurement_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to region_measurement_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help region_measurement

% Last Modified by GUIDE v2.5 19-Aug-2014 22:40:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @region_measurement_OpeningFcn, ...
    'gui_OutputFcn',  @region_measurement_OutputFcn, ...
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

% --- Executes just before region_measurement is made visible.
function region_measurement_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to region_measurement (see VARARGIN)

% Choose default command line output for region_measurement
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

setappdata(handles.figure_rs, 'regionType', 1); % 1: rectangle, 2: disc.
setappdata(handles.figure_rs, 'mainFigHandle', []);
setappdata(handles.figure_rs, 'initialVisible', 'on');

if 1 == numArgs
    
    if is_handle(args{1})
        setappdata(handles.figure_rs, 'mainFigHandle', args{1});
        setappdata(handles.figure_rs, 'initialVisible', 'off');
    else
        parse_region_data(handles, args{1});
    end
    
elseif 2 == numArgs
    
    if is_handle(args{2})
        setappdata(handles.figure_rs, 'mainFigHandle', args{2});
    else
        assert(strcmpi('Rectangle', args{2}) || strcmpi('Disc', args{2}) || strcmpi('Drawn', args{2}), ...
            'Region shape must be disc or rectangle or user-drawn.');
        if strcmpi('Disc', args{2})
            setappdata(handles.figure_rs, 'regionType', 2);
        elseif strcmpi('Drawn', args{2})
            setappdata(handles.figure_rs, 'regionType', 3);
        end
    end
    
    parse_region_data(handles, args{1});

else % if 3 == numArgs
    assert(strcmpi('Disc', args{2}) || strcmpi('Rectangle', args{2}) || strcmpi('Drawn', args{2}), ...
        'Region shape must be disc or rectangle or user-drawn.');
    if strcmpi('Disc', args{2})
        setappdata(handles.figure_rs, 'regionType', 2);
    elseif strcmpi('Drawn', args{2})
        setappdata(handles.figure_rs, 'regionType', 3);
    end
    
    assert(is_handle(args{3}), 'The third argument must be a valid GUI handle.');
    setappdata(handles.figure_rs, 'mainFigHandle', args{3});
    
    parse_region_data(handles, args{1});
    
end

update_measurement(handles);


% -------------------------------------------------------------------
function parse_region_data(handles, regionData)

regionType = getappdata(handles.figure_rs, 'regionType');

if 1 == regionType || 2 == regionType
    assert((isnumeric(regionData) || islogical(regionData)) && ismatrix(regionData), ...
        'Data must be gray or binary 2D image.');
    setappdata(handles.figure_rs, 'regionData', double(regionData));
else %3 == regionType
    assert(isa(regionData, 'cell') && isa(regionData{1}, 'double') && isa(regionData{2}, 'logical') ....
        && isequal(size(regionData{1}), size(regionData{2})));
    setappdata(handles.figure_rs, 'regionData', regionData{1});
    setappdata(handles.figure_rs, 'regionMask', regionData{2});
end


% -------------------------------------------------------------------
function update_measurement(handles)

regionData = getappdata(handles.figure_rs, 'regionData');
regionType = getappdata(handles.figure_rs, 'regionType');
pixels = [];

if 1 == regionType || 2 == regionType
    
    data = cell(6, 2);
    data{1, 1} = ' Shape';
    data{2, 1} = ' #Points';
    data{3, 1} = ' Area';
    data{4, 1} = ' Width';
    data{5, 1} = ' Height';
    data{6, 1} = ' Max';
    data{7, 1} = ' Min';
    data{8, 1} = ' Mean';
    data{9, 1} = ' SD';
    
    if 1 == regionType
        data{1, 2} = ' Rectangle';
        pixels = regionData(:);
    elseif 2 == regionType
        data{1, 2} = ' Disc';
        h = size(regionData, 1);
        w = size(regionData, 2);
        cy = h / 2.0 + 0.5;
        cx = w / 2.0 + 0.5;
        [x, y] = meshgrid(1:w, 1:h);
        dx = abs(x - cx);
        dy = abs(y - cy);
        c = ((dx .^ 2) / (cx .^ 2) + (dy .^ 2) ./ (cy .^ 2)) <= 1;
        pixels = regionData(c);
    end
    
    data{2, 2} = numel(pixels);
    
    if data{2, 2} > 0
        data{3, 2} = numel(pixels);
        data{4, 2} = size(regionData, 2);
        data{5, 2} = size(regionData, 1);
        data{6, 2} = max(pixels);
        data{7, 2} = min(pixels);
        data{8, 2} = mean(pixels);
        data{9, 2} = std(pixels);
    end
    
else %if 3 == regionType
    
    data = cell(16, 2);
    data{1, 1} = ' Shape';
    data{2, 1} = ' Area';
    data{3, 1} = ' Centroid.X';
    data{4, 1} = ' Centroid.Y';
    data{5, 1} = ' W.Cent.X';
    data{6, 1} = ' W.Cent.Y';
    data{7, 1} = ' Orientation';
    data{8, 1} = ' Major Axis';
    data{9, 1} = ' Minor Axis';
    data{10, 1} = ' Equiv Diam';
    data{11, 1} = ' Perimeter';
    data{12, 1} = ' Solidity';
    data{13, 1} = ' Max';
    data{14, 1} = ' Min';
    data{15, 1} = ' Mean';
    data{16, 1} = ' SD';

    img = getappdata(handles.figure_rs, 'regionData');
    mask = getappdata(handles.figure_rs, 'regionMask');
    stats = regionprops(mask, img, ...
        {'Area', 'Centroid', 'Orientation', 'MajorAxisLength', 'MinorAxisLength', ...
        'EquivDiameter', 'Perimeter', 'Solidity', 'MaxIntensity', 'MinIntensity', ...
        'MeanIntensity', 'WeightedCentroid'});

    if length(stats) >= 1
        stats = stats(1);
        pixels = img(mask);
        
        data{1, 2} = 'User-Drawn';
        data{2, 2} = stats.Area;
        data{3, 2} = stats.Centroid(1);
        data{4, 2} = stats.Centroid(2);
        data{5, 2} = stats.WeightedCentroid(1);
        data{6, 2} = stats.WeightedCentroid(2);
        data{7, 2} = stats.Orientation;
        data{8, 2} = stats.MajorAxisLength;
        data{9, 2} = stats.MinorAxisLength;
        data{10, 2} = stats.EquivDiameter;
        data{11, 2} = stats.Perimeter;
        data{12, 2} = stats.Solidity;
        data{13, 2} = stats.MaxIntensity;
        data{14, 2} = stats.MinIntensity;
        data{15, 2} = stats.MeanIntensity;
        data{16, 2} = std(pixels);
    end
        
end

set(handles.uitable_stats, 'Data', data);

if ~isempty(pixels)
    hist(handles.axes_hist, pixels, 100);
    box(handles.axes_hist, 'off');
    xlabel(handles.axes_hist, 'Pixel value');
    ylabel(handles.axes_hist, 'Counts');
end

check_grid(handles);


% -------------------------------------------------------------------
function check_grid(handles)

if 1 == get(handles.checkbox_grid, 'Value')
    grid(handles.axes_hist, 'on');
else
    grid(handles.axes_hist, 'off');
end


% --- Outputs from this function are returned to the command line.
function varargout = region_measurement_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_grid.
function checkbox_grid_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_grid

check_grid(handles);


% --- Executes when user attempts to close figure_rs.
function figure_rs_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_rs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainFigHandle = getappdata(handles.figure_rs, 'mainFigHandle');
if is_handle(mainFigHandle)
    hFun = getappdata(mainFigHandle, 'hFunCallbackRegionMeasurementClosed');
    if isa(hFun, 'function_handle')
        feval(hFun, mainFigHandle);
    end
else
    delete(hObject);
end
