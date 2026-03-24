%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright:
% Jun Tan
% University of Texas Southwestern Medical Center
% Department of Radiation Oncology
% Last edited: 08/19/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = image_stats(varargin)
% IMAGE_STATS MATLAB code for image_stats.fig
%      IMAGE_STATS, by itself, creates a new IMAGE_STATS or raises the existing
%      singleton*.
%
%      H = IMAGE_STATS returns the handle to a new IMAGE_STATS or the handle to
%      the existing singleton*.
%
%      IMAGE_STATS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE_STATS.M with the given input arguments.
%
%      IMAGE_STATS('Property','Value',...) creates a new IMAGE_STATS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before image_stats_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to image_stats_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help image_stats

% Last Modified by GUIDE v2.5 11-Aug-2014 20:13:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @image_stats_OpeningFcn, ...
    'gui_OutputFcn',  @image_stats_OutputFcn, ...
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

% --- Executes just before image_stats is made visible.
function image_stats_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to image_stats (see VARARGIN)

% Choose default command line output for image_stats
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

assert(1 == numArgs || 2 == numArgs, 'Must have 1 or 2 arguments');

setappdata(handles.figure_is, 'mainFigHandle', []);
setappdata(handles.figure_is, 'initialVisible', 'on');

parse_image_data(handles, args{1});

if 2 == numArgs
    assert(is_handle(args{2}));
    setappdata(handles.figure_is, 'initialVisible', 'off');
    setappdata(handles.figure_is, 'mainFigHandle', args{2});
end

update_stats(handles);


% -------------------------------------------------------------------
function parse_image_data(handles, imageData)

assert((isnumeric(imageData) || islogical(imageData)) && ismatrix(imageData), ...
    'Data must be gray or binary 2D image.');
setappdata(handles.figure_is, 'imageData', double(imageData));


% -------------------------------------------------------------------
function update_stats(handles)

imageData = getappdata(handles.figure_is, 'imageData');

data = imageData(:);
    
hist(handles.axes_hist, data, 100);
box(handles.axes_hist, 'off');
xlabel(handles.axes_hist, 'Pixel value');
ylabel(handles.axes_hist, 'Counts');

stats = cell(6, 2);

stats{1, 1} = ' Area';
stats{1, 2} = numel(data);

stats{2, 1} = ' Width';
stats{2, 2} = size(imageData, 2);

stats{3, 1} = ' Height';
stats{3, 2} = size(imageData, 1);

stats{4, 1} = ' Max';
stats{4, 2} = max(data);

stats{5, 1} = ' Min';
stats{5, 2} = min(data);

stats{6, 1} = ' Mean';
stats{6, 2} = mean(data);

stats{7, 1} = ' SD';
stats{7, 2} = std(data);

set(handles.uitable_stats, 'Data', stats);

check_grid(handles);


% -------------------------------------------------------------------
function check_grid(handles)

if 1 == get(handles.checkbox_grid, 'Value')
    grid(handles.axes_hist, 'on');
else
    grid(handles.axes_hist, 'off');
end


% --- Outputs from this function are returned to the command line.
function varargout = image_stats_OutputFcn(hObject, eventdata, handles)
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


% --- Executes when user attempts to close figure_is.
function figure_is_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_is (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainFigHandle = getappdata(handles.figure_is, 'mainFigHandle');
if is_handle(mainFigHandle)
    hFun = getappdata(mainFigHandle, 'hFunCallbackSliceStatsClosed');
    if isa(hFun, 'function_handle')
        feval(hFun, mainFigHandle);
    end
else
    delete(hObject);
end

