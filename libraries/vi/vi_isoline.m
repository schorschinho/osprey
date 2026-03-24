%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright:
% Jun Tan
% University of Texas Southwestern Medical Center
% Department of Radiation Oncology
% Last edited: 08/19/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = vi_isoline(varargin)
% VI_ISOLINE MATLAB code for vi_isoline.fig
%      VI_ISOLINE, by itself, creates a new VI_ISOLINE or raises the existing
%      singleton*.
%
%      H = VI_ISOLINE returns the handle to a new VI_ISOLINE or the handle to
%      the existing singleton*.
%
%      VI_ISOLINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VI_ISOLINE.M with the given input arguments.
%
%      VI_ISOLINE('Property','Value',...) creates a new VI_ISOLINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vi_isoline_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vi_isoline_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vi_isoline

% Last Modified by GUIDE v2.5 25-Aug-2014 11:57:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vi_isoline_OpeningFcn, ...
                   'gui_OutputFcn',  @vi_isoline_OutputFcn, ...
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


% --- Executes just before vi_isoline is made visible.
function vi_isoline_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vi_isoline (see VARARGIN)

% Choose default command line output for vi_isoline
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

numArgs = length(varargin);

assert(1 == numArgs && is_handle(varargin{1}), 'Only allow main GUI handle as argument.');

setappdata(handles.figure_isoline, 'mainFigHandle', varargin{1});

data = cell(20, 2);
data(:, 2) = {false};
set(handles.uitable_isovalues, 'Data', data);

set(hObject, 'Visible', 'off');


% -------------------------------------------------------------------
function send_isovalues_to_main_gui(handles)

mainFigHandle = getappdata(handles.figure_isoline, 'mainFigHandle');
if is_handle(mainFigHandle)
    hFun = getappdata(mainFigHandle, 'hFunCallbackIsolineUpdate');
    if isa(hFun, 'function_handle')
        data = get(handles.uitable_isovalues, 'Data');
        isovalues = [data{:, 1}];
        shown = [data{1:length(isovalues), 2}];
        isovalues = isovalues(shown);
        showCLabel = 1 == get(handles.checkbox_clabel, 'Value');
        feval(hFun, mainFigHandle, isovalues, showCLabel);
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = vi_isoline_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitable_isovalues.
function uitable_isovalues_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_isovalues (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.uitable_isovalues, 'Data');

if 1 == eventdata.Indices(2) % Only need to check data for column 1.
    if isempty(eventdata.EditData) % Delete cell content to delete a value.
        data{eventdata.Indices(1), eventdata.Indices(2)} = [];
    elseif isnan(eventdata.NewData) % Entered a non-numeric string. Revert change.
        data{eventdata.Indices(1), eventdata.Indices(2)} = eventdata.PreviousData;
        disp('Must enter a number.');
        set(handles.uitable_isovalues, 'Data', data);
        return;
    elseif length(find([data{:, 1}] == eventdata.NewData)) > 1 % Entered a existing value. Revert change.
        data{eventdata.Indices(1), eventdata.Indices(2)} = eventdata.PreviousData;
        disp('Must enter a value that did not exist.');
        set(handles.uitable_isovalues, 'Data', data);
        return;
    else % Show by default if entered a new value.
        data(eventdata.Indices(1), 2) = {true};
        set(handles.uitable_isovalues, 'Data', data);
    end
end

data = data(~cellfun(@isempty, data(:, 1)), :);
[~, idx] = sort(cell2mat(data(:, 1)));
data = data(idx, :);
data((size(data, 1)+1):20, 2) = {false};
set(handles.uitable_isovalues, 'Data', data);

send_isovalues_to_main_gui(handles);


% --- Executes when user attempts to close figure_isoline.
function figure_isoline_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_isoline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainFigHandle = getappdata(handles.figure_isoline, 'mainFigHandle');
if is_handle(mainFigHandle)
    hFun = getappdata(mainFigHandle, 'hFunCallbackIsolineClosed');
    if isa(hFun, 'function_handle')
        feval(hFun, mainFigHandle);
    end
else
    delete(hObject);
end


% --- Executes on button press in checkbox_clabel.
function checkbox_clabel_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_clabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_clabel

send_isovalues_to_main_gui(handles);


% --- Executes on button press in pushbutton_showAll.
function pushbutton_showAll_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_showAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.uitable_isovalues, 'Data');
numVals = length([data{:, 1}]);
data(1:numVals, 2) = {true};
set(handles.uitable_isovalues, 'Data', data);
send_isovalues_to_main_gui(handles);


% --- Executes on button press in pushbutton_hideAll.
function pushbutton_hideAll_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_hideAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.uitable_isovalues, 'Data');
numVals = length([data{:, 1}]);
data(1:numVals, 2) = {false};
set(handles.uitable_isovalues, 'Data', data);
send_isovalues_to_main_gui(handles);

