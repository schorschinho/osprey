%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright:
% Jun Tan
% University of Texas Southwestern Medical Center
% Department of Radiation Oncology
% Last edited: 08/19/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = vi(varargin)
% VI MATLAB code for vi.fig
%      VI, by itself, creates a new VI or raises the existing
%      singleton*.
%
%      H = VI returns the handle to a new VI or the handle to
%      the existing singleton*.
%
%      VI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VI.M with the given input arguments.
%
%      VI('Property','Value',...) creates a new VI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vi

% Last Modified by GUIDE v2.5 19-Aug-2014 13:59:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vi_OpeningFcn, ...
                   'gui_OutputFcn',  @vi_OutputFcn, ...
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



% --- Executes just before vi is made visible.
function vi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vi (see VARARGIN)

% Choose default command line output for vi
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vi wait for user response (see UIRESUME)
% uiwait(handles.figure_vi);

if nargin <= 3
    img = flipud(phantom('Modified Shepp-Logan', 512)); % Generate a phantom image if no input is given.
    img = int32(1000 * (img - 0.2) / 0.2);
else
    img = varargin{1};
end

assert(isnumeric(img) || islogical(img), 'First argument must be numeric.');

img = squeeze(img); % Remove singleton dimensions.

maxPixelVal = max(img(:));
minPixelVal = min(img(:));
assert(maxPixelVal ~= minPixelVal, 'Image must contain at least 2 different values.');

numDims = ndims(img);
assert((2 == numDims || 3 == numDims) && all(size(img) > 1), 'Image must be either 2D or 3D.');

options = varargin(2:end);
assert(0 == mod(length(options), 2), 'Argument must be name and value pairs.');

argCheckMap = containers.Map;
argCheckMap('clim') = false;
argCheckMap('aspect') = false;
setappdata(handles.figure_vi, 'argCheckMap', argCheckMap);
while ~isempty(options) && length(options) >= 2
    check_arg(handles, options{1}, options{2});
    options = options(3:end);
end
argCheckMap = getappdata(handles.figure_vi, 'argCheckMap');

imgSize = size(img);

if 2 == numDims % 2D images don't need these uicontrols.
    set([ ...
        handles.radiobutton_sagittalView, ...
        handles.radiobutton_coronalView, ...
        handles.radiobutton_3dSlice, ...
        handles.togglebutton_light, ...
        handles.togglebutton_rotate, ...
        handles.slider_xSliceNo, ...
        handles.slider_ySliceNo, ...
        handles.slider_zSliceNo, ...
        handles.edit_sliceNo, ...
        handles.slider_sliceNo], ...
        'Enable', 'off');
    set(handles.slider_sliceNo, 'Max', 1, 'Min', 0, 'SliderStep', [0.1 0.1]);
elseif 3 == numDims % 3D images need multiple views, rotation, etc.
    set([ ...
        handles.radiobutton_sagittalView, ...
        handles.radiobutton_coronalView, ...
        handles.radiobutton_3dSlice, ...
        handles.togglebutton_light, ...
        handles.togglebutton_rotate, ...
        handles.slider_xSliceNo, ...
        handles.slider_ySliceNo, ...
        handles.slider_zSliceNo, ...
        handles.edit_sliceNo, ...
        handles.slider_sliceNo], ...
        'Enable', 'on');
    set(handles.slider_xSliceNo, 'Max', imgSize(2), 'Min', 1, 'SliderStep', [1/(imgSize(2)-1) 10/(imgSize(2)-1)]);
    set(handles.slider_ySliceNo, 'Max', imgSize(1), 'Min', 1, 'SliderStep', [1/(imgSize(1)-1) 10/(imgSize(1)-1)]);
    set(handles.slider_zSliceNo, 'Max', imgSize(3), 'Min', 1, 'SliderStep', [1/(imgSize(3)-1) 10/(imgSize(3)-1)]);
    sliceNo3d = ceil(imgSize / 2); % Initially display center slices.
    sliceNo3d = sliceNo3d([2 1 3]); % Reorder as [x y z].
    setappdata(handles.figure_vi, 'sliceNo3d', sliceNo3d);
    set(handles.slider_xSliceNo, 'Value', sliceNo3d(1));
    set(handles.slider_ySliceNo, 'Value', sliceNo3d(2));
    set(handles.slider_zSliceNo, 'Value', sliceNo3d(3));
    set(handles.slider_sliceNo, 'Min', 1);
end

setappdata(handles.figure_vi, 'imgData', img);
setappdata(handles.figure_vi, 'imgSize', imgSize);
setappdata(handles.figure_vi, 'numDims', numDims);

cb = colorbar('peer', handles.axes_colorBar, 'location', 'west');
setappdata(handles.figure_vi, 'colorBar', cb);

% Get number string format.
if isinteger(img) || islogical(img)
    dataFormat = '%.0f';
else % if isfloat(img)
    maxVal = max(abs(maxPixelVal), abs(minPixelVal));
    if  maxVal > 1e5 || maxVal < 1e-1
        dataFormat = '%.3e';
    else
        dataFormat = '%.3f';
    end
end
setappdata(handles.figure_vi, 'dataFormat', dataFormat);
set(handles.text_maxPixelVal, 'String', num2str(maxPixelVal, dataFormat));
set(handles.text_minPixelVal, 'String', num2str(minPixelVal, dataFormat));

set(handles.axes_2dViewer, 'CLimMode', 'manual');
if ~argCheckMap('clim')
    set_clims(handles, minPixelVal, maxPixelVal);
end

if ~argCheckMap('aspect')
    setappdata(handles.figure_vi, 'aspectRatio', [1 1 1]);
end

setappdata(handles.figure_vi, 'buttonDown', false);
setappdata(handles.figure_vi, 'ptOnAxesBtnDown', [1 1]);

setappdata(handles.figure_vi, 'regionType', 'Rectangle');

setup_ui_menus(handles.figure_vi);

update_window_info(handles);

plot_resize_arrow(handles.axes_resizeArrow);

update_view_type(handles);


% -------------------------------------------------------------------
function setup_ui_menus(hFig)

uiMenus = struct;

guiChildren = get(hFig, 'Children');
for i = 1 : length(guiChildren)
    child = guiChildren(i);
    childTag = get(child, 'Tag');
    if strcmpi('menu_moreWindowSettings', childTag)
        uiMenus.menu_moreWindowSettings = child;
    elseif strcmpi('menu_selectRegionType', childTag)
        uiMenus.menu_selectRegionType = child;
    end
end

setappdata(hFig, 'uiMenus', uiMenus);


% -------------------------------------------------------------------
function plot_resize_arrow(hAxes)

hold(hAxes, 'on');
setappdata(hAxes, 'resizeArrow', [ ...
    plot(hAxes, [9 15], [9 15], 'b'), ...
    plot(hAxes, [15 15], [15 1], 'b'), ...
    plot(hAxes, [15 1], [15 15], 'b')]);
set(getappdata(hAxes, 'resizeArrow'), 'HitTest', 'off');
hold(hAxes, 'off');
axis(hAxes, 'off');


% -------------------------------------------------------------------
function check_arg(handles, argName, argVal)

assert(ischar(argName), 'Argument name must be characters.');

argCheckMap = getappdata(handles.figure_vi, 'argCheckMap');

if strcmpi('window', argName) && isnumeric(argVal) && 2 == length(argVal)
    set_clims(handles, argVal(2) - argVal(1) / 2, argVal(2) + argVal(1) / 2);
    argCheckMap('clim') = true;
elseif strcmpi('range', argName) && isnumeric(argVal) && 2 == length(argVal)
    set_clims(handles, argVal(1), argVal(2));
    argCheckMap('clim') = true;
elseif strcmpi('aspect', argName) && isnumeric(argVal) && 3 == length(argVal)
    setappdata(handles.figure_vi, 'aspectRatio', argVal);
    argCheckMap('aspect') = true;
else
    error([argName ' is not a valid argument name.']);
end

setappdata(handles.figure_vi, 'argCheckMap', argCheckMap);


% -------------------------------------------------------------------
function update_window_info(handles)

dataFormat = getappdata(handles.figure_vi, 'dataFormat');

cLim = get(handles.axes_2dViewer, 'CLim');

set(handles.edit_windowMax, 'String', num2str(cLim(2), dataFormat));
set(handles.edit_windowMin, 'String', num2str(cLim(1), dataFormat));
set(handles.edit_windowWidth, 'String', num2str(range(cLim), dataFormat));
set(handles.edit_windowLevel, 'String', num2str(mean(cLim), dataFormat));


% -------------------------------------------------------------------
function update_view_type(handles)

sliceNo3d = getappdata(handles.figure_vi, 'sliceNo3d');

viewType = get_view_type(handles);

if viewType == 4
    
    set(handles.uipanel_2dViewer, 'Visible', 'off');
    set(handles.uipanel_3dSlicer, 'Visible', 'on');
    set([ ...
        handles.togglebutton_light, ...
        handles.togglebutton_rotate, ...
        handles.slider_xSliceNo, ...
        handles.slider_ySliceNo, ...
        handles.slider_zSliceNo], ...
        'Enable', 'on');
    
    update_slice_3d(handles, sliceNo3d(1), sliceNo3d(2), sliceNo3d(3), true);
  
else
    
    rotate3d(handles.axes_3dSlicer, 'off'); % If in 3D rotating status, turn it off.
    
    set(handles.uipanel_2dViewer, 'Visible', 'on');
    set(handles.uipanel_3dSlicer, 'Visible', 'off');
    set([ ...
        handles.togglebutton_light, ...
        handles.togglebutton_rotate, ...
        handles.slider_xSliceNo, ...
        handles.slider_ySliceNo, ...
        handles.slider_zSliceNo], ...
        'Enable', 'off');
    
    imgSize = getappdata(handles.figure_vi, 'imgSize');
    
    switch viewType
        case 1
            sliceSize = imgSize([1 2]);
            if 2 == getappdata(handles.figure_vi, 'numDims')
                numSlices = 1;
                sliceNo = 1;
            else
                numSlices = imgSize(3);
                sliceNo = sliceNo3d(3);
            end
        case 2
            sliceSize = imgSize([3 1]);
            numSlices = imgSize(2);
            sliceNo = sliceNo3d(1);
        case 3
            sliceSize = imgSize([3 2]);
            numSlices = imgSize(1);
            sliceNo = sliceNo3d(2);
        otherwise
            error('Invalid view type!');
    end
    
    setappdata(handles.figure_vi, 'sliceSize', sliceSize);
    setappdata(handles.figure_vi, 'numSlices', numSlices);
    set(handles.text_numSlices, 'String', ['/ ' num2str(numSlices)]);
    
    if numSlices > 1
        set(handles.slider_sliceNo, ...
            'Max', numSlices, ...
            'SliderStep', [1/(numSlices-1) 10/(numSlices-1)]);
    end
    
    update_slice_2d(handles, sliceNo, true);
    
end


% -------------------------------------------------------------------
function update_slice_3d(handles, xSliceNo, ySliceNo, zSliceNo, newView)

oldSliceNo3d = getappdata(handles.figure_vi, 'sliceNo3d');

imgSize = getappdata(handles.figure_vi, 'imgSize');

xSliceNo = max(min(round(xSliceNo), imgSize(2)), 1);
ySliceNo = max(min(round(ySliceNo), imgSize(1)), 1);
zSliceNo = max(min(round(zSliceNo), imgSize(3)), 1);

setappdata(handles.figure_vi, 'sliceNo3d', [xSliceNo ySliceNo zSliceNo]);

set(handles.slider_xSliceNo, 'Value', xSliceNo);
set(handles.slider_ySliceNo, 'Value', ySliceNo);
set(handles.slider_zSliceNo, 'Value', zSliceNo);

img = getappdata(handles.figure_vi, 'imgData');

hSlices = getappdata(handles.figure_vi, 'hSlices');

if newView && isempty(hSlices)
    
    hSlices = slice(handles.axes_3dSlicer, double(img), xSliceNo, ySliceNo, zSliceNo);
    setappdata(handles.figure_vi, 'hSlices', hSlices);

    xlabel(handles.axes_3dSlicer, 'x - L/R');
    ylabel(handles.axes_3dSlicer, 'y - A/P');
    zlabel(handles.axes_3dSlicer, 'z - I/S');

    set(handles.axes_3dSlicer, ...
        'XLim', [1 imgSize(2)], 'YLim', [1 imgSize(1)], 'ZLim', [1 imgSize(3)]);

    cLim = get(handles.axes_2dViewer, 'CLim');
    set_clims(handles, cLim(1), cLim(2));
    
    shading(handles.axes_3dSlicer, 'flat');
    
    aspectRatio = getappdata(handles.figure_vi, 'aspectRatio');
    daspect(handles.axes_3dSlicer, 1 ./ [aspectRatio(1) aspectRatio(2) aspectRatio(3)]);
    
    set_colormap(handles);

    % Create light source at camera position.
    hLight = light('Parent', handles.axes_3dSlicer, ...
        'Position', get(handles.axes_3dSlicer, 'CameraPosition'));
    setappdata(handles.figure_vi, 'hLight', hLight);
    
    set_light(handles);

else
    
    if oldSliceNo3d(1) ~= xSliceNo
        set(hSlices(1), ...
            'CData', double(squeeze(img(:, xSliceNo, :))), ...
            'XData', ones(imgSize(1), imgSize(3)) * xSliceNo);
    end
    
    if oldSliceNo3d(2) ~= ySliceNo
        set(hSlices(2), ...
            'CData', double(squeeze(img(ySliceNo, :, :))), ...
            'YData', ones(imgSize(2), imgSize(3)) * ySliceNo);
    end
    
    if oldSliceNo3d(3) ~= zSliceNo
        set(hSlices(3), ...
            'CData', double(squeeze(img(:, :, zSliceNo))), ...
            'ZData', ones(imgSize(1), imgSize(2)) * zSliceNo);
    end
    
end

check_rotate(handles);


% -------------------------------------------------------------------
function set_light(handles)

axes(handles.axes_3dSlicer); % Must make axes_3dSlicer current for lighting.

if 1 == get(handles.togglebutton_light, 'Value')
    lighting phong;
else
    lighting none;
end


% -------------------------------------------------------------------
function rotate_pre_callback(obj, evd)

setappdata(obj, 'rotating', true);


% -------------------------------------------------------------------
function rotate_post_callback(obj, evd)

setappdata(obj, 'rotating', false);


% -------------------------------------------------------------------
function check_rotate(handles)

if 1 == get(handles.togglebutton_rotate, 'Value')
    rotate3d(handles.axes_3dSlicer, 'on');
    h = rotate3d(handles.figure_vi);
    set(h, ...
        'ActionPreCallback', @rotate_pre_callback, ...
        'ActionPostCallback', @rotate_post_callback, ...
        'Enable', 'on');
else
    rotate3d(handles.axes_3dSlicer, 'off');
end


% -------------------------------------------------------------------
function update_slice_2d(handles, sliceNo, newView)

numSlices = getappdata(handles.figure_vi, 'numSlices');
sliceNo = max(min(round(sliceNo), numSlices), 1);

sliceNoEdit = round(str2double(get(handles.edit_sliceNo, 'String')));
if sliceNo ~= sliceNoEdit
    set(handles.edit_sliceNo, 'String', num2str(sliceNo));
end

sliceNoSlider = round(get(handles.slider_sliceNo, 'Value'));
if sliceNo ~= sliceNoSlider
    set(handles.slider_sliceNo, 'Value', sliceNo);
end

setappdata(handles.figure_vi, 'sliceNo', sliceNo);
set(handles.edit_sliceNo, 'String', num2str(sliceNo));

img = getappdata(handles.figure_vi, 'imgData');
sliceNo3d = getappdata(handles.figure_vi, 'sliceNo3d');
viewType = get_view_type(handles);
switch viewType
    case 1
        sliceImage = squeeze(img(:, :, sliceNo));
        sliceNo3d(3) = sliceNo;
    case 2
        sliceImage = squeeze(img(:, sliceNo, :))';
        sliceNo3d(1) = sliceNo;
    case 3
        sliceImage = squeeze(img(sliceNo, :, :))';
        sliceNo3d(2) = sliceNo;
    otherwise
        error('Invalid view type!');
end
setappdata(handles.figure_vi, 'sliceNo3d', sliceNo3d);

setappdata(handles.figure_vi, 'sliceImage', sliceImage);

if newView
    cLim = get(handles.axes_2dViewer, 'CLim');
    hSliceImage = imshow(flipud(sliceImage), cLim, 'Parent', handles.axes_2dViewer);
    set_colormap(handles);
    setappdata(handles.figure_vi, 'hSliceImage', hSliceImage);
else
    hSliceImage = getappdata(handles.figure_vi, 'hSliceImage');
    set(hSliceImage, 'CData', flipud(sliceImage));
end

update_aspect_ratio(handles);

hFigSliceStats = getappdata(handles.figure_vi, 'hFigSliceStats');
if is_handle(hFigSliceStats)
    feval(getappdata(hFigSliceStats, 'entry_update_data'), ...
        hFigSliceStats, {sliceImage, handles.figure_vi});
end

hDrawRegion = getappdata(handles.figure_vi, 'hDrawRegion');
hFigRegionMeasurement = getappdata(handles.figure_vi, 'hFigRegionMeasurement');
if is_handle(hFigRegionMeasurement) && is_handle(hDrawRegion)
    if strcmpi('Drawn', getappdata(handles.figure_vi, 'regionType'))
        hDrawRegionAdd = getappdata(handles.figure_vi, 'hDrawRegionAdd');
        hFigRegionMeasurement = getappdata(handles.figure_vi, 'hFigRegionMeasurement');
        if is_handle(hFigRegionMeasurement) && is_handle(hDrawRegionAdd)
            update_user_drawn_region_data(handles);
            hDrawRegionAdd = getappdata(handles.figure_vi, 'hDrawRegionAdd');
            if is_handle(hDrawRegionAdd) && strcmpi('on', get(hDrawRegionAdd, 'Visible'))
                set(hDrawRegionAdd, 'LineStyle', '-', 'Color', 'r');
            end
        end
    else
        update_region_data(handles);
    end
end


hDrawLine = getappdata(handles.figure_vi, 'hDrawLine');
hFigLineMeasurement = getappdata(handles.figure_vi, 'hFigLineMeasurement');
if is_handle(hFigLineMeasurement) && is_handle(hDrawLine)
    update_line_data(handles);
end

if 1 == get(handles.togglebutton_isoline, 'Value')
    update_isolines(handles);
end


% -------------------------------------------------------------------
function set_clims(handles, loLim, hiLim)

set(handles.axes_2dViewer, 'CLim', [loLim hiLim]);
set(handles.axes_3dSlicer, 'CLim', [loLim hiLim]);
set(handles.axes_colorBar, 'CLim', [loLim hiLim]);


% -------------------------------------------------------------------
function set_colormap(handles)

cm = getappdata(handles.figure_vi, 'pixelCMap');

if isempty(cm)
    setappdata(handles.figure_vi, 'pixelCMap', colormap(handles.axes_2dViewer));
else
    cmInd = get(handles.listbox_colorMaps,'Value');
    switch cmInd
        case 1
            cm = getappdata(handles.figure_vi, 'pixelCMap');
            colormap(handles.axes_2dViewer, cm);
        otherwise
            contents = cellstr(get(handles.listbox_colorMaps,'String'));
            cmapFun = lower(contents{cmInd});
            colormap(handles.axes_2dViewer, feval(cmapFun, 256));
    end
end


% -------------------------------------------------------------------
function update_aspect_ratio(handles)

aspectRatio = getappdata(handles.figure_vi, 'aspectRatio');

switch get_view_type(handles)
    case 1
        set(handles.axes_2dViewer, 'DataAspectRatio', aspectRatio([1 2 3]));
    case 2
        set(handles.axes_2dViewer, 'DataAspectRatio', aspectRatio([3 2 1]));
    case 3
        set(handles.axes_2dViewer, 'DataAspectRatio', aspectRatio([3 1 2]));
    otherwise
end

daspect(handles.axes_3dSlicer, 1 ./ [aspectRatio(1) aspectRatio(2) aspectRatio(3)]);

set(handles.edit_xAspectRatio, 'String', num2str(aspectRatio(1), '%.1f'));
set(handles.edit_yAspectRatio, 'String', num2str(aspectRatio(2), '%.1f'));
set(handles.edit_zAspectRatio, 'String', num2str(aspectRatio(3), '%.1f'));


% -------------------------------------------------------------------
function viewType = get_view_type(handles)

viewType = find(get(handles.uipanel_viewType, 'SelectedObject') == ...
    [handles.radiobutton_transverseView ...
    handles.radiobutton_sagittalView ...
    handles.radiobutton_coronalView ...
    handles.radiobutton_3dSlice]);


% -------------------------------------------------------------------
function set_obj_pos(hObj, x, top, width, height)

parentPos = get(get(hObj, 'Parent'), 'Position');
y = parentPos(4) - top - height + 1;
set(hObj, 'Position', [x y width height]);


% -------------------------------------------------------------------
function [minWidth, minHeight] = get_min_gui_size(hFig)

uicontrolTypes = ...
    {'checkbox', 'edit', 'listbox', 'pushbutton', ....
    'radiobutton', 'slider', 'text', 'togglebutton', ...
    'axes', 'uitable', 'uipanel', 'uibuttongroup'};

children = get(hFig, 'Children');
children(~ismember(get(children, 'Type'), uicontrolTypes)) = [];
children(strcmpi(get(children, 'Visible'), 'off')) = [];
children(strcmpi(get(children, 'Tag'), 'uipanel_resizeArrow')) = [];

numChildren = length(children);
width = zeros(numChildren, 1);
top = zeros(numChildren, 1);
bottom = zeros(numChildren, 1) + 999;

for i = 1 : numChildren
    child = children(i);
    pos = get(child, 'Position');
    width(i) = pos(1) + pos(3);
    top(i) = pos(2) + pos(4);
    bottom(i) = pos(2);
    
    if isprop(child, 'TightInset')
        ti = get(child, 'TightInset');
        width(i) = width(i) + ti(3);
        top(i) = top(i) + ti(4);
        bottom(i) = bottom(i) - ti(2);
    end
end

minWidth = max(width);
minHeight = max(top) - min(bottom);


% -------------------------------------------------------------------
function auto_layout(handles)

guiPos = get(handles.figure_vi, 'Position');

viewSize = getappdata(handles.figure_vi, 'sliceSize');

set(handles.axes_top, 'XLim', [1 viewSize(2)]);
set(handles.axes_bottom, 'XLim', [1 viewSize(2)]);
set(handles.axes_left, 'YLim', [1 viewSize(1)]);
set(handles.axes_right, 'YLim', [1 viewSize(1)]);

aspectRatio = getappdata(handles.figure_vi, 'aspectRatio');

viewType = get_view_type(handles);

switch viewType
    case 1
        viewSize(1) = viewSize(1) * aspectRatio(1);
        viewSize(2) = viewSize(2) * aspectRatio(2);
    case 2
        viewSize(1) = viewSize(1) * aspectRatio(3);
        viewSize(2) = viewSize(2) * aspectRatio(1);
    case 3
        viewSize(1) = viewSize(1) * aspectRatio(3);
        viewSize(2) = viewSize(2) * aspectRatio(2);
    otherwise
end

leftPanelPos = get(handles.uipanel_leftControls, 'Position');

controlPanelPos2d = get(handles.uipanel_2dControls, 'Position');
tiT = get(handles.axes_top, 'TightInset');
tiTHeight = tiT(4) + 10 + tiT(2);
tiB = get(handles.axes_bottom, 'TightInset');
tiBHeight = tiB(4) + 10 + tiB(2);
tiL = get(handles.axes_left, 'TightInset');
tiLWidth = tiL(1) + 10 + tiL(3);
tiR = get(handles.axes_right, 'TightInset');
tiRWidth = tiR(1) + 10 + tiR(3);
viewerX = leftPanelPos(3) + 5;
viewerY = 1;
viewerWidth = max(controlPanelPos2d(3), tiLWidth + 2 + viewSize(2) + 2 + tiRWidth);
viewerHeight = controlPanelPos2d(4) + 2 + tiTHeight + 2 + viewSize(1) + 2 + tiBHeight;

controlPanelPos3d = get(handles.uipanel_3dControls, 'Position');
imgSize = getappdata(handles.figure_vi, 'imgSize');
if 2 == length(imgSize)
    imgSize(3) = 1;
end
slicerX = leftPanelPos(3) + 5;
slicerY = 1;
slicerWidth = ceil(max(controlPanelPos3d(3), max(imgSize .* aspectRatio) * 1.4));
slicerHeight = controlPanelPos3d(4) + slicerWidth;

set_obj_pos(handles.uipanel_leftControls, leftPanelPos(1), 1, leftPanelPos(3), leftPanelPos(4));

set_obj_pos(handles.uipanel_2dViewer, viewerX, viewerY, viewerWidth, viewerHeight);
set_obj_pos(handles.uipanel_2dControls, 1, 1, controlPanelPos2d(3), controlPanelPos2d(4));
set_obj_pos(handles.axes_top, ...
    tiLWidth+2, ...
    controlPanelPos2d(4) + 2 + tiT(4)+2, ...
    viewSize(2), 10);
set_obj_pos(handles.axes_bottom, ...
    tiLWidth+2, ...
    controlPanelPos2d(4) + 2 + tiTHeight + 2 + viewSize(1) + 2, ...
    viewSize(2), 10);
set_obj_pos(handles.axes_left, ...
    tiL(1) + 2, ...
    controlPanelPos2d(4) + 2 + tiTHeight + 2, ...
    10, viewSize(1));
set_obj_pos(handles.axes_right, ...
    tiLWidth + 2 + viewSize(2) + 2, ...
    controlPanelPos2d(4) + 2 + tiTHeight + 2, ...
    10, viewSize(1));
set_obj_pos(handles.axes_2dViewer, ...
    tiLWidth + 2, ...
    controlPanelPos2d(4) + 2 + tiTHeight + 2, ...
    viewSize(2), viewSize(1));

set_obj_pos(handles.uipanel_3dSlicer, slicerX, slicerY, slicerWidth, slicerHeight);
set_obj_pos(handles.uipanel_3dControls, 1, 1, controlPanelPos3d(3), controlPanelPos3d(4));
set_obj_pos(handles.axes_3dSlicer, 70, 80, slicerWidth - 120, slicerHeight - 130);

cb = getappdata(handles.figure_vi, 'colorBar');
if 4 == viewType
    colorBarX = slicerX + slicerWidth + 10;
else
    colorBarX = viewerX + viewerWidth + 10;
end
colorBarY = max(controlPanelPos2d(4), controlPanelPos3d(4)) + 30;
colorBarHeight = max(256, min(512, guiPos(4) - colorBarY - 30));
set_obj_pos(handles.axes_colorBar, colorBarX, colorBarY, 20, colorBarHeight);
set(cb, 'Units', 'pixels');
set_obj_pos(cb, colorBarX, colorBarY, 20, colorBarHeight);

[minWidth, minHeight] = get_min_gui_size(handles.figure_vi);

if minWidth > guiPos(3) || minHeight > guiPos(4)
    set(handles.uipanel_resizeArrow, ...
        'Position', [guiPos(3)-20 5 17 17], ...
        'Visible', 'on');
else
    set(handles.uipanel_resizeArrow, 'Visible', 'off');
end


% -------------------------------------------------------------------
function [ptOnAxesCurrent, onImage] = get_pointer_pos(handles)

if 4 == get_view_type(handles)
    imgSize = getappdata(handles.figure_vi, 'imgSize');
    ptOnAxesCurrent = round(get(handles.axes_3dSlicer, 'CurrentPoint'));
    
    onImage = all(ptOnAxesCurrent(:) >= 1) ...
        && ptOnAxesCurrent(1, 1) <= imgSize(2) ...
        && ptOnAxesCurrent(1, 2) <= imgSize(1) ...
        && ptOnAxesCurrent(1, 3) <= imgSize(3) ...
        && ptOnAxesCurrent(2, 1) <= imgSize(2) ...
        && ptOnAxesCurrent(2, 2) <= imgSize(1) ...
        && ptOnAxesCurrent(2, 3) <= imgSize(3);
    
else
    sliceSize = getappdata(handles.figure_vi, 'sliceSize');
    ptOnAxesCurrent = get(handles.axes_2dViewer, 'CurrentPoint');
    ptOnAxesCurrent = round(ptOnAxesCurrent(1, 1:2));
    ptOnAxesCurrent(2) = sliceSize(1) - ptOnAxesCurrent(2) + 1;
    
    onImage = ptOnAxesCurrent(1) >= 1 && ptOnAxesCurrent(1) <= sliceSize(2) ...
        && ptOnAxesCurrent(2) >= 1 && ptOnAxesCurrent(2) <= sliceSize(1);
    
end


% -------------------------------------------------------------------
function event_image_stats_closed(hGui)

assert(is_handle(hGui));
handles = guidata(hGui);
set(handles.togglebutton_sliceStats, 'Value', 0);
set_image_stats(handles, 0);


% -------------------------------------------------------------------
function event_region_measurement_closed(hGui)

assert(is_handle(hGui));
handles = guidata(hGui);
set(handles.togglebutton_regionMeasure, 'Value', 0);
set_region_measure(handles, 0);


% -------------------------------------------------------------------
function event_line_measurement_closed(hGui)

assert(is_handle(hGui));
handles = guidata(hGui);
set(handles.togglebutton_lineMeasure, 'Value', 0);
set_line_measure(handles, 0);


% -------------------------------------------------------------------
function event_isoline_closed(hGui)

assert(is_handle(hGui));
handles = guidata(hGui);
set(handles.togglebutton_isoline, 'Value', 0);
set_isoline(handles, 0);


% -------------------------------------------------------------------
function msg_map_isoline_update(hGui, isovalues, showCLabel)

assert(is_handle(hGui));
handles = guidata(hGui);

setappdata(handles.figure_vi, 'isovalues', isovalues);
setappdata(handles.figure_vi, 'showCLabel', showCLabel);

update_isolines(handles);


% -------------------------------------------------------------------
function update_isolines(handles)

if 0 == get(handles.togglebutton_isoline, 'Value')
    return;
end

hIsolines = getappdata(handles.figure_vi, 'hIsolines');
if is_handle(hIsolines)
    delete(hIsolines);
end

isovalues = getappdata(handles.figure_vi, 'isovalues');
if 1 == numel(isovalues)
    isovalues = [isovalues isovalues];
end

sliceImage = getappdata(handles.figure_vi, 'sliceImage');

hold(handles.axes_2dViewer, 'on');

[c, hIsolines] = contour(handles.axes_2dViewer, flipud(sliceImage), isovalues, 'color', 'r');

showCLabel = getappdata(handles.figure_vi, 'showCLabel');
if showCLabel
    clabel(c, hIsolines);
end

hold(handles.axes_2dViewer, 'off');

setappdata(handles.figure_vi, 'hIsolines', hIsolines);


% -------------------------------------------------------------------
function set_image_stats(handles, showTool)

hFigSliceStats = getappdata(handles.figure_vi, 'hFigSliceStats');

if 1 == showTool
    if ~is_handle(hFigSliceStats)
        sliceImage = getappdata(handles.figure_vi, 'sliceImage');
        hFigSliceStats = image_stats(sliceImage, handles.figure_vi);
        setappdata(handles.figure_vi, 'hFigSliceStats', hFigSliceStats);
        setappdata(handles.figure_vi, 'hFunCallbackSliceStatsClosed', @event_image_stats_closed);

        pos1 = get(handles.figure_vi, 'Position');
        pos2 = get(hFigSliceStats, 'Position');
        set(hFigSliceStats, 'Position', [pos1(1)+pos1(3)+16 pos1(2) pos2(3) pos2(4)]);
    end
    
    set(hFigSliceStats, 'Visible', 'on');
    figure(handles.figure_vi);
else
    if is_handle(hFigSliceStats)
        set(hFigSliceStats, 'Visible', 'off');
    end
end


% -------------------------------------------------------------------
function set_region_measure(handles, showTool)

% Delete region plot on image.
hDrawRegion = getappdata(handles.figure_vi, 'hDrawRegion');
if is_handle(hDrawRegion)
    delete(hDrawRegion);
end

hDrawRegionAdd = getappdata(handles.figure_vi, 'hDrawRegionAdd');
if is_handle(hDrawRegionAdd)
    delete(hDrawRegionAdd);
end

hFigRegionMeasurement = getappdata(handles.figure_vi, 'hFigRegionMeasurement');

if 1 == showTool
    if ~is_handle(hFigRegionMeasurement)
        hFigRegionMeasurement = region_measurement(handles.figure_vi);
        setappdata(handles.figure_vi, 'hFigRegionMeasurement', hFigRegionMeasurement);
        setappdata(handles.figure_vi, 'hFunCallbackRegionMeasurementClosed', @event_region_measurement_closed);
    
        pos1 = get(handles.figure_vi, 'Position');
        pos2 = get(hFigRegionMeasurement, 'Position');
        set(hFigRegionMeasurement, 'Position', [pos1(1)+pos1(3)+16 pos1(2) pos2(3) pos2(4)]);
    end
    
    set(hFigRegionMeasurement, 'Visible', 'on');
    figure(handles.figure_vi);

    if strcmpi('Drawn', getappdata(handles.figure_vi, 'regionType'))
        update_user_drawn_region_data(handles);
        hDrawRegionAdd = getappdata(handles.figure_vi, 'hDrawRegionAdd');
        if is_handle(hDrawRegionAdd) && strcmpi('on', get(hDrawRegionAdd, 'Visible'))
            set(hDrawRegionAdd, 'LineStyle', '-', 'Color', 'r');
        end
    else
        update_region_data(handles);
    end
else
    if is_handle(hFigRegionMeasurement)
        set(hFigRegionMeasurement, 'Visible', 'off');
    end
end


% -------------------------------------------------------------------
function set_line_measure(handles, showTool)

% Delete line plot on image.
hDrawLine = getappdata(handles.figure_vi, 'hDrawLine');
if is_handle(hDrawLine)
    delete(hDrawLine);
end

hFigLineMeasurement = getappdata(handles.figure_vi, 'hFigLineMeasurement');

if 1 == showTool
    if ~is_handle(hFigLineMeasurement) % Open a new line measurement GUI next to main GUI.
        hFigLineMeasurement = line_measurement(handles.figure_vi);
        setappdata(handles.figure_vi, 'hFigLineMeasurement', hFigLineMeasurement);
        setappdata(handles.figure_vi, 'hFunCallbackLineMeasurementClosed', @event_line_measurement_closed);

        pos1 = get(handles.figure_vi, 'Position');
        pos2 = get(hFigLineMeasurement, 'Position');
        set(hFigLineMeasurement, 'Position', [pos1(1)+pos1(3)+16 pos1(2) pos2(3) pos2(4)]);
    end
    
    set(hFigLineMeasurement, 'Visible', 'on');
    figure(handles.figure_vi);
    
    update_line_data(handles);
else
    if is_handle(hFigLineMeasurement)
        set(hFigLineMeasurement, 'Visible', 'off');
    end
end


% -------------------------------------------------------------------
function set_isoline(handles, showTool)

hIsolines = getappdata(handles.figure_vi, 'hIsolines');
if is_handle(hIsolines)
    delete(hIsolines);
end

hFigIsoline = getappdata(handles.figure_vi, 'hFigIsoline');

if 1 == showTool
    if ~is_handle(hFigIsoline) % Open a new isoline GUI next to main GUI.
        hFigIsoline = vi_isoline(handles.figure_vi);
        setappdata(handles.figure_vi, 'hFigIsoline', hFigIsoline);
        setappdata(handles.figure_vi, 'hFunCallbackIsolineClosed', @event_isoline_closed);
        setappdata(handles.figure_vi, 'hFunCallbackIsolineUpdate', @msg_map_isoline_update);

        pos1 = get(handles.figure_vi, 'Position');
        pos2 = get(hFigIsoline, 'Position');
        set(hFigIsoline, 'Position', [pos1(1)+pos1(3)+16 pos1(2) pos2(3) pos2(4)]);
    end
    
    set(hFigIsoline, 'Visible', 'on');
    figure(handles.figure_vi);
    
    update_isolines(handles);
else
    if is_handle(hFigIsoline)
        set(hFigIsoline, 'Visible', 'off');
    end
end


% -------------------------------------------------------------------
function update_user_drawn_region_data(handles)

userDrawRegionVertices = getappdata(handles.figure_vi, 'userDrawRegionVertices');
if isempty(userDrawRegionVertices)
    return;
end

hDrawRegion = getappdata(handles.figure_vi, 'hDrawRegion');
if ishandle(hDrawRegion)
    delete(hDrawRegion);
    setappdata(handles.figure_vi, 'hDrawRegion', []);
end

hDrawRegionAdd = getappdata(handles.figure_vi, 'hDrawRegionAdd');
if ishandle(hDrawRegionAdd)
    delete(hDrawRegionAdd);
    setappdata(handles.figure_vi, 'hDrawRegionAdd', []);
end

viewSize = getappdata(handles.figure_vi, 'sliceSize');
x = userDrawRegionVertices(:, 1);
y = userDrawRegionVertices(:, 2);

sliceImage = getappdata(handles.figure_vi, 'sliceImage');

hold(handles.axes_2dViewer, 'on');
hDrawRegion = plot(handles.axes_2dViewer, x,  viewSize(1) - y, 'r');
hDrawRegionAdd = plot(handles.axes_2dViewer, x([end 1]),  viewSize(1) - y([end 1]), ':y');
hold(handles.axes_2dViewer, 'off');
setappdata(handles.figure_vi, 'hDrawRegion', hDrawRegion);
setappdata(handles.figure_vi, 'hDrawRegionAdd', hDrawRegionAdd);

hFigRegionMeasurement = getappdata(handles.figure_vi, 'hFigRegionMeasurement');
if is_handle(hFigRegionMeasurement)
    mask = poly2mask(x, y, viewSize(1), viewSize(2));
    feval(getappdata(hFigRegionMeasurement, 'entry_update_data'), ...
        hFigRegionMeasurement, ...
        {{double(sliceImage), mask}, getappdata(handles.figure_vi, 'regionType'), handles.figure_vi});
end


% -------------------------------------------------------------------
function update_region_data(handles)

posRegion = getappdata(handles.figure_vi, 'posRegion');
if isempty(posRegion)
    return;
end

x1 = posRegion(1);
x2 = posRegion(2);
y1 = posRegion(3);
y2 = posRegion(4);

hDrawRegion = getappdata(handles.figure_vi, 'hDrawRegion');
viewSize = getappdata(handles.figure_vi, 'sliceSize');
rectX =  min(x1, x2);
rectY =  min(viewSize(1) - [y1 y2]) + 1;
rectWidth = range([x1, x2]);
rectHeight = range([y1, y2]);
if rectWidth > 0 && rectHeight > 0
    if is_handle(hDrawRegion)
        set(hDrawRegion, 'Position', [rectX, rectY, rectWidth, rectHeight]);
    else
        if strcmpi('Rectangle', getappdata(handles.figure_vi, 'regionType'))
            curvature = [0 0];
        else
            curvature = [1 1];
        end
        
        hDrawRegion = rectangle( ...
            'Position', [rectX, rectY, rectWidth, rectHeight], ...
            'Parent', handles.axes_2dViewer, ...
            'EdgeColor', 'r', ...
            'Curvature', curvature);
        
        setappdata(handles.figure_vi, 'hDrawRegion', hDrawRegion);
    end
    
end

sliceImage = getappdata(handles.figure_vi, 'sliceImage');

hFigRegionMeasurement = getappdata(handles.figure_vi, 'hFigRegionMeasurement');
if is_handle(hFigRegionMeasurement)
    feval(getappdata(hFigRegionMeasurement, 'entry_update_data'), ...
        hFigRegionMeasurement, ...
        {sliceImage(min(y1, y2) : max(y1, y2), min(x1, x2) : max(x1, x2)), ...
        getappdata(handles.figure_vi, 'regionType'), handles.figure_vi});
end


% -------------------------------------------------------------------
function update_line_data(handles)

posLine = getappdata(handles.figure_vi, 'posLine');
if isempty(posLine)
    return;
end

x1 = posLine(1);
x2 = posLine(2);
y1 = posLine(3);
y2 = posLine(4);

lineLen = ceil(sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2));
if lineLen < 2
    return;
end

hDrawLine = getappdata(handles.figure_vi, 'hDrawLine');
viewSize = getappdata(handles.figure_vi, 'sliceSize');
if is_handle(hDrawLine)
    set(hDrawLine, ...
        'XData', [x1 x2], ...
        'YData', viewSize(1) - [y1 y2] + 1);
else
    hDrawLine = line( ...
        [x1 x2], viewSize(1) - [y1 y2] + 1, ...
        'Parent', handles.axes_2dViewer, ...
        'Color', 'r');
    
    setappdata(handles.figure_vi, 'hDrawLine', hDrawLine);
end

sliceImage = getappdata(handles.figure_vi, 'sliceImage');
viewSize = getappdata(handles.figure_vi, 'sliceSize');

if ~isfloat(sliceImage)
    sliceImage = double(sliceImage);
end
[cx, cy, c] = improfile(sliceImage, [x1 x2], [y1 y2], lineLen);

hFigLineMeasurement = getappdata(handles.figure_vi, 'hFigLineMeasurement');
if is_handle(hFigLineMeasurement)
    feval(getappdata(hFigLineMeasurement, 'entry_update_data'), ...
        hFigLineMeasurement, {[c cx cy], [1 viewSize(2) 1 viewSize(1)], handles.figure_vi});
end


% -------------------------------------------------------------------
function [pt, nPlane] = pick_3d_point(handles)

assert(4 == get_view_type(handles));

ptOnAxesCurrent = round(get(handles.axes_3dSlicer, 'CurrentPoint'));

hSlices = getappdata(handles.figure_vi, 'hSlices');
xdata = get(hSlices(1), 'XData');
ydata = get(hSlices(2), 'YData');
zdata = get(hSlices(3), 'ZData');

ptPos = zeros(3, 3);

ptPos(1, 1) = xdata(1);
ptPos(1, 2) = interp1nosort(ptOnAxesCurrent(:, 1), ptOnAxesCurrent(:, 2), xdata(1));
ptPos(1, 3) = interp1nosort(ptOnAxesCurrent(:, 1), ptOnAxesCurrent(:, 3), xdata(1));

ptPos(2, 1) = interp1nosort(ptOnAxesCurrent(:, 2), ptOnAxesCurrent(:, 1), ydata(1));
ptPos(2, 2) = ydata(1);
ptPos(2, 3) = interp1nosort(ptOnAxesCurrent(:, 2), ptOnAxesCurrent(:, 3), ydata(1));

ptPos(3, 1) = interp1nosort(ptOnAxesCurrent(:, 3), ptOnAxesCurrent(:, 1), zdata(1));
ptPos(3, 2) = interp1nosort(ptOnAxesCurrent(:, 3), ptOnAxesCurrent(:, 2), zdata(1));
ptPos(3, 3) = zdata(1);

ptPos = round(ptPos);

d2 = zeros(3, 1);

if any(isnan(ptPos(1, :)))
    d2(1) = inf;
else
    d2(1) = sum((ptOnAxesCurrent(1, :) - ptPos(1, :)) .^ 2);
end

if any(isnan(ptPos(2, :)))
    d2(2) = inf;
else
    d2(2) = sum((ptOnAxesCurrent(1, :) - ptPos(2, :)) .^ 2);
end

if any(isnan(ptPos(3, :)))
    d2(3) = inf;
else
    d2(3) = sum((ptOnAxesCurrent(1, :) - ptPos(3, :)) .^ 2);
end

if ~all(isinf(d2))
    [~, nPlane] = min(d2);
    pt = ptPos(nPlane, :);
else
    nPlane = 0;
    pt = [];
end


% --- Outputs from this function are returned to the command line.
function varargout = vi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse motion over figure - except title and menu.
function figure_vi_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure_vi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ptOnAxesCurrent, onImage] = get_pointer_pos(handles);

is2d = 4 ~= get_view_type(handles);

cursorType = get(handles.figure_vi, 'Pointer');

if is2d && onImage && ~strcmpi(cursorType, 'cross')
    set(handles.figure_vi, 'Pointer', 'crosshair');
elseif ~onImage && ~strcmpi(cursorType, 'arrow')
    set(handles.figure_vi, 'Pointer', 'arrow');
end

if ~onImage
    set(handles.text_pixelInfo, 'String', '');
    set(handles.text_voxelInfo, 'String', '');
    return;
end

buttonDown = getappdata(handles.figure_vi, 'buttonDown');
selectionType = get(handles.figure_vi, 'SelectionType');
leftButtonDown = buttonDown && strcmpi(selectionType, 'normal');
rightButtonDown = buttonDown && strcmpi(selectionType, 'alt');

ptOnAxesBtnDown = getappdata(handles.figure_vi, 'ptOnAxesBtnDown');

dataFormat = getappdata(handles.figure_vi, 'dataFormat');

if is2d
    sliceImage = getappdata(handles.figure_vi, 'sliceImage');
    set(handles.text_pixelInfo, 'String', ...
        ['(' num2str(round(ptOnAxesCurrent(1))) ', ' num2str(round(ptOnAxesCurrent(2))) ')  '...
        num2str(sliceImage(ptOnAxesCurrent(2), ptOnAxesCurrent(1)), dataFormat)]);

else
    [pt, nPlane] = pick_3d_point(handles);
    
    if 0 == nPlane
        set(handles.text_voxelInfo, 'String', '');
        
    else
        hSlices = getappdata(handles.figure_vi, 'hSlices');
        sliceImage = get(hSlices(nPlane), 'CData');
        if 1 == nPlane
            v = sliceImage(pt(2), pt(3));
        elseif 2 == nPlane
            v = sliceImage(pt(1), pt(3));
        elseif 3 == nPlane
            v = sliceImage(pt(2), pt(1));
        end
        set(handles.text_voxelInfo, 'String', ...
            ['(' num2str(pt(1)) ', ' num2str(pt(2)) ', ' num2str(pt(3)) ')  '...
            num2str(v, dataFormat)]);
       
    end
    
end

ptOnFigCurrent = get(handles.figure_vi, 'CurrentPoint');

if buttonDown
    ptOnFigBtnDown = getappdata(handles.figure_vi, 'ptOnFigBtnDown');
    ptShift = ptOnFigCurrent - ptOnFigBtnDown;
end

measuringLine = 1 == get(handles.togglebutton_lineMeasure, 'Value');
measuringRegion = 1 == get(handles.togglebutton_regionMeasure, 'Value');

if leftButtonDown && ~measuringLine && ~measuringRegion % Change window.
    % Move mouse down to INCREASE level to DECREASE the brightness.
    % Move mouse left and right to change width.
    referenceCLim = getappdata(handles.figure_vi, 'referenceCLim');
    windowLevel = mean(referenceCLim);
    windowWidth = range(referenceCLim);
    windowLevel = windowLevel - windowWidth * ptShift(2) / 1000;
    windowWidth = max(eps('single'), windowWidth + windowWidth * ptShift(1) / 500);
    windowMax = windowLevel + windowWidth / 2;
    windowMin = windowLevel - windowWidth / 2;
    set_clims(handles, windowMin, windowMax);
    update_window_info(handles);
    
elseif is2d && leftButtonDown && measuringLine && ~measuringRegion % Line measurement.
    setappdata(handles.figure_vi, 'posLine', ...
        [ptOnAxesBtnDown(1), ptOnAxesCurrent(1), ptOnAxesBtnDown(2), ptOnAxesCurrent(2)]);
    update_line_data(handles);
    
elseif is2d && leftButtonDown && ~measuringLine && measuringRegion % Region measurement.
    if strcmpi('Drawn', getappdata(handles.figure_vi, 'regionType'))
        userDrawRegionVertices = getappdata(handles.figure_vi, 'userDrawRegionVertices');
        assert(~isempty(userDrawRegionVertices));
        userDrawRegionVertices = [userDrawRegionVertices; ptOnAxesCurrent];
        setappdata(handles.figure_vi, 'userDrawRegionVertices', userDrawRegionVertices);
        update_user_drawn_region_data(handles);
    else
        setappdata(handles.figure_vi, ...
            'posRegion', [ptOnAxesBtnDown(1), ptOnAxesCurrent(1), ptOnAxesBtnDown(2), ptOnAxesCurrent(2)]);
        update_region_data(handles);
    end
    
elseif ~is2d && 1 == get(handles.togglebutton_rotate, 'Value')
    hLight = getappdata(handles.figure_vi, 'hLight');
    set(hLight, 'Position', get(handles.axes_3dSlicer, 'CameraPosition'));
    
elseif ~is2d && rightButtonDown
    referenceNPlane = getappdata(handles.figure_vi, 'referenceNPlane');
    referenceSliceNo3d = getappdata(handles.figure_vi, 'referenceSliceNo3d');
    
    if referenceNPlane == 1
        update_slice_3d(handles, referenceSliceNo3d(1) + ptShift(2), referenceSliceNo3d(2), referenceSliceNo3d(3), false);
        
    elseif referenceNPlane == 2
        update_slice_3d(handles, referenceSliceNo3d(1), referenceSliceNo3d(2) + ptShift(2), referenceSliceNo3d(3), false);
        
    elseif referenceNPlane == 3
        update_slice_3d(handles, referenceSliceNo3d(1), referenceSliceNo3d(2), referenceSliceNo3d(3) + ptShift(2), false);
        
    else
    end
    
end


% --- Executes on scroll wheel click while the figure is in focus.
function figure_vi_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure_vi (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

if 4 == get_view_type(handles)
    return;
end
    
[~, onImage] = get_pointer_pos(handles);
if ~onImage
    return;
end

numSlices = getappdata(handles.figure_vi, 'numSlices');
sliceNo = getappdata(handles.figure_vi, 'sliceNo');
sliceNo = sliceNo - eventdata.VerticalScrollCount(1);
sliceNo = max(min(sliceNo, numSlices), 1);

setappdata(handles.figure_vi, 'sliceNo', sliceNo);
set(handles.edit_sliceNo, 'String', num2str(sliceNo));

update_slice_2d(handles, sliceNo, false);


function edit_windowMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowMax as text
%        str2double(get(hObject,'String')) returns contents of edit_windowMax as a double

cLim = get(handles.axes_2dViewer, 'CLim');
windowMin = cLim(1);
windowMax = str2double(get(hObject,'String'));

if windowMax < windowMin
    windowMax = windowMin;
end
set_clims(handles, windowMin, windowMax);

update_window_info(handles);


% --- Executes during object creation, after setting all properties.
function edit_windowMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_windowMin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowMin as text
%        str2double(get(hObject,'String')) returns contents of edit_windowMin as a double

cLim = get(handles.axes_2dViewer, 'CLim');
windowMax = cLim(2);
windowMin = str2double(get(hObject,'String'));

if windowMin > windowMax
    windowMin = windowMax;
end
set_clims(handles, windowMin, windowMax);

update_window_info(handles);


% --- Executes during object creation, after setting all properties.
function edit_windowMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_windowWidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowWidth as text
%        str2double(get(hObject,'String')) returns contents of edit_windowWidth as a double

windowWidth = str2double(get(hObject,'String'));
if windowWidth < 1
    windowWidth = 1;
end
cLim = get(handles.axes_2dViewer, 'CLim');
windowLevel = mean(cLim);
windowMax = windowLevel + windowWidth / 2;
windowMin = windowLevel - windowWidth / 2;

set_clims(handles, windowMin, windowMax);

update_window_info(handles);


% --- Executes during object creation, after setting all properties.
function edit_windowWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_windowLevel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowLevel as text
%        str2double(get(hObject,'String')) returns contents of edit_windowLevel as a double

windowLevel = str2double(get(hObject,'String'));
cLim = get(handles.axes_2dViewer, 'CLim');
windowWidth = range(cLim);
windowMax = windowLevel + windowWidth / 2;
windowMin = windowLevel - windowWidth / 2;

set_clims(handles, windowMin, windowMax);

update_window_info(handles);


% --- Executes during object creation, after setting all properties.
function edit_windowLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_sliceNo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sliceNo as text
%        str2double(get(hObject,'String')) returns contents of edit_sliceNo as a double

update_slice_2d(handles, round(str2double(get(hObject,'String'))), false);


% --- Executes during object creation, after setting all properties.
function edit_sliceNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_sliceNo_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

update_slice_2d(handles, get(hObject, 'Value'), false);


% --- Executes during object creation, after setting all properties.
function slider_sliceNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when selected object is changed in uipanel_viewType.
function uipanel_viewType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_viewType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

update_view_type(handles);
auto_layout(handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure_vi_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure_vi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setappdata(handles.figure_vi, 'buttonDown', true);
setappdata(handles.figure_vi, 'ptOnFigBtnDown', get(handles.figure_vi, 'CurrentPoint'));

[ptOnAxesCurrent, onImage] = get_pointer_pos(handles);
if ~onImage
    return;
end

setappdata(handles.figure_vi, 'ptOnAxesBtnDown', ptOnAxesCurrent);

selectionType = get(handles.figure_vi, 'SelectionType');

if strcmpi(selectionType, 'normal')
    hDrawRegion = getappdata(handles.figure_vi, 'hDrawRegion');
    if is_handle(hDrawRegion)
        delete(hDrawRegion);
        setappdata(handles.figure_vi, 'hDrawRegion', []);
    end
    
    hDrawLine = getappdata(handles.figure_vi, 'hDrawLine');
    if is_handle(hDrawLine)
        delete(hDrawLine);
        setappdata(handles.figure_vi, 'hDrawLine', []);
    end
    
    hDrawRegion = getappdata(handles.figure_vi, 'hDrawRegion');
    if ishandle(hDrawRegion)
        delete(hDrawRegion);
        setappdata(handles.figure_vi, 'hDrawRegion', []);
    end
    
    hDrawRegionAdd = getappdata(handles.figure_vi, 'hDrawRegionAdd');
    if ishandle(hDrawRegionAdd)
        delete(hDrawRegionAdd);
        setappdata(handles.figure_vi, 'hDrawRegionAdd', []);
    end
    
    cLim = get(handles.axes_2dViewer, 'CLim');
    setappdata(handles.figure_vi, 'referenceCLim', cLim);
    
    if strcmpi('Drawn', getappdata(handles.figure_vi, 'regionType')) ...
            && 1 == get(handles.togglebutton_regionMeasure, 'Value')
        setappdata(handles.figure_vi, 'userDrawRegionVertices', ptOnAxesCurrent);
    end
    
elseif strcmpi(selectionType, 'alt')
    [~, nPlane] = pick_3d_point(handles);
    
    setappdata(handles.figure_vi, 'referenceNPlane', nPlane);
    
    sliceNo3d = getappdata(handles.figure_vi, 'sliceNo3d');
    setappdata(handles.figure_vi, 'referenceSliceNo3d', sliceNo3d);
    
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure_vi_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure_vi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setappdata(handles.figure_vi, 'buttonDown', false);

hDrawRegionAdd = getappdata(handles.figure_vi, 'hDrawRegionAdd');
if is_handle(hDrawRegionAdd) && strcmpi('on', get(hDrawRegionAdd, 'Visible'))
    set(hDrawRegionAdd, 'LineStyle', '-', 'Color', 'r');
end


function edit_zAspectRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zAspectRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_zAspectRatio as a double

aspectRatio = getappdata(handles.figure_vi, 'aspectRatio');
aspectRatio(3) = str2double(get(hObject, 'String'));
setappdata(handles.figure_vi, 'aspectRatio', aspectRatio);
update_aspect_ratio(handles);
auto_layout(handles);


% --- Executes during object creation, after setting all properties.
function edit_zAspectRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yAspectRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yAspectRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_yAspectRatio as a double

aspectRatio = getappdata(handles.figure_vi, 'aspectRatio');
aspectRatio(2) = str2double(get(hObject, 'String'));
setappdata(handles.figure_vi, 'aspectRatio', aspectRatio);
update_aspect_ratio(handles);
auto_layout(handles);


% --- Executes during object creation, after setting all properties.
function edit_yAspectRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xAspectRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xAspectRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_xAspectRatio as a double

aspectRatio = getappdata(handles.figure_vi, 'aspectRatio');
aspectRatio(1) = str2double(get(hObject, 'String'));
setappdata(handles.figure_vi, 'aspectRatio', aspectRatio);
update_aspect_ratio(handles);
auto_layout(handles);


% --- Executes during object creation, after setting all properties.
function edit_xAspectRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_resetAspectRatio.
function pushbutton_resetAspectRatio_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_resetAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

aspectRatio = [1 1 1];
setappdata(handles.figure_vi, 'aspectRatio', aspectRatio);
update_aspect_ratio(handles);
auto_layout(handles);


% --- Executes on button press in pushbutton_moreWindows.
function pushbutton_moreWindows_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_moreWindows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiMenus = getappdata(handles.figure_vi, 'uiMenus');
pos1 = get(handles.pushbutton_moreWindows, 'Position');
pos2 = get(handles.uipanel_windowSetting, 'Position');
pos3 = get(handles.uipanel_leftControls, 'Position');
set(uiMenus.menu_moreWindowSettings, 'Position', pos1(1:2)+pos2(1:2)+pos3(1:2), 'Visible', 'on');


% --- Executes on key press with focus on figure_vi and none of its controls.
function figure_vi_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure_vi (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(eventdata.Character, '+')
    aspectRatio = getappdata(handles.figure_vi, 'aspectRatio');
    if isempty(eventdata.Modifier)
        aspectRatio = aspectRatio * 1.1;
    elseif strcmpi(eventdata.Modifier{1}, 'control')
        aspectRatio = aspectRatio * 1.5;
    end
    setappdata(handles.figure_vi, 'aspectRatio', aspectRatio);
    update_aspect_ratio(handles);
    auto_layout(handles);
elseif strcmpi(eventdata.Character, '-')
    aspectRatio = getappdata(handles.figure_vi, 'aspectRatio');
    if isempty(eventdata.Modifier)
        aspectRatio = aspectRatio / 1.1;
    elseif strcmpi(eventdata.Modifier{1}, 'control')
        aspectRatio = aspectRatio / 1.5;
    end
    if all(aspectRatio > 0.2)
        setappdata(handles.figure_vi, 'aspectRatio', aspectRatio);
        update_aspect_ratio(handles);
        auto_layout(handles);
    end
elseif strcmpi(eventdata.Character, '*')
    aspectRatio = [1 1 1];
    setappdata(handles.figure_vi, 'aspectRatio', aspectRatio);
    update_aspect_ratio(handles);
    auto_layout(handles);
else
end


% --- Executes on selection change in listbox_colorMaps.
function listbox_colorMaps_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_colorMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_colorMaps contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_colorMaps

set_colormap(handles);


% --- Executes during object creation, after setting all properties.
function listbox_colorMaps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_colorMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_xSliceNo_Callback(hObject, eventdata, handles)
% hObject    handle to slider_xSliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliceNo3d = getappdata(handles.figure_vi, 'sliceNo3d');
update_slice_3d(handles, get(hObject, 'Value'), sliceNo3d(2), sliceNo3d(3), false);


% --- Executes during object creation, after setting all properties.
function slider_xSliceNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_xSliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_ySliceNo_Callback(hObject, eventdata, handles)
% hObject    handle to slider_ySliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliceNo3d = getappdata(handles.figure_vi, 'sliceNo3d');
update_slice_3d(handles, sliceNo3d(1), get(hObject, 'Value'), sliceNo3d(3), false);


% --- Executes during object creation, after setting all properties.
function slider_ySliceNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_ySliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_zSliceNo_Callback(hObject, eventdata, handles)
% hObject    handle to slider_zSliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliceNo3d = getappdata(handles.figure_vi, 'sliceNo3d');
update_slice_3d(handles, sliceNo3d(1), sliceNo3d(2), get(hObject, 'Value'), false);


% --- Executes during object creation, after setting all properties.
function slider_zSliceNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_zSliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function mi_windowPelvis_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowPelvis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -160, 240);
update_window_info(handles);


% --------------------------------------------------------------------
function mi_windowBone_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowBone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -750, 1750);
update_window_info(handles);


% --------------------------------------------------------------------
function mi_windowBreast_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowBreast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -250, 150);
update_window_info(handles);


% --------------------------------------------------------------------
function mi_windowAbdomen_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowAbdomen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -125, 225);
update_window_info(handles);


% --------------------------------------------------------------------
function mi_windowLung_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowLung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -1000, 250);
update_window_info(handles);


% --------------------------------------------------------------------
function mi_windowLiver_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowLiver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -25, 125);
update_window_info(handles);


% --------------------------------------------------------------------
function mi_windowCerebellum_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowCerebellum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -20, 100);
update_window_info(handles);


% --------------------------------------------------------------------
function mi_windowReset_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img = getappdata(handles.figure_vi, 'imgData');
maxPixelVal = max(img(:));
minPixelVal = min(img(:));
set_clims(handles, minPixelVal, maxPixelVal);
update_window_info(handles);


% --------------------------------------------------------------------
function menu_moreWindowSettings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_moreWindowSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mi_windowCustom1_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowCustom1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -1000, 1000);
update_window_info(handles);


% --------------------------------------------------------------------
function mi_windowCustom2_Callback(hObject, eventdata, handles)
% hObject    handle to mi_windowCustom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_clims(handles, -1500, 1500);
update_window_info(handles);


% --- Executes on button press in togglebutton_rotate.
function togglebutton_rotate_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_rotate

check_rotate(handles);


% --- Executes on button press in togglebutton_light.
function togglebutton_light_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_light (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_light

set_light(handles);


% --- Executes on button press in togglebutton_sliceStats.
function togglebutton_sliceStats_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_sliceStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_sliceStats

if 1 == get(hObject, 'Value')
    set([ ...
        handles.togglebutton_regionMeasure, ...
        handles.togglebutton_lineMeasure, ...
        handles.togglebutton_isoline], 'Value', 0);
    set_region_measure(handles, 0);
    set_line_measure(handles, 0);
    set_isoline(handles, 0);
    
    set_image_stats(handles, 1);
else
    set_image_stats(handles, 0);
end


% --- Executes on button press in togglebutton_regionMeasure.
function togglebutton_regionMeasure_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_regionMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_regionMeasure

if 1 == get(hObject, 'Value')
    set([ ...
        handles.togglebutton_sliceStats, ...
        handles.togglebutton_lineMeasure, ...
        handles.togglebutton_isoline], 'Value', 0);
    set_image_stats(handles, 0);
    set_line_measure(handles, 0);
    set_isoline(handles, 0);

    set_region_measure(handles, 1);
else
    set_region_measure(handles, 0);
end


% --- Executes on button press in togglebutton_lineMeasure.
function togglebutton_lineMeasure_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_lineMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_lineMeasure

if 1 == get(hObject, 'Value')
    set([ ...
        handles.togglebutton_sliceStats, ...
        handles.togglebutton_regionMeasure, ...
        handles.togglebutton_isoline], 'Value', 0);
    set_image_stats(handles, 0);
    set_region_measure(handles, 0);
    set_isoline(handles, 0);

    set_line_measure(handles, 1);
else
    set_line_measure(handles, 0);
end



% --- Executes on button press in togglebutton_isoline.
function togglebutton_isoline_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_isoline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_isoline

if 1 == get(hObject, 'Value')
    set([ ...
        handles.togglebutton_sliceStats, ...
        handles.togglebutton_regionMeasure, ...
        handles.togglebutton_lineMeasure], 'Value', 0);
    set_image_stats(handles, 0);
    set_region_measure(handles, 0);
    set_line_measure(handles, 0);
    
    set_isoline(handles, 1);
else
    set_isoline(handles, 0);
end


% --- Executes when figure_vi is resized.
function figure_vi_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure_vi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

auto_layout(handles);


% --- Executes when user attempts to close figure_vi.
function figure_vi_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_vi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFigSliceStats = getappdata(handles.figure_vi, 'hFigSliceStats');
if is_handle(hFigSliceStats)
    delete(hFigSliceStats);
    setappdata(handles.figure_vi, 'hFigSliceStats', []);
end

hFigRegionMeasurement = getappdata(handles.figure_vi, 'hFigRegionMeasurement');
if is_handle(hFigRegionMeasurement)
    delete(hFigRegionMeasurement);
    setappdata(handles.figure_vi, 'hFigRegionMeasurement', []);
end

hFigLineMeasurement = getappdata(handles.figure_vi, 'hFigLineMeasurement');
if is_handle(hFigLineMeasurement)
    delete(hFigLineMeasurement);
    setappdata(handles.figure_vi, 'hFigLineMeasurement', []);
end

hFigIsoline = getappdata(handles.figure_vi, 'hFigIsoline');
if is_handle(hFigIsoline)
    delete(hFigIsoline);
    setappdata(handles.figure_vi, 'hFigIsoline', []);
end

delete(hObject);


% --------------------------------------------------------------------
function mi_regionRect_Callback(hObject, eventdata, handles)
% hObject    handle to mi_regionRect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

select_region_type(handles, 'Rectangle');


% --------------------------------------------------------------------
function mi_regionDisc_Callback(hObject, eventdata, handles)
% hObject    handle to mi_regionDisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

select_region_type(handles, 'Disc');


% -------------------------------------------------------------------
function select_region_type(handles, regionType)

uiMenus = getappdata(handles.figure_vi, 'uiMenus');
menuItems = get(uiMenus.menu_selectRegionType, 'Children');
for i = 1 : length(menuItems)
    mi = menuItems(i);
    miLabel = get(mi, 'Label');
    if strcmpi(regionType, miLabel) || strcmpi(regionType, 'Drawn') && strcmpi(miLabel, 'User Drawn')
        set(mi, 'Checked', 'on');
        setappdata(handles.figure_vi, 'regionType', regionType);
    else
        set(mi, 'Checked', 'off');
    end
end


% --------------------------------------------------------------------
function menu_selectRegionType_Callback(hObject, eventdata, handles)
% hObject    handle to menu_selectRegionType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over togglebutton_regionMeasure.
function togglebutton_regionMeasure_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to togglebutton_regionMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiMenus = getappdata(handles.figure_vi, 'uiMenus');
pos1 = get(handles.togglebutton_regionMeasure, 'Position');
pos2 = get(handles.uipanel_2dControls, 'Position');
pos3 = get(handles.uipanel_2dViewer, 'Position');
set(uiMenus.menu_selectRegionType, 'Position', pos1(1:2)+pos2(1:2)+pos3(1:2), 'Visible', 'on');


% --------------------------------------------------------------------
function mi_regionDrawn_Callback(hObject, eventdata, handles)
% hObject    handle to mi_regionDrawn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

select_region_type(handles, 'Drawn');
