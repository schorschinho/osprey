function  imshow3D( Img, disprange, initS )
%IMSHOW3D displays 3D grayscale or RGB images in a slice by slice fashion
%with mouse-based slice browsing and window and level adjustment control,
%and auto slice browsing control.
%
% Usage:
% imshow3D ( Image )
% imshow3D ( Image , [] )
% imshow3D ( Image , [LOW HIGH] )
% imshow3D ( Image , [] , initsn )
%   
%    Image:      3D image MxNxKxC (K slices of MxN images) C is either 1
%                (for grayscale images) or 3 (for RGB images)  
%    [LOW HIGH]: display range that controls the display intensity range of
%                a grayscale image (default: the broadest available range)
%    initsn:     The slice number to be displayed initially (default:
%                mid-slice number) 
%
% Use the scroll bar or mouse scroll wheel to switch between slices. To
% adjust window and level values keep the mouse right button pressed, and
% drag the mouse up and down (for level adjustment) or right and left (for
% window adjustment). Window and level adjustment control works only for
% grayscale images.
% "Play" button displays all the slices as a sequence of frames. The time
% interval value can also be adjusted (default time interval is 100 ms).
% 
% "Auto W/L" button adjust the window and level automatically for grayscale
% images.
%
% While "Fine Tune" checkbox is checked the window/level adjustment gets 16
% times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse-based window and level adjustment is set
% based on the user-defined display intensity range; the wider the range,
% the more sensitivity to mouse drag.
% 
% Note: IMSHOW3DFULL is a newer version of IMSHOW3D (also available on
% MathWorks) that displays 3D grayscale or RGB images from three
% perpendicular views (i.e., axial, sagittal, and coronal).
% 
%   Example
%   --------
%       % To display an image (MRI example)
%       load mri 
%       Image = squeeze(D); 
%       figure, 
%       imshow3D(Image) 
%
%       % To display the image, and adjust the display range
%       figure,
%       imshow3D(Image,[20 100]);
%
%       % To define the initial slice number
%       figure,
%       imshow3D(Image,[],5);
%
%   See also IMSHOW.

%
% - Maysam Shahedi (mshahedi@gmail.com)
% - Released: 1.0.0   Date: 2013/04/15
% - Revision: 1.1.0   Date: 2013/04/19
% - Revision: 1.5.0   Date: 2016/09/22
% - Revision: 1.6.0   Date: 2018/06/07
% - Revision: 1.6.1   Date: 2018/10/29
% 

sno = size(Img,3);  % number of slices
S = round(sno/2);

PlayFlag = false;   % Play flag, playing when it is 'True'
Tinterv = 100;

global InitialCoord;

MinV = 0;
MaxV = max(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end    

SFntSz = 9;
txtFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;

if (nargin < 3)
    S = round(sno/2);
else
    S = initS;
    if S > sno
        S = sno;
        warning('Initial slice number out of range');
    elseif S < 1
        S = 1;
        warning('Initial slice number out of range');
    end
end

if (nargin < 2)
    [Rmin Rmax] = WL2R(Win, LevV);
elseif numel(disprange) == 0
    [Rmin Rmax] = WL2R(Win, LevV);
else
    LevV = (double(disprange(2)) + double(disprange(1))) / 2;
    Win = double(disprange(2)) - double(disprange(1));
    WLAdjCoe = (Win + 1)/1024;
    [Rmin Rmax] = WL2R(Win, LevV);
end

clf
axes('position',[0,0.2,1,0.8]), imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])

FigPos = get(gcf,'Position');
S_Pos = [30 45 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [30 65 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [20 18 60 20];
Wval_Pos = [75 20 50 20];
Ltxt_Pos = [130 18 45 20];
Lval_Pos = [170 20 50 20];
Btn_Pos = [240 20 70 20];
ChBx_Pos = [320 20 80 20];
Play_Pos = [uint16(FigPos(3)-100)+40 45 30 20];
Time_Pos = [uint16(FigPos(3)-100)+35 20 40 20];
Ttxt_Pos = [uint16(FigPos(3)-100)-50 18 90 20];

% W/L Button styles:
WL_BG = ones(Btn_Pos(4),Btn_Pos(3),3)*0.85;
WL_BG(1,:,:) = 1; WL_BG(:,1,:) = 1; WL_BG(:,end-1,:) = 0.6; WL_BG(:,end,:) = 0.4; WL_BG(end,:,:) = 0.4;

% Play Button styles:
Play_BG = ones(Play_Pos(4),Play_Pos(3),3)*0.85;
Play_BG(1,:,:) = 1; Play_BG(:,1,:) = 1; Play_BG(:,end-1,:) = 0.6; Play_BG(:,end,:) = 0.4; Play_BG(end,:,:) = 0.4;
Play_Symb = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1;...
             0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1;...
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1;...
             0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;...
             0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
Play_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-7:floor(Play_Pos(3)/2)+6,:) = ...
    repmat(Play_Symb,1,1,3) .* Play_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-7:floor(Play_Pos(3)/2)+6,:);
Pause_BG = ones(Play_Pos(4),Play_Pos(3),3)*0.85;
Pause_BG(1,:,:) = 1; Pause_BG(:,1,:) = 1; Pause_BG(:,end-1,:) = 0.6; Pause_BG(:,end,:) = 0.4; Pause_BG(end,:,:) = 0.4;
Pause_Symb = repmat([0, 0, 0, 1, 1, 1, 1, 0, 0, 0],13,1);
Pause_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-5:floor(Play_Pos(3)/2)+4,:) = ...
    repmat(Pause_Symb,1,1,3) .* Pause_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-5:floor(Play_Pos(3)/2)+4,:);


if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'FontSize', SFntSz);
    playhand = uicontrol('Style', 'pushbutton','Position', Play_Pos, 'Callback' , @Play);
    set(playhand, 'cdata', Play_BG)
    ttxthand = uicontrol('Style', 'text','Position', Ttxt_Pos,'String','Interval (ms): ',  'FontSize', txtFntSz);
    timehand = uicontrol('Style', 'edit','Position', Time_Pos,'String',sprintf('%d',Tinterv), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @TimeChanged);
else
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'FontSize', SFntSz);
end    
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ',  'FontSize', txtFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ',  'FontSize', txtFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
set(Btnhand, 'cdata', WL_BG)
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine-tune', 'FontSize', txtFntSz);

set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [30 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [30 65 uint16(FigPos(3)-100)+1 15];
        Play_Pos = [uint16(FigPos(3)-100)+40 45 30 20];
        Time_Pos = [uint16(FigPos(3)-100)+35 20 40 20];
        Ttxt_Pos = [uint16(FigPos(3)-100)-50 18 90 20];
        if sno > 1
            set(shand,'Position', S_Pos);
            set(playhand, 'Position', Play_Pos)
            set(ttxthand, 'Position', Ttxt_Pos)
            set(timehand, 'Position', Time_Pos)
        end
        set(stxthand,'Position', Stxt_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', '')
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end

% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;

        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)

        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end

% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end

% -=< Play button callback function >=-
    function Play (hObj,event)
        PlayFlag = ~PlayFlag;
        if PlayFlag
            set(playhand, 'cdata', Pause_BG)
        else
            set(playhand, 'cdata', Play_BG)
        end            
        while PlayFlag
            S = S + 1;
            if (S > sno)
                S = 1;
            end
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
            set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
            pause(Tinterv/1000)
        end
    end

% -=< Time interval adjustment callback function>=-
    function TimeChanged(varargin)
        Tinterv = str2double(get(timehand, 'string'));
    end
    
end
% -=< Maysam Shahedi (mshahedi@gmail.com), October 29, 2018>=-