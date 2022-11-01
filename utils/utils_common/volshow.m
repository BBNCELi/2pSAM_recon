% VOLSHOW is a flexible function for displaying multiple volume images tightly
%         on the same figure. Mouse scroll wheel scrolls the z dimension.
%         Padding between images, grid dimensions, contrast scale, and
%         colourmaps can be specified. Attributes apply to all images. Best
%         results with same sized images. Grayscale or colour images.
% 
% Input arguments: (any order)
%    images(s) - any number of 3D grayscale or colour images. Rendered in the
%                order they are presented, top to bottom, left to right. 
%                * The x-dimension of any image should not have a size of 3,
%                  else it will be confused for a colourmap.
%                * multiple volumes should have the same 3rd dimension.
%                * if using colour volumes, the colour channels should be the
%                  last dimension.
% 
%    padval    - decimal value on the interval (0, 0.5) dictating the relative
%                padded spacing between images.
%                Default: 0.005
% 
%    gridstr   - string like "5x2", specifying the number of images to tile
%                horizontally (5) and vertically (2)
%                Default: as square as possible, wider bias
% 
%    minmax    - minmax specification for contrast scaling, as in imshow(I,[]).
%                array of size: 1 by 2, or a empty array: []
%                Default: []
% 
%    colourmap - colourmap used for displaying images:
%                array of size: M by 3 or a colourmap function
%                Default: curent default figure colormap
% 
%              * if 2+ non-image arguments are given, only the last one is used.
% 
% Examples:
% 
%    timshow(I1, I2, I3, I4, hot, 0, [0,1], '4x1');
%                Show volumes I1, I2, I3, I4 using the hot colourmap, with no
%                space between, contrast from 0 to 1, and in a horizontal line.
% 
%    timshow(DB(:).I);
%                Show all volume fields .I in the struct array DB using the
%                default figure colourmap, automatic contrast scaling per image,
%                with 0.5% of total figure size padded between, and arranged as
%                close to square as possible.
% 
% Jesse Knight 2015

function varargout = volshow(varargin)
% -- do not edit: MATLAB GUIDE 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @volshow_OpeningFcn, ...
    'gui_OutputFcn',  @volshow_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end % -- end do not edit

% --- Executes just before volshow is made visible.
function volshow_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output  = hObject;
handles.badcall = 0;
set(0,'defaultTextFontName','Courier New');

try % if any errors: abort and give invalid input warning
    % default values
    handles.img       = [];
    handles.minmax    = [];
    handles.colourmap = get(0,'defaultfigurecolormap');
    handles.pad       = 0.005;
    
    % handle input arguments based on their dimensions
    for v = 1:numel(varargin)
        sizev = size(varargin{v});
        % image volume (possibly including colour in 3rd dimension)
        if numel(sizev) >= 3
            handles.img(end+1).data  = varargin{v};
            handles.img(end).size    = size(handles.img(end).data);
            handles.img(end).frame   = round(handles.img(end).size(end)/2);
            handles.img(end).textpos = [round(handles.img(end).size(1)/20 + 1),...
                                        round(handles.img(end).size(2)/20 + 1)];
        % padval
        elseif all(sizev == [1,1])
            handles.pad = varargin{v};
        % minmax (numerical)
        elseif all(sizev == [1,2])
            handles.minmax = varargin{v};
        % minmax ([])
        elseif sizev(1) == 0
            handles.minmax = [];
        % colourmap
        elseif sizev(2) == 3
            handles.colourmap = varargin{v};
        % gridstr
        elseif ischar(varargin{v}) && numel(sscanf(varargin{v},'%dx%d')) == 2
            xy = sscanf(varargin{v},'%dx%d');
            handles.nSubx = xy(1);
            handles.nSuby = xy(2);
        % argument not recognized: ignoring
        else
            warning(['Ignoring argument number ',num2str(v),'.']);
        end
    end
    % optimize display grid square-ish if not user specified
    handles.N = numel(handles.img);
    if ~all(isfield(handles,{'nSubx','nSuby'}))
        nSubx = ceil(sqrt(handles.N));
        nSuby = ceil(handles.N/nSubx);
    end
    % subplot handles initialization & spacing
    for a = 1:handles.N
        y = ceil(a / nSubx);
        x = mod(a, nSubx);
        x(~x) = nSubx;
        
        handles.ax(a) = subplot(nSuby,nSubx,a);
        set(handles.ax(a),'position',[(x - 1) / nSubx + 0.5*handles.pad,  ...
                                       1 - (y / nSuby - 0.5*handles.pad), ...
                                            1 / nSubx - handles.pad,      ...
                                            1 / nSuby - handles.pad]);
    end
    % optimize figure display size for the current monitor and first image size
    % centres the figure in onscreen too.
    screensize = get(0,'screensize');
    imgSize = min(500, (0.4*screensize(3)) / nSubx);
    set(gcf,'position',...
        [(screensize(3) - (imgSize*nSubx))/2,...
         (screensize(4) - (imgSize*nSuby))/2,...
         (imgSize*nSubx),...
         (imgSize*nSuby)]);
    % render the middle frame of each volume to start
    imupdate(handles);
% input argument parsing failed: exit (could be more graceful)
catch
    warning('Error using input arguments.');
    handles.badcall = 1;
end

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = volshow_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
if handles.badcall
    figure1_CloseRequestFcn(hObject, eventdata, handles);
end

% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% for all volumes
for i = 1:numel(handles.img)
    % adjust the frame index by the scroll count
    handles.img(i).frame = handles.img(i).frame + eventdata.VerticalScrollCount;
    % wrap around if z is less than 1 or larger than img.size
    if handles.img(i).frame > handles.img(i).size(3)
        handles.img(i).frame = 1;
    elseif handles.img(i).frame < 1
        handles.img(i).frame = handles.img(i).size(3);
    end
end
% update the frames onscreen
guidata(hObject, handles);
imupdate(handles);

% --- Called by other functions on WindowScrollWheelFcn movement.
function imupdate(handles)
% for all volumes
for i = 1:handles.N
    % show the current frame
    imshow(squeeze(handles.img(i).data(:,:,handles.img(i).frame,:)),handles.minmax,...
        'parent',handles.ax(i));
    % print the current frame number in the top left corner
    text(handles.img(i).textpos(2),handles.img(i).textpos(1),...
        num2str(handles.img(i).frame),'color','r','parent',handles.ax(i));
end
colormap(handles.colourmap);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);
