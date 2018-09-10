function varargout = INK(varargin)
% INK MATLAB code for INK.fig
%      INK, by itself, creates a new INK or raises the existing
%      singleton*.
%
%      H = INK returns the handle to a new INK or the handle to
%      the existing singleton*.
%
%      INK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INK.M with the given input arguments.
%
%      INK('Property','Value',...) creates a new INK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before INK_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to INK_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help INK

% Last Modified by GUIDE v2.5 02-Mar-2016 17:44:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @INK_OpeningFcn, ...
                   'gui_OutputFcn',  @INK_OutputFcn, ...
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


% --- Executes just before INK is made visible.
function INK_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to INK (see VARARGIN)

% Choose default command line output for INK
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes INK wait for user response (see UIRESUME)
% uiwait(handles.figure1);


global drawing;
drawing =0;
set(gcf,'WindowButtonDownFcn',@mouseDown)
set(gcf,'WindowButtonMotionFcn',@mouseMove)
set(gcf,'WindowButtonUpFcn',@mouseUp)

global pnt
global Npnt
pnt = [];
Npnt = 1;
tic

% --- Outputs from this function are returned to the command line.
function varargout = INK_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla
global pnt
global Npnt
pnt = [];
Npnt = 1;
set(handles.Digit_Class, 'String', '');

% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pnt
global Npnt

if Npnt<1000 
pnt(Npnt+1:end,:) =[];
end

dlmwrite('InkData.txt', pnt, 'delimiter', '\t')
type('InkData.txt')
set(handles.Digit_Class, 'String', '');

function mouseDown(hObject, eventdata, handles) 
global drawing
drawing = 1;

function mouseUp(hObject, eventdata, handles) 
global drawing
drawing = 0;

function mouseMove(hObject, eventdata, handles) 
global drawing
global Npnt
global pnt

if drawing
    C = get(gca,'CurrentPoint');
    if C(1,1)<1 && C(1,1)>0 && C(1,2)<1 && C(1,2)>0
        pnt(Npnt,1) = C(1,1);
        pnt(Npnt,2) = C(1,2);
        pnt(Npnt,3) = toc;
        Npnt = Npnt+1;
        plot(C(1,1),C(1,2),'k','marker','o','MarkerFaceColor','r')
        hold on
        xlim([0 1]); ylim([0 1]);
        set(gca,'XTick',[],'YTick',[])
        box on
    end
end


% --- Executes on button press in Classify_Digit.
function Classify_Digit_Callback(hObject, eventdata, handles)
% hObject    handle to Classify_Digit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Npnt
global pnt
global class

Npnt = length(pnt(:,3))-1;
x = pnt(:,1);
y = pnt(:,2);
t = pnt(:,3);

f(1) = (x(3)-x(1))/sqrt((x(3)-x(1))^2+(y(3)-y(1))^2);
f(2) = (y(3)-y(1))/sqrt((x(3)-x(1))^2+(y(3)-y(1))^2);
f(3) = sqrt((max(x)-min(x))^2+(max(y)-min(y))^2);
f(4) = atan2((max(y)-min(y))/(max(x)-min(x)),1);
f(5) = sqrt((x(Npnt)-x(1))^2+(y(Npnt)-y(1))^2);
f(6) = (x(Npnt)-x(1))/f(5);
f(7) = (y(Npnt)-y(1))/f(5);
f(8) = 0;
f(9) = 0;
f(10) = 0;
f(11) = 0;
for p = 1:Npnt-1
    dx(p) = x(p+1)-x(p);
    dy(p) = y(p+1)-y(p);
    dt(p) = t(p+1)-t(p);
    f(8) = f(8) + sqrt(dx(p)^2+dy(p)^2); 
end
for p = 2:Npnt-1
    f(9) = f(9) + atan2((dx(p)*dy(p-1)-dx(p-1)*dy(p))/(dx(p)*dx(p-1)+dy(p)*dy(p-1)),1);
    f(10) = f(10) + abs(atan2((dx(p)*dy(p-1)-dx(p-1)*dy(p))/(dx(p)*dx(p-1)+dy(p)*dy(p-1)),1));
    f(11) = f(11) + (atan2((dx(p)*dy(p-1)-dx(p-1)*dy(p))/(dx(p)*dx(p-1)+dy(p)*dy(p-1)),1))^2;
end    

f(12) = max((dx(1:Npnt-1).^2+dy(1:Npnt-1).^2)./dt(1:Npnt-1).^2);
f(13) = t(Npnt)-t(1);

load('w.mat');
load('w0.mat');

[v c] = max(w0 + w*f');
class = mod(c,10);

set(handles.Digit_Class, 'String', class);



% --- Executes during object creation, after setting all properties.
function Digit_Class_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Digit_Class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
