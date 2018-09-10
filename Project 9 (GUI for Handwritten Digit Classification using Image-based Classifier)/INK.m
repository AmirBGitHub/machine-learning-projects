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
global TrainingData
pnt = [];
Npnt = 1;

TrainingDataRead = xlsread('Digit Templates');
    
for c = 1:10
    p = 1;
    if sum(isnan(TrainingDataRead(:,2*(c-1)+1))) == 0
        TrainingData{1,c} = TrainingDataRead(:,2*(c-1)+1:2*(c-1)+2);
    elseif sum(isnan(TrainingDataRead(:,2*(c-1)+1))) ~= 0
        while  isnan(TrainingDataRead(p,2*(c-1)+1)) == 0
            p = p + 1;
        end
        TrainingData{1,c} = TrainingDataRead(1:p-1,2*(c-1)+1:2*(c-1)+2);
    end    
end






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

global pnt
global TrainingData

INK_Data = [pnt(:,1) pnt(:,2)];

[index, score] = Image_Based_Classifier(INK_Data,TrainingData);

set(handles.Digit_Class, 'String', mod(index,10));



% --- Executes during object creation, after setting all properties.
function Digit_Class_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Digit_Class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
