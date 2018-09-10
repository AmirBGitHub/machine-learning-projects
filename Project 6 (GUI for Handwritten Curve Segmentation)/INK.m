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

% Last Modified by GUIDE v2.5 20-Feb-2016 17:41:56

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
global pnt
global Npnt
drawing =0;

set(gcf,'WindowButtonDownFcn',@mouseDown)
set(gcf,'WindowButtonMotionFcn',@mouseMove)
set(gcf,'WindowButtonUpFcn',@mouseUp)

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


function mouseDown(hObject, eventdata, handles) 
global drawing
drawing = 1;

function mouseUp(hObject, eventdata, handles) 
global drawing
global Npnt
global pnt
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
        plot(C(1,1),C(1,2),'k','marker','.','MarkerFaceColor','r')
        hold on
        xlim([0 1]); ylim([0 1]);
        set(gca,'XTick',[],'YTick',[])
        box on
    end
end


% --- Executes on button press in Segment.
function Segment_Callback(hObject, eventdata, handles)
% hObject    handle to Segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Npnt
global pnt

Npnt = length(pnt(:,3));

dist = zeros(Npnt,1);
spd = zeros(Npnt,1);
spdavg = zeros(Npnt,1);
Smooth_Speed = zeros(Npnt,1);                                                           % Finding distance passed and raw speed
for i = 2:Npnt-1
    dist(i+1,1) = dist(i,1) + norm([pnt(i,1)-pnt(i-1,1), pnt(i,2)-pnt(i-1,2)],2);
    spd(i,1) = abs((dist(i+1) - dist(i-1))/(pnt(i+1,3) - pnt(i-1,3))); 
end

for i = 3:Npnt-2
    spdavg(i:Npnt-3,1) = mean(spd(i-2:i+2,1));
end
                                                                                        % Smoothing the speed
Smooth_Speed = [mean([spd(1,1) spd(Npnt,1)]); mean([spd(2,1) spd(Npnt-1,1)]); mean([spd(3,1) spd(Npnt-2,1)]); spdavg(4:Npnt-3,1); mean([spd(3,1) spd(Npnt-2,1)]); mean([spd(2,1) spd(Npnt-1,1)]); mean([spd(1,1) spd(Npnt,1)])];  
Average_Speed = dist(Npnt)/(pnt(Npnt,3)-pnt(1,3));

SegPoint_Set1 = [];
index1 = [];
j = 7;
while j <= Npnt-7
    if Smooth_Speed(j,1) < 0.7*Average_Speed                                            % Speed based segment points
        index1 = [index1; j];
        SegPoint_Set1 = [SegPoint_Set1; [pnt(j,1) pnt(j,2) pnt(j,3)]];
        j = j+7;
    end    
    j = j+1;
end
SegPoint_Set1 = [[pnt(1,1) pnt(1,2) pnt(1,3)]; SegPoint_Set1; [pnt(Npnt,1) pnt(Npnt,2) pnt(Npnt,3)]];
index1 = [1; index1; Npnt];
%plot(pnt(1:end,1),pnt(1:end,2), SegPoint_Set1(:,1),SegPoint_Set1(:,2), 'ro')
%hold on

AB = zeros(Npnt,2);
abc = zeros(Npnt,3);
error_line = zeros(Npnt,1);
error_circle = zeros(Npnt,1);
arclength = zeros(Npnt,1);
C_line(:,1) = zeros(Npnt,1);
C_circle(:,1) = zeros(Npnt,1);
for i = 6:Npnt-5                                                                        % Curvature and error finding
    AB(i,:) = (inv([sum(pnt(i-5:i+5,1).^2), sum(pnt(i-5:i+5,1)); sum(pnt(i-5:i+5,1)), 11])*[sum(pnt(i-5:i+5,1).*pnt(i-5:i+5,2)); sum(pnt(i-5:i+5,2))])';
    abc(i,:) = (inv([2*sum(pnt(i-5:i+5,1).^2), 2*sum(pnt(i-5:i+5,1).*pnt(i-5:i+5,2)), sum(pnt(i-5:i+5,1)); 2*sum(pnt(i-5:i+5,1).*pnt(i-5:i+5,2)), 2*sum(pnt(i-5:i+5,2).^2), sum(pnt(i-5:i+5,2)); 2*sum(pnt(i-5:i+5,1)), 2*sum(pnt(i-5:i+5,2)), 11])*[-sum((pnt(i-5:i+5,1).^2 + pnt(i-5:i+5,2).^2).*pnt(i-5:i+5,1)); -sum((pnt(i-5:i+5,1).^2 + pnt(i-5:i+5,2).^2).*pnt(i-5:i+5,2)); -sum(pnt(i-5:i+5,1).^2 + pnt(i-5:i+5,2).^2)])';
    error_line(i,1) = (1/11)*sum(abs(AB(i,1)*pnt(i-5:i+5,1) + AB(i,2) - pnt(i-5:i+5,2)));
    error_circle(i,1) = (1/11)*sum(abs(sqrt((pnt(i-5:i+5,1) + abc(i,1)).^2 + (pnt(i-5:i+5,2) + abc(i,2)).^2) - sqrt(abc(i,1)^2 + abc(i,2)^2 - abc(i,3))));
    arclength(i,1) = dist(i+5) - dist(i-4);
    C_line(i,1) = abs(diff(atan2d(AB(i-1:i,1),1))./ diff(arclength(i-1:i,1)));
    C_circle(i,1) = abs(diff(atan2d(pnt(i-1:i,2)+abc(i-1:i,2), pnt(i-1:i,1)+abc(i-1:i,1)))./ diff(arclength(i-1:i,1)));
end

SegPoint_Set2 = [];
SegPoint_Set3 = [];
index2 = [];
index3 = [];
k = 1;
while k <= Npnt-7                                                                       % Curvature based segment points
    if C_line(k,1) > 1e5 && C_line(k,1) < 1e+8
        index2 = [index2; k];
        SegPoint_Set2 = [SegPoint_Set2; [pnt(k,1) pnt(k,2) pnt(k,3)]];
        k = k+7;
        %plot(SegPoint_Set2(:,1), SegPoint_Set2(:,2), 'go')
        %hold on
    elseif k < Npnt-20 && sign(diff(atan2d(AB(k+4:k+5,1),1))) == -sign(diff(atan2d(AB(k+5:k+6,1),1)))  
        sgn = sign(diff(atan2d(AB(k+5:k+6,1),1)));                                      % Curvature sign based segment points
        flag = 0;
        for j = k+6:k+15
            if sign(diff(atan2d(AB(j:j+1,1),1))) == sgn
                flag = flag + 1;
            end
            if flag >=10
                index3 = [index3; k];
                SegPoint_Set3 = [SegPoint_Set3; [pnt(k+5,1) pnt(k+5,2) pnt(k+5,3)]];
                k = k+7;
                %plot(SegPoint_Set3(:,1), SegPoint_Set3(:,2), 'ko')
                %hold on
            end
        end
    end
    k = k+1;
end

index = sort([index1; index2; index3]);
dif = diff(index)
for n = 1:length(dif)
    if dif(n) < 7
        index(n) = [];
    end
end
SegPoints = pnt(index,:);
%plot(pnt(index,1),pnt(index,2), 'b*')
SegPoints = pnt(index,:);

for i = 1:length(SegPoints(:,1))-1                                                       % Segmenting
    A = (SegPoints(i+1,2)-SegPoints(i,2))/(SegPoints(i+1,1)-SegPoints(i,1));
    B = SegPoints(i,2) - A*SegPoints(i,1);
    error_segline(i,1) = (1/(index(i+1)-index(i)+1))*sum(abs(A*pnt(index(i):index(i+1),1) + B - pnt(index(i):index(i+1),2)));
    seg_arclength(i,1) = dist(index(i+1)) - dist(index(i));
    if error_segline(i,1) < 0.1*seg_arclength(i,1)
        %fprintf('line\n')
        plot([SegPoints(i,1) SegPoints(i+1,1)],[SegPoints(i,2) SegPoints(i+1,2)],'-r','LineWidth',2)
        hold on
    elseif error_segline(i,1) >= 0.1*seg_arclength(i,1)
        %fprintf('circle\n')
        abc = (inv([2*sum(pnt(index(i):index(i+1),1).^2), 2*sum(pnt(index(i):index(i+1),1).*pnt(index(i):index(i+1),2)),...
            sum(pnt(index(i):index(i+1),1)); 2*sum(pnt(index(i):index(i+1),1).*pnt(index(i):index(i+1),2)), 2*sum(pnt(index(i):index(i+1),2).^2),...
            sum(pnt(index(i):index(i+1),2)); 2*sum(pnt(index(i):index(i+1),1)), 2*sum(pnt(index(i):index(i+1),2)), (index(i+1)-index(i)+1)])...
            *[-sum((pnt(index(i):index(i+1),1).^2 + pnt(index(i):index(i+1),2).^2).*pnt(index(i):index(i+1),1)); -sum((pnt(index(i):index(i+1),1).^2 ...
            + pnt(index(i):index(i+1),2).^2).*pnt(index(i):index(i+1),2)); -sum(pnt(index(i):index(i+1),1).^2 + pnt(index(i):index(i+1),2).^2)])';       
                    
        th1 = atan2d((pnt(index(i),2)+abc(2)), (pnt(index(i),1)+abc(1)));
        th2 = atan2d((pnt(index(i+1),2)+abc(2)), (pnt(index(i+1),1)+abc(1)));
        th_test = atan2d((pnt(index(i)+1,2)+abc(2)), (pnt(index(i)+1,1)+abc(1)));
        
        if th1 < 0 
            th1 = th1 + 360;
        end
        if th2 < 0
            th2 = th2 + 360;
        end
        if th_test < 0
            th_test = th_test + 360;
        end
        if th_test > th1 && th2 < th1
            th1 = th1 - 360;
        end
        if th_test < th1 && th2 > th1
            th2 = th2 - 360;
        end
        if abs(th2 - th1) > 340
            th1 = 0;
            th2 = 360;
        end
        plot(-abc(1)+sqrt(abc(1)^2 + abc(2)^2 - abc(3))*cosd(linspace(th1,th2,50)), -abc(2)+sqrt(abc(1)^2 + abc(2)^2 - abc(3))*sind(linspace(th1,th2,50)),'-r','LineWidth',2);
    end

end

pnt = [];
Npnt = 1;
