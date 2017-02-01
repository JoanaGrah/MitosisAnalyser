function varargout = examine(varargin)
% EXAMINE MATLAB code for examine.fig
%      EXAMINE, by itself, creates a new EXAMINE or raises the existing
%      singleton*.
%
%      H = EXAMINE returns the handle to a new EXAMINE or the handle to
%      the existing singleton*.
%
%      EXAMINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXAMINE.M with the given input arguments.
%
%      EXAMINE('Property','Value',...) creates a new EXAMINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before examine_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to examine_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help examine

% Last Modified by GUIDE v2.5 18-Aug-2015 10:53:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @examine_OpeningFcn, ...
                   'gui_OutputFcn',  @examine_OutputFcn, ...
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


% --- Executes just before examine is made visible.
function examine_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to examine (see VARARGIN)

% Choose default command line output for examine
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes examine wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Get parameters from main MiA GUI
hMiA = getappdata(0,'hMiA');

handles.img = getappdata(hMiA,'img');
handles.imageName = getappdata(hMiA,'imageName');
handles.imageSize = getappdata(hMiA,'imageSize');
handles.pixelSize = getappdata(hMiA,'pixelSize');

axes(handles.axesImage), imshow(handles.img);

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = examine_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonCrop.
function pushbuttonCrop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imcrop(handles.axesImage);


% --- Executes on button press in pushbuttonMeasureDistance.
function pushbuttonMeasureDistance_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMeasureDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imdistline(handles.axesImage);


% --- Executes on button press in pushbuttonAddScalebar.
function pushbuttonAddScalebar_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddScalebar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

scalebarWidth = handles.pixelSize*100;
m = handles.imageSize(1);
n = handles.imageSize(2);
axes(handles.axesImage),
line([n-80,n-180],[m-80,m-80],'Color','m','LineWidth',3)
line([n-80,n-80],[m-90,m-70],'Color','m','LineWidth',3)
line([n-180,n-180],[m-90,m-70],'Color','m','LineWidth',3)
axes(handles.axesImage),
text(n-170,m-50,[num2str(scalebarWidth),' \mum'],'Color','m','FontSize',14,'FontWeight','bold')


% --- Executes on button press in pushbuttonContrastEnhancement.
function pushbuttonContrastEnhancement_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonContrastEnhancement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imcontrast(handles.axesImage);


% --- Executes on button press in pushbuttonHistogramEqualisation.
function pushbuttonHistogramEqualisation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonHistogramEqualisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

imgHistEq = histeq(handles.img);
axes(handles.axesImage), imshow(imgHistEq);


% --- Executes on button press in pushbuttonAdaptiveHistogramEqualisation.
function pushbuttonAdaptiveHistogramEqualisation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdaptiveHistogramEqualisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

imgAdaptHistEq = adapthisteq(handles.img);
axes(handles.axesImage), imshow(imgAdaptHistEq);


% --- Executes on button press in pushbuttonClear.
function pushbuttonClear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axesImage), imshow(handles.img);


% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Please enter a file name for the image to be saved:'};
dlgTitle = 'File name';
numLines = 1;
defAns = {[handles.imageName,'_edited']};
answer = inputdlg(prompt,dlgTitle,numLines,defAns);
frame = getframe(gca);
saveImage = frame.cdata;
imwrite(saveImage,[answer{1},'.png']);
disp(['Image "',answer{1},'.png" saved.'])
