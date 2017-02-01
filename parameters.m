function varargout = parameters(varargin)
% PARAMETERS MATLAB code for parameters.fig
%      PARAMETERS, by itself, creates a new PARAMETERS or raises the existing
%      singleton*.
%
%      H = PARAMETERS returns the handle to a new PARAMETERS or the handle to
%      the existing singleton*.
%
%      PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETERS.M with the given input arguments.
%
%      PARAMETERS('Property','Value',...) creates a new PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help parameters

% Last Modified by GUIDE v2.5 13-Aug-2015 14:25:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @parameters_OutputFcn, ...
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


% --- Executes just before parameters is made visible.
function parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to parameters (see VARARGIN)

% Choose default command line output for parameters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes parameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Get parameters from main MiA GUI
hMiA = getappdata(0,'hMiA');

handles.nameDataSet = getappdata(hMiA,'nameDataSet');
handles.images = getappdata(hMiA,'images');
handles.currentImage = getappdata(hMiA,'currentImage');
handles.pixelSize = getappdata(hMiA,'pixelSize');
handles.timeResolution = getappdata(hMiA,'timeResolution');
handles.parameters.mitosisThreshold = getappdata(hMiA,'mitosisThreshold');
handles.parameters.radiusMin = getappdata(hMiA,'radiusMin');
handles.parameters.radiusMax = getappdata(hMiA,'radiusMax');
handles.parameters.sensitivity = getappdata(hMiA,'sensitivity');
handles.parameters.lambda1 = getappdata(hMiA,'lambda1');
handles.parameters.lambda2 = getappdata(hMiA,'lambda2');
handles.parameters.mu = getappdata(hMiA,'mu');
handles.parameters.nu = getappdata(hMiA,'nu');
handles.parameters.g_adjust_low = getappdata(hMiA,'g_adjust_low');
handles.parameters.g_adjust_high = getappdata(hMiA,'g_adjust_high');
handles.parameters.omega = getappdata(hMiA,'omega');
handles.parameters.timeStep = getappdata(hMiA,'timeStep');
handles.parameters.maxIterations = getappdata(hMiA,'maxIterations');
handles.parameters.phiUpdate = getappdata(hMiA,'phiUpdate');
handles.parameters.epsilonNormGradReg = getappdata(hMiA,'epsilonNormGradReg');
handles.parameters.epsilonDeltaReg = getappdata(hMiA,'epsilonDeltaReg');

set(handles.textNameDataSet,'String',handles.nameDataSet);
set(handles.editPixelSize,'String',num2str(handles.pixelSize));
set(handles.editTimeRes,'String',num2str(handles.timeResolution));
set(handles.editMitosisThreshold,'String',num2str(handles.parameters.mitosisThreshold));
set(handles.editRadiusMin,'String',num2str(handles.parameters.radiusMin));
set(handles.editRadiusMax,'String',num2str(handles.parameters.radiusMax));
set(handles.editSensitivity,'String',num2str(handles.parameters.sensitivity));
set(handles.editLambda1,'String',num2str(handles.parameters.lambda1));
set(handles.editLambda2,'String',num2str(handles.parameters.lambda2));
set(handles.editMu,'String',num2str(handles.parameters.mu));
set(handles.editNu,'String',num2str(handles.parameters.nu));
set(handles.editGAdjustLow,'String',num2str(handles.parameters.g_adjust_low));
set(handles.editGAdjustHigh,'String',num2str(handles.parameters.g_adjust_high));
set(handles.editOmega,'String',num2str(handles.parameters.omega));
set(handles.editTimeStep,'String',num2str(handles.parameters.timeStep));
set(handles.editMaxIterations,'String',num2str(handles.parameters.maxIterations));
set(handles.editPhiUpdate,'String',num2str(handles.parameters.phiUpdate));
set(handles.editEpsNormGradReg,'String',num2str(handles.parameters.epsilonNormGradReg));
set(handles.editEpsDeltaReg,'String',num2str(handles.parameters.epsilonDeltaReg));

% Create edge detector function g
img = handles.images(:,:,handles.currentImage);
if verLessThan('matlab','8.5')
    gaussfilt = fspecial('gaussian',5,1);
    smoothed = imfilter(img,gaussfilt,'replicate','conv');
else
    smoothed = imgaussfilt(img,1,'FilterSize',5);
end
locstd = stdfilt(smoothed);
scaled = (locstd - min(locstd(:))) ./ (max(locstd(:)) - min(locstd(:)));
g = imadjust(scaled,[handles.parameters.g_adjust_low handles.parameters.g_adjust_high],[1 0],3);

axes(handles.axesG);
imagesc(g);axis off;axis image;colorbar;title('g');

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA SET PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbuttonOpen.
function pushbuttonOpen_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

listing = dir('presetDataSets/*.txt');
numberOfPresets = numel(listing);
presets = cell(numberOfPresets,1);
displayPresets = cell(numberOfPresets,1);
for i = 1:numberOfPresets
    presets{i} = listing(i).name;
    [~,displayPresets{i},~] = fileparts(listing(i).name);
end
positionDefault = find(~isempty(strfind(displayPresets,'NEW_DATA_SET')));
if positionDefault~=1
    presets(positionDefault,:)=[];
    displayPresets(positionDefault,:)=[];
    presets_new = cell(numberOfPresets,1);
    displayPresets_new = cell(numberOfPresets,1);
    presets_new{1} = 'NEW_DATA_SET.txt';
    displayPresets_new{1} = 'NEW_DATA_SET';
    presets_new(2:numberOfPresets) = presets;
    displayPresets_new(2:numberOfPresets) = displayPresets;
    clear presets displayPresets
    presets = presets_new;
    displayPresets = displayPresets_new;
end
[nrSelectedPreset,ok] = listdlg('ListString',displayPresets,...
                              'SelectionMode','single',...
                              'ListSize',[200 150],...
                              'Name','Data Set',...
                              'PromptString','Please select your type of data set.');
if ok==0
    errordlg('No data set selected.');
    return
end

fileID = fopen(['presetDataSets/',presets{nrSelectedPreset}]);
selectedPreset = textscan(fileID,'%s');
fclose(fileID);

handles.nameDataSet                   = displayPresets(nrSelectedPreset);
handles.parameters.lambda1            = str2double(selectedPreset{1}{1});
handles.parameters.lambda2            = str2double(selectedPreset{1}{2});
handles.parameters.mu                 = str2double(selectedPreset{1}{3});
handles.parameters.nu                 = str2double(selectedPreset{1}{4});
handles.parameters.g_adjust_low       = str2double(selectedPreset{1}{5});
handles.parameters.g_adjust_high      = str2double(selectedPreset{1}{6});
handles.parameters.omega              = str2double(selectedPreset{1}{7});
handles.parameters.timeStep           = str2double(selectedPreset{1}{8});
handles.parameters.maxIterations      = str2double(selectedPreset{1}{9});
handles.parameters.phiUpdate          = str2double(selectedPreset{1}{10});
handles.parameters.epsilonNormGradReg = str2double(selectedPreset{1}{11});
handles.parameters.epsilonDeltaReg    = str2double(selectedPreset{1}{12});
handles.parameters.mitosisThreshold   = str2double(selectedPreset{1}{13});
handles.parameters.radiusMin          = str2double(selectedPreset{1}{14});
handles.parameters.radiusMax          = str2double(selectedPreset{1}{15});
handles.parameters.sensitivity        = str2double(selectedPreset{1}{16});
handles.pixelSize          = str2double(selectedPreset{1}{17});
handles.timeResolution     = str2double(selectedPreset{1}{18});

set(handles.textNameDataSet,'String',handles.nameDataSet);
set(handles.editPixelSize,'String',num2str(handles.pixelSize));
set(handles.editTimeRes,'String',num2str(handles.timeResolution));
set(handles.editMitosisThreshold,'String',num2str(handles.parameters.mitosisThreshold));
set(handles.editRadiusMin,'String',num2str(handles.parameters.radiusMin));
set(handles.editRadiusMax,'String',num2str(handles.parameters.radiusMax));
set(handles.editSensitivity,'String',num2str(handles.parameters.sensitivity));
set(handles.editLambda1,'String',num2str(handles.parameters.lambda1));
set(handles.editLambda2,'String',num2str(handles.parameters.lambda2));
set(handles.editMu,'String',num2str(handles.parameters.mu));
set(handles.editNu,'String',num2str(handles.parameters.nu));
set(handles.editGAdjustLow,'String',num2str(handles.parameters.g_adjust_low));
set(handles.editGAdjustHigh,'String',num2str(handles.parameters.g_adjust_high));
set(handles.editOmega,'String',num2str(handles.parameters.omega));
set(handles.editTimeStep,'String',num2str(handles.parameters.timeStep));
set(handles.editMaxIterations,'String',num2str(handles.parameters.maxIterations));
set(handles.editPhiUpdate,'String',num2str(handles.parameters.phiUpdate));
set(handles.editEpsNormGradReg,'String',num2str(handles.parameters.epsilonNormGradReg));
set(handles.editEpsDeltaReg,'String',num2str(handles.parameters.epsilonDeltaReg));

% Create edge detector function g
img = handles.images(:,:,handles.currentImage);
if verLessThan('matlab','8.5')
    gaussfilt = fspecial('gaussian',5,1);
    smoothed = imfilter(img,gaussfilt,'replicate','conv');
else
    smoothed = imgaussfilt(img,1,'FilterSize',5);
end
locstd = stdfilt(smoothed);
scaled = (locstd - min(locstd(:))) ./ (max(locstd(:)) - min(locstd(:)));
g = imadjust(scaled,[handles.parameters.g_adjust_low handles.parameters.g_adjust_high],[1 0],3);

axes(handles.axesG);
imagesc(g);axis off;axis image;colorbar;title('g');

guidata(hObject, handles);


% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

handles.nameDataSet = get(handles.textNameDataSet,'String');
handles.pixelSize = str2double(get(handles.editPixelSize,'String'));
handles.timeResolution = str2double(get(handles.editTimeRes,'String'));
handles.parameters.mitosisThreshold = str2double(get(handles.editMitosisThreshold,'String'));
handles.parameters.radiusMin = str2double(get(handles.editRadiusMin,'String'));
handles.parameters.radiusMax = str2double(get(handles.editRadiusMax,'String'));
handles.parameters.sensitivity = str2double(get(handles.editSensitivity,'String'));
handles.parameters.lambda1 = str2double(get(handles.editLambda1,'String'));
handles.parameters.lambda2 = str2double(get(handles.editLambda2,'String'));
handles.parameters.mu = str2double(get(handles.editMu,'String'));
handles.parameters.nu = str2double(get(handles.editNu,'String'));
handles.parameters.g_adjust_low = str2double(get(handles.editGAdjustLow,'String'));
handles.parameters.g_adjust_high = str2double(get(handles.editGAdjustHigh,'String'));
handles.parameters.omega = str2double(get(handles.editOmega,'String'));
handles.parameters.timeStep = str2double(get(handles.editTimeStep,'String'));
handles.parameters.maxIterations = str2double(get(handles.editMaxIterations,'String'));
handles.parameters.phiUpdate = str2double(get(handles.editPhiUpdate,'String'));
handles.parameters.epsilonNormGradReg = str2double(get(handles.editEpsNormGradReg,'String'));
handles.parameters.epsilonDeltaReg = str2double(get(handles.editEpsDeltaReg,'String'));

hMiA = getappdata(0,'hMiA');

setappdata(hMiA,'nameDataSet',handles.nameDataSet);
setappdata(hMiA,'pixelSize',handles.pixelSize);
setappdata(hMiA,'timeResolution',handles.timeResolution);
setappdata(hMiA,'mitosisThreshold',handles.parameters.mitosisThreshold);
setappdata(hMiA,'radiusMin',handles.parameters.radiusMin);
setappdata(hMiA,'radiusMax',handles.parameters.radiusMax);
setappdata(hMiA,'sensitivity',handles.parameters.sensitivity);
setappdata(hMiA,'lambda1',handles.parameters.lambda1);
setappdata(hMiA,'lambda2',handles.parameters.lambda2);
setappdata(hMiA,'mu',handles.parameters.mu);
setappdata(hMiA,'nu',handles.parameters.nu);
setappdata(hMiA,'g_adjust_low',handles.parameters.g_adjust_low);
setappdata(hMiA,'g_adjust_high',handles.parameters.g_adjust_high);
setappdata(hMiA,'omega',handles.parameters.omega);
setappdata(hMiA,'timeStep',handles.parameters.timeStep);
setappdata(hMiA,'maxIterations',handles.parameters.maxIterations);
setappdata(hMiA,'phiUpdate',handles.parameters.phiUpdate);
setappdata(hMiA,'epsilonNormGradReg',handles.parameters.epsilonNormGradReg);
setappdata(hMiA,'epsilonDeltaReg',handles.parameters.epsilonDeltaReg);

contentTextFile = [handles.parameters.lambda1;
                   handles.parameters.lambda2;
                   handles.parameters.mu;
                   handles.parameters.nu;
                   handles.parameters.g_adjust_low;
                   handles.parameters.g_adjust_high;
                   handles.parameters.omega;
                   handles.parameters.timeStep;
                   handles.parameters.maxIterations;
                   handles.parameters.phiUpdate;
                   handles.parameters.epsilonNormGradReg;
                   handles.parameters.epsilonDeltaReg;
                   handles.parameters.mitosisThreshold;
                   handles.parameters.radiusMin;
                   handles.parameters.radiusMax;
                   handles.parameters.sensitivity;
                   handles.pixelSize;
                   handles.timeResolution];

prompt = {'Please enter a file name:'};
dlgTitle = 'File name';
numLines = 1;
answer = inputdlg(prompt,dlgTitle,numLines);
fileName = ['presetDataSets/',answer{1},'.txt'];
fileID = fopen(fileName,'w');
formatSpec = '%f \n';
fprintf(fileID,formatSpec,contentTextFile);
fclose(fileID);

guidata(hObject, handles);


% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

if strcmp(handles.nameDataSet,'NEW_DATA_SET')
    errordlg('Default values in NEW_DATA_SET cannot be changed! Please use the ''Save as...'' button and rename your data set.');
    return
end

handles.nameDataSet = get(handles.textNameDataSet,'String');
handles.pixelSize = str2double(get(handles.editPixelSize,'String'));
handles.timeResolution = str2double(get(handles.editTimeRes,'String'));
handles.parameters.mitosisThreshold = str2double(get(handles.editMitosisThreshold,'String'));
handles.parameters.radiusMin = str2double(get(handles.editRadiusMin,'String'));
handles.parameters.radiusMax = str2double(get(handles.editRadiusMax,'String'));
handles.parameters.sensitivity = str2double(get(handles.editSensitivity,'String'));
handles.parameters.lambda1 = str2double(get(handles.editLambda1,'String'));
handles.parameters.lambda2 = str2double(get(handles.editLambda2,'String'));
handles.parameters.mu = str2double(get(handles.editMu,'String'));
handles.parameters.nu = str2double(get(handles.editNu,'String'));
handles.parameters.g_adjust_low = str2double(get(handles.editGAdjustLow,'String'));
handles.parameters.g_adjust_high = str2double(get(handles.editGAdjustHigh,'String'));
handles.parameters.omega = str2double(get(handles.editOmega,'String'));
handles.parameters.timeStep = str2double(get(handles.editTimeStep,'String'));
handles.parameters.maxIterations = str2double(get(handles.editMaxIterations,'String'));
handles.parameters.phiUpdate = str2double(get(handles.editPhiUpdate,'String'));
handles.parameters.epsilonNormGradReg = str2double(get(handles.editEpsNormGradReg,'String'));
handles.parameters.epsilonDeltaReg = str2double(get(handles.editEpsDeltaReg,'String'));

hMiA = getappdata(0,'hMiA');

setappdata(hMiA,'nameDataSet',handles.nameDataSet);
setappdata(hMiA,'pixelSize',handles.pixelSize);
setappdata(hMiA,'timeResolution',handles.timeResolution);
setappdata(hMiA,'mitosisThreshold',handles.parameters.mitosisThreshold);
setappdata(hMiA,'radiusMin',handles.parameters.radiusMin);
setappdata(hMiA,'radiusMax',handles.parameters.radiusMax);
setappdata(hMiA,'sensitivity',handles.parameters.sensitivity);
setappdata(hMiA,'lambda1',handles.parameters.lambda1);
setappdata(hMiA,'lambda2',handles.parameters.lambda2);
setappdata(hMiA,'mu',handles.parameters.mu);
setappdata(hMiA,'nu',handles.parameters.nu);
setappdata(hMiA,'g_adjust_low',handles.parameters.g_adjust_low);
setappdata(hMiA,'g_adjust_high',handles.parameters.g_adjust_high);
setappdata(hMiA,'omega',handles.parameters.omega);
setappdata(hMiA,'timeStep',handles.parameters.timeStep);
setappdata(hMiA,'maxIterations',handles.parameters.maxIterations);
setappdata(hMiA,'phiUpdate',handles.parameters.phiUpdate);
setappdata(hMiA,'epsilonNormGradReg',handles.parameters.epsilonNormGradReg);
setappdata(hMiA,'epsilonDeltaReg',handles.parameters.epsilonDeltaReg);

contentTextFile = [handles.parameters.lambda1;
                   handles.parameters.lambda2;
                   handles.parameters.mu;
                   handles.parameters.nu;
                   handles.parameters.g_adjust_low;
                   handles.parameters.g_adjust_high;
                   handles.parameters.omega;
                   handles.parameters.timeStep;
                   handles.parameters.maxIterations;
                   handles.parameters.phiUpdate;
                   handles.parameters.epsilonNormGradReg;
                   handles.parameters.epsilonDeltaReg;
                   handles.parameters.mitosisThreshold;
                   handles.parameters.radiusMin;
                   handles.parameters.radiusMax;
                   handles.parameters.sensitivity;
                   handles.pixelSize;
                   handles.timeResolution];

fileName = ['presetDataSets/',char(handles.nameDataSet),'.txt'];
fileID = fopen(fileName,'w');
formatSpec = '%f \n';
fprintf(fileID,formatSpec,contentTextFile);
fclose(fileID);

guidata(hObject, handles);


function editTimeRes_Callback(hObject, eventdata, handles)
% hObject    handle to editTimeRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTimeRes as text
%        str2double(get(hObject,'String')) returns contents of editTimeRes as a double


% --- Executes during object creation, after setting all properties.
function editTimeRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTimeRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editPixelSize_Callback(hObject, eventdata, handles)
% hObject    handle to editPixelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPixelSize as text
%        str2double(get(hObject,'String')) returns contents of editPixelSize as a double


% --- Executes during object creation, after setting all properties.
function editPixelSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPixelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PARAMETERS MITOSIS DETECTION PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editMitosisThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editMitosisThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMitosisThreshold as text
%        str2double(get(hObject,'String')) returns contents of editMitosisThreshold as a double


% --- Executes during object creation, after setting all properties.
function editMitosisThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMitosisThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editRadiusMin_Callback(hObject, eventdata, handles)
% hObject    handle to editRadiusMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRadiusMin as text
%        str2double(get(hObject,'String')) returns contents of editRadiusMin as a double


% --- Executes during object creation, after setting all properties.
function editRadiusMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRadiusMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editRadiusMax_Callback(hObject, eventdata, handles)
% hObject    handle to editRadiusMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRadiusMax as text
%        str2double(get(hObject,'String')) returns contents of editRadiusMax as a double


% --- Executes during object creation, after setting all properties.
function editRadiusMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRadiusMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSensitivity_Callback(hObject, eventdata, handles)
% hObject    handle to editSensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSensitivity as text
%        str2double(get(hObject,'String')) returns contents of editSensitivity as a double


% --- Executes during object creation, after setting all properties.
function editSensitivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PARAMETERS CELL TRACKING PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editLambda1_Callback(hObject, eventdata, handles)
% hObject    handle to editLambda1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLambda1 as text
%        str2double(get(hObject,'String')) returns contents of editLambda1 as a double


% --- Executes during object creation, after setting all properties.
function editLambda1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLambda1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editLambda2_Callback(hObject, eventdata, handles)
% hObject    handle to editLambda2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLambda2 as text
%        str2double(get(hObject,'String')) returns contents of editLambda2 as a double


% --- Executes during object creation, after setting all properties.
function editLambda2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLambda2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editMu_Callback(hObject, eventdata, handles)
% hObject    handle to editMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMu as text
%        str2double(get(hObject,'String')) returns contents of editMu as a double


% --- Executes during object creation, after setting all properties.
function editMu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editNu_Callback(hObject, eventdata, handles)
% hObject    handle to editNu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNu as text
%        str2double(get(hObject,'String')) returns contents of editNu as a double


% --- Executes during object creation, after setting all properties.
function editNu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editGAdjustLow_Callback(hObject, eventdata, handles)
% hObject    handle to editGAdjustLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGAdjustLow as text
%        str2double(get(hObject,'String')) returns contents of editGAdjustLow as a double

guidata(hObject, handles);

g_adjust_low = str2double(get(hObject,'String'));
g_adjust_high = str2double(get(handles.editGAdjustHigh,'String'));

img = handles.images(:,:,handles.currentImage);
if verLessThan('matlab','8.5')
    gaussfilt = fspecial('gaussian',5,1);
    smoothed = imfilter(img,gaussfilt,'replicate','conv');
else
    smoothed = imgaussfilt(img,1,'FilterSize',5);
end
locstd = stdfilt(smoothed);
scaled = (locstd - min(locstd(:))) ./ (max(locstd(:)) - min(locstd(:)));
g = imadjust(scaled,[g_adjust_low g_adjust_high],[1 0],3);

axes(handles.axesG);
imagesc(g);axis off;axis image;colorbar;title('g');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editGAdjustLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGAdjustLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editGAdjustHigh_Callback(hObject, eventdata, handles)
% hObject    handle to editGAdjustHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGAdjustHigh as text
%        str2double(get(hObject,'String')) returns contents of editGAdjustHigh as a double

guidata(hObject, handles);

g_adjust_low = str2double(get(handles.editGAdjustLow,'String'));
g_adjust_high = str2double(get(hObject,'String'));

img = handles.images(:,:,handles.currentImage);
if verLessThan('matlab','8.5')
    gaussfilt = fspecial('gaussian',5,1);
    smoothed = imfilter(img,gaussfilt,'replicate','conv');
else
    smoothed = imgaussfilt(img,1,'FilterSize',5);
end
locstd = stdfilt(smoothed);
scaled = (locstd - min(locstd(:))) ./ (max(locstd(:)) - min(locstd(:)));
g = imadjust(scaled,[g_adjust_low g_adjust_high],[1 0],3);

axes(handles.axesG);
imagesc(g);axis off;axis image;colorbar;title('g');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editGAdjustHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGAdjustHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOmega_Callback(hObject, eventdata, handles)
% hObject    handle to editOmega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOmega as text
%        str2double(get(hObject,'String')) returns contents of editOmega as a double


% --- Executes during object creation, after setting all properties.
function editOmega_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOmega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editTimeStep_Callback(hObject, eventdata, handles)
% hObject    handle to editTimeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTimeStep as text
%        str2double(get(hObject,'String')) returns contents of editTimeStep as a double


% --- Executes during object creation, after setting all properties.
function editTimeStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTimeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editMaxIterations_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxIterations as text
%        str2double(get(hObject,'String')) returns contents of editMaxIterations as a double


% --- Executes during object creation, after setting all properties.
function editMaxIterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editPhiUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to editPhiUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPhiUpdate as text
%        str2double(get(hObject,'String')) returns contents of editPhiUpdate as a double


% --- Executes during object creation, after setting all properties.
function editPhiUpdate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPhiUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editEpsNormGradReg_Callback(hObject, eventdata, handles)
% hObject    handle to editEpsNormGradReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEpsNormGradReg as text
%        str2double(get(hObject,'String')) returns contents of editEpsNormGradReg as a double


% --- Executes during object creation, after setting all properties.
function editEpsNormGradReg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEpsNormGradReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editEpsDeltaReg_Callback(hObject, eventdata, handles)
% hObject    handle to editEpsDeltaReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEpsDeltaReg as text
%        str2double(get(hObject,'String')) returns contents of editEpsDeltaReg as a double


% --- Executes during object creation, after setting all properties.
function editEpsDeltaReg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEpsDeltaReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
