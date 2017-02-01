function varargout = MiA(varargin)
% MIA MATLAB code for MiA.fig
%      MIA, by itself, creates a new MIA or raises the existing
%      singleton*.
%
%      H = MIA returns the handle to a new MIA or the handle to
%      the existing singleton*.
%
%      MIA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MIA.M with the given input arguments.
%
%      MIA('Property','Value',...) creates a new MIA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MiA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MiA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MiA

% Last Modified by GUIDE v2.5 12-Jul-2015 14:55:12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Code for main GUI "MiA", created with MATLAB, R2015a                   %
%  Required files: examine.fig, examine.m, parameters.fig, parameters.m,  %
%  statistics.fig, statistics.m, mitosisDetection.m,                      %
%  backwardsTracking.m, forwardsTracking.m, normalVelocity.m,             %
%  createbandmapping2D.m, topologyCheck2D.m, getLevelSetFromCoords.m,     %
%  calculateStatistics.m                                                  %
%  Required sub-folders: presetDataSets, results                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  See "documentation.pdf" for detailed explanation                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Copyright: Joana Grah, jg704@cam.ac.uk                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MiA_OpeningFcn, ...
                   'gui_OutputFcn',  @MiA_OutputFcn, ...
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


% --- Executes just before MiA is made visible.
function MiA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MiA (see VARARGIN)

% Choose default command line output for MiA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MiA wait for user response (see UIRESUME)
% uiwait(handles.figureMiA);

setappdata(0,'hMiA',gcf);

% OPTION 1: Select individual images within a folder
% [fileName,pathName] = uigetfile({'*.tif';'*.tiff'},...
%     'Please select the sequence of images to be analysed.',...
%     'MultiSelect','on');
% numberOfFrames = numel(fileName);
% OPTION 2: Select a whole folder of images
pathName = uigetdir('',...
    'Please select the folder with the image sequence to be analysed.');
searchFiles = strcat(pathName,'/*.tif');
listFiles = dir(searchFiles);
numberOfFrames = numel(listFiles);
fileName = cell(numberOfFrames,1);
for i = 1:numberOfFrames
    fileName{i} = listFiles(i).name;
end
if numberOfFrames==0
    errordlg('No images selected.');
    return
end

% Choose type of data set
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

disp('Loading data...');

% Create waitbar
wb1 = waitbar(0, 'Loading data... Please wait.');
set(findobj(wb1,'type','patch'),'edgecolor','k','facecolor','b');

% Read text file with preset parameters
fileID = fopen(['presetDataSets/',presets{nrSelectedPreset}]);
selectedPreset = textscan(fileID,'%s');
fclose(fileID);

% Read image sequence
image1 = imread([pathName,'/',fileName{1}]);
m = size(image1,1);
n = size(image1,2);
images = zeros(m,n,numberOfFrames);
for i = 1:numberOfFrames
    waitbar(i/numberOfFrames);
    drawnow;
    imgtmp = imread([pathName,'/',fileName{i}]);
    if ismatrix(imgtmp)
        images(:,:,i) = im2double(imgtmp);
    elseif ndims(imgtmp)==3
        images(:,:,i) = im2double(rgb2gray(imgtmp));
    else
        errordlg('Image dimensions are too small or too big (should be 2D or 3D).');
    end
end
clear imgtmp

set(handles.sliderImages, 'Min',        1, ...
                          'Max',        numberOfFrames, ...
                          'Value',      1, ...
                          'SliderStep', [1 5]/(numberOfFrames-1));

axes(handles.axesImages), imshow(images(:,:,1));
set(handles.textData,     'String', fileName{1});
set(handles.textImgTotal, 'String', ['  / ', num2str(numberOfFrames)]);
set(handles.editImgNr,    'String', '1');

if verLessThan('matlab','8.3')
    disp('MATLAB version is less than R2014a. GPU computation not possible.');
    handles.gpu = 0;
else
    try
        gpuArray(0);
        handles.gpu = 1;
        disp('GPU computation with CUDA initialised...');
    catch err
        handles.gpu = 0;
        disp('GPU computation not possible.');
        disp(err.message);
    end
end

handles.pathName = pathName;
handles.fileName = fileName;
handles.numberOfFrames = numberOfFrames;
handles.imageSize = [m,n];
handles.images = images;
handles.currentImage = 1;
handles.presets = presets;
handles.displayPresets = displayPresets;
handles.selectedPreset = selectedPreset;
handles.nameDataSet = displayPresets(nrSelectedPreset);
handles.mitosisDetectionDone = 0;
handles.cellTrackingDone = 0;
handles.displayMitosisDetectionResults = 0;
handles.displayCellTrackingResults = 0;

% Set preset parameter values
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

close(wb1);

guidata(hObject,handles);

% --- Outputs from this function are returned to the command line.
function varargout = MiA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listboxResults.
function listboxResults_Callback(hObject, eventdata, handles)
% hObject    handle to listboxResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxResults contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxResults

guidata(hObject,handles);

if ~isfield(handles,'eventFrames') || isempty(handles.eventFrames)
    return
end

event = get(hObject, 'Value'); %[]

if isempty(event)
    return
end

framesEvents = vertcat(handles.eventData.frames);
selectedFrame = framesEvents(event);
ind = find(vertcat(framesEvents)==selectedFrame);
nrEventsInSelectedFrame = length(ind);

axes(handles.axesImages); imshow(handles.images(:,:,selectedFrame));

if handles.mitosisDetectionDone && handles.displayMitosisDetectionResults
    axes(handles.axesImages); hold on;
    for i=1:nrEventsInSelectedFrame
        viscircles(handles.eventData(ind(i)).centres,handles.eventData(ind(i)).radii,'DrawBackgroundCircle',0,'Color',[0.75 0 0]);
    end
    viscircles(handles.eventData(event).centres,handles.eventData(event).radii,'DrawBackgroundCircle',0);
    hold off;
end

if handles.cellTrackingDone && handles.displayCellTrackingResults
    axes(handles.axesImages); hold on;
    for i=1:nrEventsInSelectedFrame
        contours = fliplr(handles.frameData(selectedFrame).contours{i});
        plot(contours(:,1),contours(:,2),'m','LineWidth',2);
    end
    hold off;
end

handles.currentImage = selectedFrame;
set(handles.textData,'String', handles.fileName{selectedFrame});
set(handles.editImgNr, 'String', num2str(selectedFrame));
set(handles.sliderImages, 'Value', selectedFrame);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function listboxResults_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxDisplayMitosisDetectionResults.
function checkboxDisplayMitosisDetectionResults_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDisplayMitosisDetectionResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDisplayMitosisDetectionResults

guidata(hObject,handles);

handles.displayMitosisDetectionResults = get(hObject,'Value');

guidata(hObject,handles);


% --- Executes on button press in checkboxDisplayCellTrackingResults.
function checkboxDisplayCellTrackingResults_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDisplayCellTrackingResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDisplayCellTrackingResults

guidata(hObject,handles);

handles.displayCellTrackingResults = get(hObject,'Value');

guidata(hObject,handles);


% --- Executes on button press in pushbuttonExamineImg.
function pushbuttonExamineImg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExamineImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

hMiA = getappdata(0,'hMiA');

setappdata(hMiA, 'img',       handles.images(:,:,handles.currentImage));
setappdata(hMiA, 'imageName', handles.fileName{handles.currentImage});
setappdata(hMiA, 'imageSize', handles.imageSize);
setappdata(hMiA, 'pixelSize', handles.pixelSize);

examine

uiwait

guidata(hObject,handles);


% --- Executes on button press in pushbuttonEdit.
function pushbuttonEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

hMiA = getappdata(0,'hMiA');

setappdata(hMiA,'nameDataSet',handles.nameDataSet);
setappdata(hMiA,'images',handles.images);
setappdata(hMiA,'currentImage',handles.currentImage);
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

parameters

uiwait

hMiA = getappdata(0,'hMiA');

handles.nameDataSet = getappdata(hMiA,'nameDataSet');
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

guidata(hObject,handles);


% --- Executes on button press in pushbuttonFullAnalysis.
function pushbuttonFullAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFullAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

% 1) Mitosis Detection

% Run mitosis detection algorithm
handles.eventData = mitosisDetection(handles.images,handles.numberOfFrames,handles.parameters.mitosisThreshold,handles.parameters.radiusMin,handles.parameters.radiusMax,handles.parameters.sensitivity);

handles.mitosisDetectionDone = 1;

handles.nrEvents = size(handles.eventData,2);
eventFrames = vertcat(handles.eventData.frames);
nrEventFrames = size(unique(eventFrames),1);
handles.eventFrames = zeros(nrEventFrames,1);
count = 1;
for i=1:nrEventFrames
    handles.eventFrames(i,1) = eventFrames(count);
    nrEventsInFrame = sum(eventFrames==eventFrames(count));
    handles.eventFrames(i,2) = nrEventsInFrame;
    count = count + nrEventsInFrame;
end

% Display overview of mitosis detection results in listbox
set(handles.textResults, 'String', 'Event Nr. | Frame Nr.');
handles.displayResults = {};
for i=1:handles.nrEvents
    handles.displayResults{i} = ['Event ',num2str(i),' | Frame ',num2str(handles.eventData(i).frames)];
end
set(handles.listboxResults, 'String', handles.displayResults);
set(handles.listboxResults, 'Value', 1);

% Display total number of detected events
set(handles.textTotalNrDetectedEvents,'String',['Total Number of Detected Events: ',num2str(handles.nrEvents)]);

set(handles.checkboxDisplayMitosisDetectionResults, 'Value', 1);
handles.displayMitosisDetectionResults = 1;
checkboxDisplayMitosisDetectionResults_Callback(hObject, eventdata, handles)
listboxResults_Callback(hObject, eventdata, handles)

pause(3)

% 2) Cell Tracking

% Create waitbar
wb2 = waitbar(0, 'TRACKING: Please wait...');
set(findobj(wb2,'type','patch'),'edgecolor','k','facecolor','b');

% Calculation of finite differences kernels

% Forwards differences
Dxpk = [1; -1; 0]; % D^x_+ (corresponds to [0; -1; 1])
% x corresponds to second dimension in image
Dypk = Dxpk'; % D^y_+ (corresponds to [0 -1 1])
% y corresponds to first dimension in image

% Backwards differences
Dxmk = [0; 1; -1]; % D^x_- (corresponds to [-1; 1; 0])
%x corresponds to second dimension in image
Dymk = Dxmk'; % D^y_- (corresponds to [-1 1 0])
%y corresponds to first dimension in image

% Central differences
Dxck = [1; 0; -1]/2; % D^x_c kernel (corresponds to [-1/2; 0; 1/2])
%x corresponds to second dimension in image
Dyck = Dxck'; % D^y_c kernel (corresponds to [-1/2 0 1/2])
%y corresponds to first dimension in image

m = handles.imageSize(1);
n = handles.imageSize(2);

% Run forwards and backwards tracking functions
for i=1:handles.nrEvents
    
    disp(['Starting backwards tracking for event ',num2str(i),'/',num2str(handles.nrEvents),'...'])
    [handles.eventData(i).phiBack,handles.eventData(i).contoursBack,handles.eventData(i).beginFrame,handles.eventData(i).statisticsBack,err] = backwardsTracking(handles.images,m,n,handles.eventData(i).frames,handles.numberOfFrames,handles.eventData(i).boundaries,handles.parameters.mitosisThreshold,handles.parameters.lambda1,handles.parameters.lambda2,handles.parameters.mu,handles.parameters.nu,handles.parameters.g_adjust_low,handles.parameters.g_adjust_high,handles.parameters.omega,handles.parameters.timeStep,handles.parameters.maxIterations,handles.parameters.phiUpdate,handles.parameters.epsilonNormGradReg,handles.parameters.epsilonDeltaReg,Dxpk,Dypk,Dxmk,Dymk,Dxck,Dyck,handles.gpu);
    disp(['Backwards tracking for event ',num2str(i),'/',num2str(handles.nrEvents),' done.'])
    waitbar((i*2-1)/(handles.nrEvents*2));
    
    disp(['Starting forwards tracking for event ',num2str(i),'/',num2str(handles.nrEvents),'...'])
    [handles.eventData(i).phiFor,handles.eventData(i).phiForSave,handles.eventData(i).contoursFor,handles.eventData(i).contoursForSave,handles.eventData(i).endFrame,handles.eventData(i).outcome,handles.eventData(i).statisticsFor] = forwardsTracking(handles.images,m,n,handles.eventData(i).frames,handles.numberOfFrames,handles.eventData(i).phiBack,handles.eventData(i).radii,handles.parameters.mitosisThreshold,handles.parameters.lambda1,handles.parameters.lambda2,handles.parameters.mu,handles.parameters.nu,handles.parameters.g_adjust_low,handles.parameters.g_adjust_high,handles.parameters.omega,handles.parameters.timeStep,handles.parameters.maxIterations,handles.parameters.phiUpdate,handles.parameters.epsilonNormGradReg,handles.parameters.epsilonDeltaReg,Dxpk,Dypk,Dxmk,Dymk,Dxck,Dyck,handles.gpu,err);
    disp(['Forwards tracking for event ',num2str(i),'/',num2str(handles.nrEvents),' done.'])
    waitbar((i*2)/(handles.nrEvents*2));
    
    handles.eventData(i).length = handles.eventData(i).endFrame - handles.eventData(i).beginFrame + 1;
    
    nrBack = size(handles.eventData(i).statisticsBack,2);
    count = 0;
    for j=nrBack:-1:1
        count = count+1;
        handles.eventData(i).contours(count) = handles.eventData(i).contoursBack(j);
        handles.eventData(i).masks{count} = (handles.eventData(i).phiBack(:,:,j)<0);
        handles.eventData(i).statistics(count) = handles.eventData(i).statisticsBack(j);
    end
    nrFor = size(handles.eventData(i).statisticsFor,2);
    count = nrBack;
    for j=1:nrFor
        count = count+1;
        handles.eventData(i).contours(count) = handles.eventData(i).contoursForSave(j);
        handles.eventData(i).masks{count} = (handles.eventData(i).phiForSave(:,:,j)<0);
    end
    if nrFor>2
        count = nrBack;
        for j=1:nrFor-1
            count = count+1;
            handles.eventData(i).statistics(count) = handles.eventData(i).statisticsFor(j);
        end
    end
    
end

for i=handles.nrEvents:-1:1
    if handles.eventData(i).beginFrame==1 || strcmp(handles.eventData(i).outcome,'unknown')
        handles.eventData(i) = [];
    end
end
handles.nrEvents = size(handles.eventData,2);

handles.cellTrackingDone = 1;

handles.frameData = struct([]);
for i = 1:handles.numberOfFrames
    handles.frameData(i).contours = [];
    handles.frameData(i).masks = [];
end
for i = 1:handles.nrEvents
    beginFrame = handles.eventData(i).beginFrame;
    endFrame = handles.eventData(i).endFrame;
    count = 0;
    for j = beginFrame:endFrame
        count = count+1;
        existingEvents = size(handles.frameData(j).contours,2);
        handles.frameData(j).contours{existingEvents+1} = handles.eventData(i).contours{count};
        handles.frameData(j).masks{existingEvents+1} = handles.eventData(i).masks{count};
    end
end

% Calculate and display average mitosis length
handles.averageLength = round(sum(vertcat(handles.eventData.length))/handles.nrEvents);
handles.averageLengthTime = handles.averageLength * handles.timeResolution;
set(handles.textAverageMitosisLength, 'String', ['Average Length of Mitotic Phase: ',num2str(handles.averageLengthTime),' min']);

% Display overview of tracking results in listbox
set(handles.textResults, 'String', 'Event Nr. | Frames | Duration (frames) | Duration (minutes) | Cell Fate');
handles.displayResults = {};
for i=1:handles.nrEvents
    handles.displayResults{i} = ['Event ',num2str(i),' | Frames ',num2str(handles.eventData(i).beginFrame),' - ',num2str(handles.eventData(i).endFrame),' | ',num2str(handles.eventData(i).length),' | ',num2str(handles.eventData(i).length * handles.timeResolution),' min | ',handles.eventData(i).outcome];
end
set(handles.listboxResults, 'String', handles.displayResults);
set(handles.listboxResults, 'Value', 1);

close(wb2);

set(handles.textTotalNrDetectedEvents,'String',['Total Number of Detected Events: ',num2str(handles.nrEvents)]);
set(handles.checkboxDisplayCellTrackingResults, 'Value', 1);
handles.displayCellTrackingResults = 1;
checkboxDisplayCellTrackingResults_Callback(hObject, eventdata, handles)
listboxResults_Callback(hObject, eventdata, handles)

guidata(hObject,handles);


% --- Executes on button press in pushbuttonMitosisDetection.
function pushbuttonMitosisDetection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMitosisDetection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

% Run mitosis detection algorithm
handles.eventData = mitosisDetection(handles.images,handles.numberOfFrames,handles.parameters.mitosisThreshold,handles.parameters.radiusMin,handles.parameters.radiusMax,handles.parameters.sensitivity);

handles.mitosisDetectionDone = 1;

handles.nrEvents = size(handles.eventData,2);
eventFrames = vertcat(handles.eventData.frames);
nrEventFrames = size(unique(eventFrames),1);
handles.eventFrames = zeros(nrEventFrames,1);
count = 1;
for i=1:nrEventFrames
    handles.eventFrames(i,1) = eventFrames(count);
    nrEventsInFrame = sum(eventFrames==eventFrames(count));
    handles.eventFrames(i,2) = nrEventsInFrame;
    count = count + nrEventsInFrame;
end

% Display overview of mitosis detection results in listbox
set(handles.textResults, 'String', 'Event Nr. | Frame Nr.');
handles.displayResults = {};
for i=1:handles.nrEvents
    handles.displayResults{i} = ['Event ',num2str(i),' | Frame ',num2str(handles.eventData(i).frames)];
end
set(handles.listboxResults, 'String', handles.displayResults);
set(handles.listboxResults, 'Value', 1);

% Display total number of detected events
set(handles.textTotalNrDetectedEvents,'String',['Total Number of Detected Events: ',num2str(handles.nrEvents)]);

set(handles.checkboxDisplayMitosisDetectionResults, 'Value', 1);
handles.displayMitosisDetectionResults = 1;
checkboxDisplayMitosisDetectionResults_Callback(hObject, eventdata, handles)
listboxResults_Callback(hObject, eventdata, handles)

guidata(hObject,handles);


% --- Executes on button press in pushbuttonCellTracking.
function pushbuttonCellTracking_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCellTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

if ~isfield(handles,'eventFrames')
    msgbox([{'Tracking cannot be executed.'},{'Please run mitosis detection first!'}]);
    return
end

% Create waitbar
wb2 = waitbar(0, 'TRACKING: Please wait...');
set(findobj(wb2,'type','patch'),'edgecolor','k','facecolor','b');

% Calculation of finite differences kernels

% Forwards differences
Dxpk = [1; -1; 0]; % D^x_+ (corresponds to [0; -1; 1])
% x corresponds to second dimension in image
Dypk = Dxpk'; % D^y_+ (corresponds to [0 -1 1])
% y corresponds to first dimension in image

% Backwards differences
Dxmk = [0; 1; -1]; % D^x_- (corresponds to [-1; 1; 0])
%x corresponds to second dimension in image
Dymk = Dxmk'; % D^y_- (corresponds to [-1 1 0])
%y corresponds to first dimension in image

% Central differences
Dxck = [1; 0; -1]/2; % D^x_c kernel (corresponds to [-1/2; 0; 1/2])
%x corresponds to second dimension in image
Dyck = Dxck'; % D^y_c kernel (corresponds to [-1/2 0 1/2])
%y corresponds to first dimension in image

m = handles.imageSize(1);
n = handles.imageSize(2);

% Run forwards and backwards tracking functions
for i=1:handles.nrEvents
    
    disp(['Starting backwards tracking for event ',num2str(i),'/',num2str(handles.nrEvents),'...'])
    [handles.eventData(i).phiBack,handles.eventData(i).contoursBack,handles.eventData(i).beginFrame,handles.eventData(i).statisticsBack,err] = backwardsTracking(handles.images,m,n,handles.eventData(i).frames,handles.numberOfFrames,handles.eventData(i).boundaries,handles.parameters.mitosisThreshold,handles.parameters.lambda1,handles.parameters.lambda2,handles.parameters.mu,handles.parameters.nu,handles.parameters.g_adjust_low,handles.parameters.g_adjust_high,handles.parameters.omega,handles.parameters.timeStep,handles.parameters.maxIterations,handles.parameters.phiUpdate,handles.parameters.epsilonNormGradReg,handles.parameters.epsilonDeltaReg,Dxpk,Dypk,Dxmk,Dymk,Dxck,Dyck,handles.gpu);
    disp(['Backwards tracking for event ',num2str(i),'/',num2str(handles.nrEvents),' done.'])
    waitbar((i*2-1)/(handles.nrEvents*2));
    
    disp(['Starting forwards tracking for event ',num2str(i),'/',num2str(handles.nrEvents),'...'])
    [handles.eventData(i).phiFor,handles.eventData(i).phiForSave,handles.eventData(i).contoursFor,handles.eventData(i).contoursForSave,handles.eventData(i).endFrame,handles.eventData(i).outcome,handles.eventData(i).statisticsFor] = forwardsTracking(handles.images,m,n,handles.eventData(i).frames,handles.numberOfFrames,handles.eventData(i).phiBack,handles.eventData(i).radii,handles.parameters.mitosisThreshold,handles.parameters.lambda1,handles.parameters.lambda2,handles.parameters.mu,handles.parameters.nu,handles.parameters.g_adjust_low,handles.parameters.g_adjust_high,handles.parameters.omega,handles.parameters.timeStep,handles.parameters.maxIterations,handles.parameters.phiUpdate,handles.parameters.epsilonNormGradReg,handles.parameters.epsilonDeltaReg,Dxpk,Dypk,Dxmk,Dymk,Dxck,Dyck,handles.gpu,err);
    disp(['Forwards tracking for event ',num2str(i),'/',num2str(handles.nrEvents),' done.'])
    waitbar((i*2)/(handles.nrEvents*2));
    
    handles.eventData(i).length = handles.eventData(i).endFrame - handles.eventData(i).beginFrame + 1;

    nrBack = size(handles.eventData(i).statisticsBack,2);
    count = 0;
    for j=nrBack:-1:1
        count = count+1;
        handles.eventData(i).contours(count) = handles.eventData(i).contoursBack(j);
        handles.eventData(i).masks{count} = (handles.eventData(i).phiBack(:,:,j)<0);
        handles.eventData(i).statistics(count) = handles.eventData(i).statisticsBack(j);
    end
    nrFor = size(handles.eventData(i).statisticsFor,2);
    count = nrBack;
    for j=1:nrFor
        count = count+1;
        handles.eventData(i).contours(count) = handles.eventData(i).contoursForSave(j);
        handles.eventData(i).masks{count} = (handles.eventData(i).phiForSave(:,:,j)<0);
    end
    if nrFor>2
        count = nrBack;
        for j=1:nrFor-1
            count = count+1;
            handles.eventData(i).statistics(count) = handles.eventData(i).statisticsFor(j);
        end
    end
    
end

for i=handles.nrEvents:-1:1
    if handles.eventData(i).beginFrame==1 || strcmp(handles.eventData(i).outcome,'unknown')
        handles.eventData(i) = [];
    end
end
handles.nrEvents = size(handles.eventData,2);

handles.cellTrackingDone = 1;

handles.frameData = struct([]);
for i = 1:handles.numberOfFrames
    handles.frameData(i).contours = [];
    handles.frameData(i).masks = [];
end
for i = 1:handles.nrEvents
    beginFrame = handles.eventData(i).beginFrame;
    endFrame = handles.eventData(i).endFrame;
    count = 0;
    for j = beginFrame:endFrame
        count = count+1;
        existingEvents = size(handles.frameData(j).contours,2);
        handles.frameData(j).contours{existingEvents+1} = handles.eventData(i).contours{count};
        handles.frameData(j).masks{existingEvents+1} = handles.eventData(i).masks{count};
    end
end

% Calculate and display average mitosis length
handles.averageLength = round(sum(vertcat(handles.eventData.length))/handles.nrEvents);
handles.averageLengthTime = handles.averageLength * handles.timeResolution;
set(handles.textAverageMitosisLength, 'String', ['Average Length of Mitotic Phase: ',num2str(handles.averageLengthTime),' min']);

% Display overview of tracking results in listbox
set(handles.textResults, 'String', 'Event Nr. | Frames | Duration (frames) | Duration (minutes) | Cell Fate');
handles.displayResults = {};
for i=1:handles.nrEvents
    handles.displayResults{i} = ['Event ',num2str(i),' | Frames ',num2str(handles.eventData(i).beginFrame),' - ',num2str(handles.eventData(i).endFrame),' | ',num2str(handles.eventData(i).length),' | ',num2str(handles.eventData(i).length * handles.timeResolution),' min | ',handles.eventData(i).outcome];
end
set(handles.listboxResults, 'String', handles.displayResults);
set(handles.listboxResults, 'Value', 1);

close(wb2);

set(handles.textTotalNrDetectedEvents,'String',['Total Number of Detected Events: ',num2str(handles.nrEvents)]);
set(handles.checkboxDisplayCellTrackingResults, 'Value', 1);
handles.displayCellTrackingResults = 1;
checkboxDisplayCellTrackingResults_Callback(hObject, eventdata, handles)
listboxResults_Callback(hObject, eventdata, handles)

guidata(hObject,handles);


% --- Executes on button press in pushbuttonDeleteEvent.
function pushbuttonDeleteEvent_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDeleteEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

if handles.mitosisDetectionDone==0
    errordlg('Nothing to delete. Please start Mitosis Detection first.');
    return
end
if handles.nrEvents==0
    errordlg('There are no events to delete.');
    return
end

deleteEvent = get(handles.listboxResults, 'Value');

if isempty(deleteEvent)
   errordlg('All events have been deleted. No further analysis possible.');
   return
end

handles.nrEvents = handles.nrEvents - 1;
set(handles.textTotalNrDetectedEvents,'String',['Total Number of Detected Events: ',num2str(handles.nrEvents)]);
if handles.nrEvents==0
    msgbox('All events have been deleted. No further analysis possible.');
    handles.displayResults = {};
    set(handles.listboxResults, 'String', handles.displayResults);
    axes(handles.axesImages); imshow(handles.images(:,:,handles.eventFrames(deleteEvent,1)));
    return
end
deleteFrame = handles.eventData(deleteEvent).frames;
ind = find(handles.eventFrames(:,1)==deleteFrame);
if handles.eventFrames(ind,2)==1
    handles.eventFrames(ind,:) = [];
else
    handles.eventFrames(ind,2) = handles.eventFrames(ind,2) - 1;
end

handles.eventData(deleteEvent) = [];

handles.displayResults = {};
if handles.cellTrackingDone
    for i=1:handles.nrEvents
        handles.displayResults{i} = ['Event ',num2str(i),' | Frames ',num2str(handles.eventData(i).beginFrame),' - ',num2str(handles.eventData(i).endFrame),' | ',num2str(handles.eventData(i).length),' | ',num2str(handles.eventData(i).length * handles.timeResolution),' min | ',handles.eventData(i).outcome];
    end
    handles.averageLength = round(sum(vertcat(handles.eventData.length))/handles.nrEvents);
    handles.averageLengthTime = handles.averageLength * handles.timeResolution;
    set(handles.textAverageMitosisLength, 'String', ['Average Length of Mitotic Phase: ',num2str(handles.averageLengthTime),' min']);
else
    for i=1:handles.nrEvents
        handles.displayResults{i} = ['Event ',num2str(i),' | Frame ',num2str(handles.eventData(i).frames)];
    end
end
set(handles.listboxResults, 'String', handles.displayResults);
set(handles.listboxResults, 'Value', 1);

guidata(hObject,handles);


% --- Executes on button press in pushbuttonAddBwd.
function pushbuttonAddBwd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddBwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

if (~isfield(handles,'eventFrames') || isempty(handles.eventFrames)) || handles.cellTrackingDone == 0;
    errordlg('No events to analyse / Cell Tracking not yet started.');
    return
end

addBwdEvent = get(handles.listboxResults, 'Value');

if handles.eventData(addBwdEvent).beginFrame == 1
    errordlg('Cannot add backwards frame since beginning frame is already first frame of whole image sequence!');
    return
else
    handles.eventData(addBwdEvent).beginFrame = handles.eventData(addBwdEvent).beginFrame - 1;
    nrBwdFrames = size(handles.eventData(addBwdEvent).phiBack,3) + 1;
    tempPhiBack = handles.eventData(addBwdEvent).phiBack;
    handles.eventData(addBwdEvent).phiBack = zeros(handles.imageSize(1),handles.imageSize(2),nrBwdFrames);
    handles.eventData(addBwdEvent).phiBack(:,:,2:nrBwdFrames) = tempPhiBack;
    tempContoursBack = handles.eventData(addBwdEvent).contoursBack;
    handles.eventData(addBwdEvent).contoursBack = cell(1,nrBwdFrames);
    for i = 2:nrBwdFrames
        handles.eventData(addBwdEvent).contoursBack{i} = tempContoursBack{i-1};
    end
    handles.eventData(addBwdEvent).length = handles.eventData(addBwdEvent).length + 1;
    handles.averageLength = round(sum(vertcat(handles.eventData.length))/handles.nrEvents);
    handles.averageLengthTime = handles.averageLength * handles.timeResolution;
    set(handles.textAverageMitosisLength, 'String', ['Average Length of Mitotic Phase: ',num2str(handles.averageLengthTime),' min']);
    handles.displayIndividualEvents{addBwdEvent} = ['Event ',num2str(addBwdEvent),': ',handles.eventData(addBwdEvent).outcome,' , ',num2str(handles.eventData(addBwdEvent).length * handles.timeResolution),' min'];
end    

guidata(hObject,handles);


% --- Executes on button press in pushbuttonAddFwd.
function pushbuttonAddFwd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddFwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

if (~isfield(handles,'eventFrames') || isempty(handles.eventFrames)) || handles.cellTrackingDone == 0;
    errordlg('No events to analyse / Cell Tracking not yet started.');
    return
end

addFwdEvent = get(handles.listboxResults, 'Value');

if handles.eventData(addFwdEvent).endFrame == handles.numberOfFrames
    errordlg('Cannot add forwards frame since end frame is already last frame of whole image sequence!');
    return
else
    handles.eventData(addFwdEvent).endFrame = handles.eventData(addFwdEvent).endFrame + 1;
    nrFwdFrames = size(handles.eventData(addFwdEvent).phiFor,3) + 1;
    handles.eventData(addFwdEvent).phiFor(:,:,nrFwdFrames) = zeros(handles.imageSize(1),handles.imageSize(2));
    handles.eventData(addFwdEvent).contoursFor{nrFwdFrames} = [];
    handles.eventData(addFwdEvent).length = handles.eventData(addFwdEvent).length + 1;
    handles.averageLength = round(sum(vertcat(handles.eventData.length))/handles.nrEvents);
    handles.averageLengthTime = handles.averageLength * handles.timeResolution;
    set(handles.textAverageMitosisLength, 'String', ['Average Length of Mitotic Phase: ',num2str(handles.averageLengthTime),' min']);
    handles.displayIndividualEvents{addFwdEvent} = ['Event ',num2str(addFwdEvent),': ',handles.eventData(addFwdEvent).outcome,' , ',num2str(handles.eventData(addFwdEvent).length * handles.timeResolution),' min'];
end    

guidata(hObject,handles);


% --- Executes on button press in pushbuttonDeleteBwd.
function pushbuttonDeleteBwd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDeleteBwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

if (~isfield(handles,'eventFrames') || isempty(handles.eventFrames)) || handles.cellTrackingDone == 0;
    errordlg('No events to analyse / Cell Tracking not yet started.');
    return
end

deleteBwdEvent = get(handles.listboxResults, 'Value');

if size(handles.eventData(deleteBwdEvent).phiBack,3) == 1
    errordlg('Cannot delete backwards frame since minimum number of 1 has already been reached!');
    return
else
    handles.eventData(deleteBwdEvent).beginFrame = handles.eventData(deleteBwdEvent).beginFrame + 1;
    handles.eventData(deleteBwdEvent).phiBack(:,:,1) = [];
    handles.eventData(deleteBwdEvent).contoursBack(1) = [];
    handles.eventData(deleteBwdEvent).length = handles.eventData(deleteBwdEvent).length - 1;
    handles.averageLength = round(sum(vertcat(handles.eventData.length))/handles.nrEvents);
    handles.averageLengthTime = handles.averageLength * handles.timeResolution;
    set(handles.textAverageMitosisLength, 'String', ['Average Length of Mitotic Phase: ',num2str(handles.averageLengthTime),' min']);
    handles.displayIndividualEvents{deleteBwdEvent} = ['Event ',num2str(deleteBwdEvent),': ',handles.eventData(deleteBwdEvent).outcome,' , ',num2str(handles.eventData(deleteBwdEvent).length * handles.timeResolution),' min'];
end    

guidata(hObject,handles);


% --- Executes on button press in pushbuttonDeleteFwd.
function pushbuttonDeleteFwd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDeleteFwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

if (~isfield(handles,'eventFrames') || isempty(handles.eventFrames)) || handles.cellTrackingDone == 0;
    errordlg('No events to analyse / Cell Tracking not yet started.');
    return
end

deleteFwdEvent = get(handles.listboxResults, 'Value');

size(handles.eventData(deleteFwdEvent).phiFor,3)

if size(handles.eventData(deleteFwdEvent).phiFor,3) == 1
    errordlg('Cannot delete forwards frame since minimum number of 1 has already been reached!');
    return
else
    nrFwdFrames = size(handles.eventData(deleteFwdEvent).phiFor,3) - 1;
    handles.eventData(deleteFwdEvent).endFrame = handles.eventData(deleteFwdEvent).endFrame - 1;
    handles.eventData(deleteFwdEvent).phiFor(:,:,nrFwdFrames) = [];
    handles.eventData(deleteFwdEvent).contoursFor(nrFwdFrames) = [];
    handles.eventData(deleteFwdEvent).length = handles.eventData(deleteFwdEvent).length - 1;
    handles.averageLength = round(sum(vertcat(handles.eventData.length))/handles.nrEvents);
    handles.averageLengthTime = handles.averageLength * handles.timeResolution;
    set(handles.textAverageMitosisLength, 'String', ['Average Length of Mitotic Phase: ',num2str(handles.averageLengthTime),' min']);
    handles.displayIndividualEvents{deleteFwdEvent} = ['Event ',num2str(deleteFwdEvent),': ',handles.eventData(deleteFwdEvent).outcome,' , ',num2str(handles.eventData(deleteFwdEvent).length * handles.timeResolution),' min'];
end  

guidata(hObject,handles);


% --- Executes on button press in pushbuttonStatistics.
function pushbuttonStatistics_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStatistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

if handles.mitosisDetectionDone == 0 || handles.cellTrackingDone==0
    errordlg('No statistics available.');
    return
end

hMiA = getappdata(0,'hMiA');

setappdata(hMiA, 'numberOfFrames',    handles.numberOfFrames);
setappdata(hMiA, 'imageSize',         handles.imageSize);
setappdata(hMiA, 'images',            handles.images);
setappdata(hMiA, 'nameDataSet',       handles.nameDataSet);
setappdata(hMiA, 'pixelSize',         handles.pixelSize);
setappdata(hMiA, 'timeResolution',    handles.timeResolution);
setappdata(hMiA, 'eventData',         handles.eventData);
setappdata(hMiA, 'nrEvents',          handles.nrEvents);
setappdata(hMiA, 'averageLength',     handles.averageLength);
setappdata(hMiA, 'averageLengthTime', handles.averageLengthTime);

statistics

uiwait

guidata(hObject,handles);


% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

if handles.mitosisDetectionDone == 0
    errordlg({'No data available.','Please start Mitosis Detection first.'});
    return
end

parameters = handles.parameters;
frameData = handles.frameData;
eventData = handles.eventData;
eventsData = rmfield(eventData,{'centres','radii','metrics',...
    'boundaries','phiBack','contoursBack','statisticsBack','phiFor',...
    'phiForSave','contoursFor','contoursForSave','statisticsFor',...
    'masks','statistics'});
try
    save results/resultData.mat parameters eventsData frameData
catch
    save results/resultData.mat parameters eventsData
end

t1 = struct2table(parameters);
writetable(t1,'results/tableParameters.csv')
fates = cell(handles.nrEvents,1);
for i=1:handles.nrEvents
    fates{i} = eventData(i).outcome;
end
t2 = table((1:handles.nrEvents)',vertcat(eventData.beginFrame),...
    vertcat(eventData.endFrame),vertcat(eventData.length),...
    vertcat(eventData.length)*handles.timeResolution,fates,...
    'VariableNames',{'Event','Start_frame','End_frame',...
    'Duration_frames','Duration_minutes','Fate'});
writetable(t2,'results/tableData.csv')

guidata(hObject,handles);


% --- Executes on slider movement.
function sliderImages_Callback(hObject, eventdata, handles)
% hObject    handle to sliderImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

guidata(hObject,handles);

value = round(get(hObject,'Value'));

set(handles.textData, 'String', handles.fileName{value});
set(handles.editImgNr, 'String', num2str(value));
handles.currentImage = value;

axes(handles.axesImages), imshow(handles.images(:,:,value));

guidata(hObject,handles);

if handles.mitosisDetectionDone
    ind = find(vertcat(handles.eventData.frames)==value);
else
    return
end

if ~isempty(ind)
    nrEvents = length(ind);
    if handles.mitosisDetectionDone && handles.displayMitosisDetectionResults
        axes(handles.axesImages); hold on;
        for i=1:nrEvents
            viscircles(handles.eventData(ind(i)).centres,handles.eventData(ind(i)).radii,'DrawBackgroundCircle',0);
        end
        hold off;
    end
end

if handles.cellTrackingDone
    if ~isempty(handles.frameData(value).contours)
        nrEvents = size(handles.frameData(value).contours,2);
        if handles.displayCellTrackingResults
            axes(handles.axesImages); hold on;
            for i=1:nrEvents
                contours = fliplr(handles.frameData(value).contours{i});
                plot(contours(:,1),contours(:,2),'m','LineWidth',2);
            end
            hold off;
        end
    end
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sliderImages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editImgNr_Callback(hObject, eventdata, handles)
% hObject    handle to editImgNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editImgNr as text
%        str2double(get(hObject,'String')) returns contents of editImgNr as a double

guidata(hObject,handles);

imgNr = str2double(get(hObject,'String'));
set(handles.sliderImages, 'Value', imgNr);

set(handles.textData, 'String', handles.fileName{imgNr});
handles.currentImage = imgNr;

axes(handles.axesImages), imshow(handles.images(:,:,imgNr));

guidata(hObject,handles);

if handles.mitosisDetectionDone
    ind = find(vertcat(handles.eventData.frames)==imgNr);
else
    return
end

if ~isempty(ind)
    nrEvents = length(ind);
    if handles.mitosisDetectionDone && handles.displayMitosisDetectionResults
        axes(handles.axesImages); hold on;
        for i=1:nrEvents
            viscircles(handles.eventData(ind(i)).centres,handles.eventData(ind(i)).radii,'DrawBackgroundCircle',0);
        end
        hold off;
    end
end

if handles.cellTrackingDone
    if ~isempty(handles.frameData(imgNr).contours)
        nrEvents = size(handles.frameData(imgNr).contours,2);
        if handles.displayCellTrackingResults
            axes(handles.axesImages); hold on;
            for i=1:nrEvents
                contours = fliplr(handles.frameData(imgNr).contours{i});
                plot(contours(:,1),contours(:,2),'m','LineWidth',2);
            end
            hold off;
        end
    end
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editImgNr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editImgNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
