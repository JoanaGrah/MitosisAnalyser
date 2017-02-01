function varargout = statistics(varargin)
% STATISTICS MATLAB code for statistics.fig
%      STATISTICS, by itself, creates a new STATISTICS or raises the existing
%      singleton*.
%
%      H = STATISTICS returns the handle to a new STATISTICS or the handle to
%      the existing singleton*.
%
%      STATISTICS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STATISTICS.M with the given input arguments.
%
%      STATISTICS('Property','Value',...) creates a new STATISTICS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before statistics_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to statistics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help statistics

% Last Modified by GUIDE v2.5 18-Aug-2015 12:38:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @statistics_OpeningFcn, ...
                   'gui_OutputFcn',  @statistics_OutputFcn, ...
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


% --- Executes just before statistics is made visible.
function statistics_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to statistics (see VARARGIN)

% Choose default command line output for statistics
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes statistics wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Get parameters from main MiA GUI
hMiA = getappdata(0,'hMiA');

handles.numberOfFrames = getappdata(hMiA, 'numberOfFrames');
handles.imageSize = getappdata(hMiA, 'imageSize');
handles.images = getappdata(hMiA, 'images');
handles.nameDataSet = getappdata(hMiA, 'nameDataSet');
handles.pixelSize = getappdata(hMiA, 'pixelSize');
handles.timeResolution = getappdata(hMiA, 'timeResolution');
handles.eventData = getappdata(hMiA, 'eventData');
handles.nrEvents = getappdata(hMiA, 'nrEvents');
handles.averageLength = getappdata(hMiA, 'averageLength');
handles.averageLengthTime = getappdata(hMiA, 'averageLengthTime');

set(handles.textDataSetDescriptionResult,'String',handles.nameDataSet);
set(handles.textTotalNumberOfImagesResult,'String',num2str(handles.numberOfFrames));
m = handles.imageSize(1);
n = handles.imageSize(2);
stringImageSizeResult = [num2str(m),' x ',num2str(n),' pixels / ',num2str(m*handles.pixelSize),' x ',num2str(n*handles.pixelSize),' microns'];
set(handles.textImageSizeResult,'String',stringImageSizeResult);
if handles.pixelSize==1
    stringPixelSizeResult = [num2str(handles.pixelSize),' micron'];
else
    stringPixelSizeResult = [num2str(handles.pixelSize),' microns'];
end
set(handles.textPixelSizeResult,'String',stringPixelSizeResult);
stringTimeIntervalResult = [num2str(handles.timeResolution),' minutes'];
set(handles.textTimeIntervalResult,'String',stringTimeIntervalResult);
set(handles.textTotalNumberOfEventsResult,'String',handles.nrEvents);
stringAverageMitosisDurationResult = [num2str(handles.averageLength),' frames / ',num2str(handles.averageLengthTime),' minutes'];
set(handles.textAverageMitosisDurationResult,'String',stringAverageMitosisDurationResult);
stringChooseEvent = cell(handles.nrEvents,1);
for i=1:handles.nrEvents
    stringChooseEvent{i} = ['Event ',num2str(i)];
end
set(handles.popupmenuChooseEvent,'String',stringChooseEvent);
set(handles.popupmenuChooseEvent,'Value',1);

eventNr = 1;
eventLength = handles.eventData(eventNr).length;

stringEventDuration = ['Duration: ',num2str(eventLength), ' frms / ',num2str(eventLength*handles.timeResolution),' mins'];
set(handles.textEventDuration,'String',stringEventDuration);
stringEventFate = ['Fate: ',handles.eventData(eventNr).outcome];
set(handles.textEventFate,'String',stringEventFate);

areas = cell(eventLength-1,1);
for i=1:eventLength-1
    areas{i} = handles.eventData(eventNr).statistics(i).area;
end
areas = cell2mat(areas);
perimeters = cell(eventLength-1,1);
for i=1:eventLength-1
    perimeters{i} = handles.eventData(eventNr).statistics(i).perimeter;
end
perimeters = cell2mat(perimeters);
circularities = cell(eventLength-1,1);
for i=1:eventLength-1
    circularities{i} = handles.eventData(eventNr).statistics(i).circularity;
end
circularities = cell2mat(circularities);
centroids = cell(eventLength-1,1);
for i=1:eventLength-1
    centroids{i} = handles.eventData(eventNr).statistics(i).centroid;
end
centroids = cell2mat(centroids);
axes(handles.axesArea), plot(1:eventLength-1,areas,'--.','MarkerSize',25); title('Area');
axes(handles.axesPerimeter), plot(1:eventLength-1,perimeters,'--.','Color',[0.8500 0.3250 0.0980],'MarkerSize',25); title('Perimeter');
axes(handles.axesCircularity), plot(1:eventLength-1,circularities,'--.','Color',[0.9290 0.6940 0.1250],'MarkerSize',25); title('Circularity');
axes(handles.axesCentroid), scatter(centroids(:,1),centroids(:,2),100,'filled','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560]);...
    hold on;plot(centroids(:,1),centroids(:,2),'Color',[0.4940 0.1840 0.5560]); title('Centroid Movement');hold off;

outcomes = cell(handles.nrEvents,1);
for i=1:handles.nrEvents
    outcomes{i} = handles.eventData(i).outcome;
end
pieData(1) = sum(strcmp(outcomes,'no division'));
pieData(2) = sum(strcmp(outcomes,'regular'));
pieData(3) = sum(strcmp(outcomes,'3 daughter cells'));
pieData(4) = sum(strcmp(outcomes,'apoptosis'));
pieLabels = {'1DC','2DC','3DC','CD'};
warning('off','MATLAB:pie:NonPositiveData');
axes(handles.axesPieChart), pie(pieData,pieLabels);
warning('on','MATLAB:pie:NonPositiveData');

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = statistics_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuChooseEvent.
function popupmenuChooseEvent_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuChooseEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuChooseEvent contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuChooseEvent

guidata(hObject, handles);

eventNr = get(hObject,'Value');
eventLength = handles.eventData(eventNr).length;

stringEventDuration = ['Duration: ',num2str(eventLength), ' frms / ',num2str(eventLength*handles.timeResolution),' mins'];
set(handles.textEventDuration,'String',stringEventDuration);
stringEventFate = ['Fate: ',handles.eventData(eventNr).outcome];
set(handles.textEventFate,'String',stringEventFate);

areas = cell(eventLength-1,1);
for i=1:eventLength-1
    areas{i} = handles.eventData(eventNr).statistics(i).area;
end
areas = cell2mat(areas);
perimeters = cell(eventLength-1,1);
for i=1:eventLength-1
    perimeters{i} = handles.eventData(eventNr).statistics(i).perimeter;
end
perimeters = cell2mat(perimeters);
circularities = cell(eventLength-1,1);
for i=1:eventLength-1
    circularities{i} = handles.eventData(eventNr).statistics(i).circularity;
end
circularities = cell2mat(circularities);
centroids = cell(eventLength-1,1);
for i=1:eventLength-1
    centroids{i} = handles.eventData(eventNr).statistics(i).centroid;
end
centroids = cell2mat(centroids);
axes(handles.axesArea), plot(1:eventLength-1,areas,'--.','MarkerSize',25); title('Area');
axes(handles.axesPerimeter), plot(1:eventLength-1,perimeters,'--.','Color',[0.8500 0.3250 0.0980],'MarkerSize',25); title('Perimeter');
axes(handles.axesCircularity), plot(1:eventLength-1,circularities,'--.','Color',[0.9290 0.6940 0.1250],'MarkerSize',25); title('Circularity');
axes(handles.axesCentroid), scatter(centroids(:,1),centroids(:,2),100,'filled','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560]);...
    hold on;plot(centroids(:,1),centroids(:,2),'Color',[0.4940 0.1840 0.5560]); title('Centroid Movement');hold off;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenuChooseEvent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuChooseEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonShowTrackingResultsContours.
function pushbuttonShowTrackingResultsContours_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShowTrackingResultsContours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

eventNr = get(handles.popupmenuChooseEvent,'Value');
eventLength = handles.eventData(eventNr).length;

nrRowsColumns = ceil(sqrt(eventLength));
nrSubplots = nrRowsColumns^2;

frameNr = handles.eventData(eventNr).beginFrame-1;

figure;
for i=1:nrSubplots
    if i<=eventLength
        frameNr = frameNr + 1;
        contours = fliplr(handles.eventData(eventNr).contours{i});
        subplot(nrRowsColumns,nrRowsColumns,i);imshow(handles.images(:,:,frameNr));hold on;...
            plot(contours(:,1),contours(:,2),'m','LineWidth',1);title(['Frame nr ',num2str(frameNr)]);hold off;
    else
        subplot(nrRowsColumns,nrRowsColumns,i,'replace');axis off;axis image;
    end
end     

guidata(hObject, handles);


% --- Executes on button press in pushbuttonShowTrackingResultsMasks.
function pushbuttonShowTrackingResultsMasks_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShowTrackingResultsMasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

eventNr = get(handles.popupmenuChooseEvent,'Value');
eventLength = handles.eventData(eventNr).length;

nrRowsColumns = ceil(sqrt(eventLength));
nrSubplots = nrRowsColumns^2;

frameNr = handles.eventData(eventNr).beginFrame-1;

figure;
for i=1:nrSubplots
    if i<=eventLength
        frameNr = frameNr + 1;
        subplot(nrRowsColumns,nrRowsColumns,i);imshow(handles.eventData(eventNr).masks{i});title(['Frame nr ',num2str(frameNr)]);
    else
        subplot(nrRowsColumns,nrRowsColumns,i,'replace');axis off;axis image;
    end
end     

guidata(hObject, handles);
