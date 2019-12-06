function varargout = GUI_example(varargin)
% GUI_EXAMPLE MATLAB code for GUI_example.fig
%      GUI_EXAMPLE, by itself, creates a new GUI_EXAMPLE or raises the existing
%      singleton*.
%
%      H = GUI_EXAMPLE returns the handle to a new GUI_EXAMPLE or the handle to
%      the existing singleton*.
%
%      GUI_EXAMPLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EXAMPLE.M with the given input arguments.
%
%      GUI_EXAMPLE('Property','Value',...) creates a new GUI_EXAMPLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_example_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_example_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_example

% Last Modified by GUIDE v2.5 01-Dec-2019 19:34:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_example_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_example_OutputFcn, ...
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


% --- Executes just before GUI_example is made visible.
function GUI_example_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_example (see VARARGIN)

% Choose default command line output for GUI_example
handles.output = hObject;
handles.sep = '\';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_example wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_example_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in preprocessStatusSetButton.
function preprocessStatusSetButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessStatusSetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
targetVariable = get(handles.preprocessStatusMenu,'String');
targetVariable = targetVariable{get(handles.preprocessStatusMenu,'Value')};
if ~strcmp(targetVariable, 'None')
    configTable = get(handles.uitable1,'Data');
    configTableColumnName = get(handles.uitable1,'ColumnName');
    configTable = cell2table(configTable);
    configTable.Properties.VariableNames = configTableColumnName;

    configTable.(targetVariable) = ones(size(configTable,1),1);
    set(handles.uitable1,'Data', table2cell(configTable));
end

% --- Executes on button press in loadAnimalConfigButton.
function loadAnimalConfigButton_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
% hObject    handle to loadAnimalConfigButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mouse = get(handles.editAnimal,'String');
masterInfo = func_loadMasterConfig('Mouse',mouse,'Root',handles.configRoot);
masterInfo = masterInfo.(mouse);
%[csvname, csvpath] = uigetfile('*.csv');
handles.sessionConfigName = [handles.configRoot handles.sep 'mouse' handles.sep masterInfo.sessionConfig];
handles.imagingConfigName = [handles.configRoot handles.sep 'imaging' handles.sep masterInfo.imagingConfig];
handles.behaviorConfigName = [handles.configRoot handles.sep 'behavior' handles.sep masterInfo.behaviorConfig];
configTable = readtable(handles.sessionConfigName);
set(handles.uitable1,'Data',table2cell(configTable));
set(handles.uitable1,'ColumnName', configTable.Properties.VariableNames);
handles.configTable = configTable;
%handles.csvname = csvname;
%handles.csvpath = csvpath;
handles.mouse = mouse; 
handles.masterInfo = masterInfo;
guidata(hObject,handles);

% --- Executes on selection change in preprocessStatusMenu.
function preprocessStatusMenu_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to preprocessStatusMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns preprocessStatusMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from preprocessStatusMenu


% --- Executes during object creation, after setting all properties.
function preprocessStatusMenu_CreateFcn(hObject, ~, handles)
% hObject    handle to preprocessStatusMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.listboxSelectedIndex = eventdata.Source.Value;
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in preprocessMenu.
function preprocessMenu_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns preprocessMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from preprocessMenu
allOpt = get(hObject,'String');
selectedOpt = allOpt{eventdata.Source.Value};
%configTable = get(handles.uitable1,'Data');
%configTableColumnName = get(handles.uitable1,'ColumnName');
%onfigTable = cell2table(configTable);
%configTable.Properties.VariableNames = configTableColumnName;
configTable = getConfigTable(handles);

switch selectedOpt
    case 'TC extraction'
        filenames = configTable.('ImagingFile');
    case 'Deconvolve'
        filenames = configTable.('ImagingFile');
        
    case 'ROI Tracking'
        filenames = configTable.('ImagingFile');
    case 'Tuning Curve'
        filenames = 'SelectFile';
    case 'None'
        filenames = 'SelectFile';
end

set(handles.listbox1,'String',filenames);

% --- Executes during object creation, after setting all properties.
function preprocessMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to preprocessMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function configTable = getConfigTable(handles)
configTable = get(handles.uitable1,'Data');
configTableColumnName = get(handles.uitable1,'ColumnName');
configTable = cell2table(configTable);
configTable.Properties.VariableNames = configTableColumnName;
    


% --- Executes on button press in preprocessPreviewButton.
function preprocessPreviewButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessPreviewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imagingConfigDisplay = func_loadImagingConfig(handles.imagingConfigName,'Split', false);
imagingConfig = func_loadImagingConfig(handles.imagingConfigName,'Split', true);
paramValueDisplay = {};
paramValue = {};
param = [fieldnames(handles.masterInfo); fieldnames(imagingConfig)];
temp = fieldnames(handles.masterInfo);
for i = 1:length(temp)
    paramValue = [paramValue;temp{i};handles.masterInfo.(temp{i})];
    paramValueDisplay = [paramValueDisplay; handles.masterInfo.(temp{i})];
end
temp = fieldnames(imagingConfig);
tempDisplay = fieldnames(imagingConfigDisplay);
for i = 1:length(temp)
    paramValue = [paramValue;temp{i};{imagingConfig.(temp{i})}];
    paramValueDisplay = [paramValueDisplay; {imagingConfigDisplay.(tempDisplay{i})}];
end


filenames = get(handles.listbox1,'String');
if ~iscell(filenames)
    filenames = {filenames};
end
% add sbxpath h5path and suite2ppath into the parameter
filename = filenames{1};
configTable = getConfigTable(handles);
sbxpath = configTable.sbxpath(strcmp(configTable.ImagingFile,filename));
h5path = configTable.h5path(strcmp(configTable.ImagingFile,filename));
suite2ppath = configTable.suite2ppath(strcmp(configTable.ImagingFile,filename));
param = [param;'sbxpath';'h5path';'suite2ppath';'root';'savepath'];
paramValue = [paramValue;'sbxpath';sbxpath;'h5path';h5path;'suite2ppath';suite2ppath;...
    'root';handles.configRoot; 'savepath'; handles.savePath];
paramValueDisplay = [paramValueDisplay;sbxpath;h5path;suite2ppath;handles.configRoot;handles.savePath];
% add filename into the parameter
paramValue = [paramValue;'filename';filenames];
% add nFrames_onePlane
nFrames_oneplaneStr = configTable.nFrames_onePlane(strcmp(configTable.ImagingFile,filename));
if ~iscell(nFrames_oneplaneStr); nFrames_oneplaneStr = {nFrames_oneplaneStr}; end
for i = 1:length(nFrames_oneplaneStr)
    nFrames_oneplaneSplit = strsplit(nFrames_oneplaneStr{i});
    for j = 1:str2double(imagingConfig.nPlanes)
        nFrames_oneplane(i,j) = str2double(nFrames_oneplaneSplit{j});
    end
end
for j = 1:str2double(imagingConfig.nPlanes)
    nFrames_oneplane(:,j) = cumsum(nFrames_oneplane(:,j));
end
nFrames_oneplane = [zeros(1,str2double(imagingConfig.nPlanes));nFrames_oneplane];
paramValue = [paramValue;'nFrames_oneplane';nFrames_oneplane];
% add roi, TC file names
paramValue = [paramValue;'datapath';configTable.dataPath(strcmp(configTable.ImagingFile,filename))];
paramValue = [paramValue;'roiFile';configTable.roiFile(strcmp(configTable.ImagingFile,filename))];
paramValue = [paramValue;'tcFile';configTable.tcFile(strcmp(configTable.ImagingFile,filename))];

% display parameters and save it to the handle
paramTableDisplay = table(param,paramValueDisplay,'VariableNames',{'name','value'});
set(handles.uitable2,'Data',table2cell(paramTableDisplay));
set(handles.uitable2,'ColumnName', paramTableDisplay.Properties.VariableNames);
handles.param = paramValue;
guidata(hObject,handles);

% --- Executes on button press in fileAddButton.
function fileAddButton_Callback(hObject, eventdata, handles)
% hObject    handle to fileAddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listBoxStr = get(handles.listbox1,'String');
if ischar(listBoxStr)
    listBoxStr = {listBoxStr};
end
configTable = get(handles.uitable1,'Data');
for i = 1:size(handles.uitableSelectedIndex,1)
    tempFilename = configTable{handles.uitableSelectedIndex(i,1), handles.uitableSelectedIndex(i,2)};
    if ~any(strcmp(listBoxStr, tempFilename))
        listBoxStr = [listBoxStr; tempFilename];
    end
end
%if strcmp(listBoxStr,'SelectFile') || strcmp(listBoxStr{1},'SelectFile')
listBoxStr(strcmp(listBoxStr,'SelectFile')) = [];
%end
set(handles.listbox1,'String',listBoxStr);

% --- Executes on button press in fileDropButton.
function fileDropButton_Callback(hObject, eventdata, handles)
% hObject    handle to fileDropButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listBoxStr = get(handles.listbox1,'String');
if ischar(listBoxStr)
   listBoxStr = {listBoxStr}; 
end
selectBoxStr = listBoxStr(get(handles.listbox1,'Value'));

if length(listBoxStr) >2
    listBoxStr(strcmp(selectBoxStr,listBoxStr)) = [];
elseif length(listBoxStr)==2
    listBoxStr(strcmp(selectBoxStr,listBoxStr)) = [];
    if length(listBoxStr)==1
        listBoxStr = listBoxStr{1};
    end
elseif strcmp(listBoxStr,selectBoxStr)
    listBoxStr = 'SelectFile';
else
    disp('CHECK THIS!!!')
end

set(handles.listbox1,'String',listBoxStr);

% --- Executes on button press in preprocessStatusClearButton.
function preprocessStatusClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessStatusClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
targetVariable = get(handles.preprocessStatusMenu,'String');
targetVariable = targetVariable{get(handles.preprocessStatusMenu,'Value')};
if ~strcmp(targetVariable, 'None')
    configTable = get(handles.uitable1,'Data');
    configTableColumnName = get(handles.uitable1,'ColumnName');
    configTable = cell2table(configTable);
    configTable.Properties.VariableNames = configTableColumnName;

    configTable.(targetVariable) = zeros(size(configTable,1),1);
    set(handles.uitable1,'Data', table2cell(configTable));
end

% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
handles.uitableSelectedIndex = eventdata.Indices;
guidata(hObject,handles);



% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in preprocessStatusEditButton.
function preprocessStatusEditButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessStatusEditButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value==1
    set(handles.uitable1,'ColumnEditable',true(1,length(get(handles.uitable1,'ColumnName'))));
else
    set(handles.uitable1,'ColumnEditable',false(1,length(get(handles.uitable1,'ColumnName'))));
end
guidata(hObject,handles);

% --- Executes on button press in preprocessStatusSaveButton.
function preprocessStatusSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessStatusSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
configTable = get(handles.uitable1,'Data');
configTableColumnName = get(handles.uitable1,'ColumnName');
configTable = cell2table(configTable);
configTable.Properties.VariableNames = configTableColumnName;
writetable(configTable,handles.sessionConfigName);

% --- Executes on key press with focus on listbox1 and none of its controls.
function listbox1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in rootButton.
function rootButton_Callback(hObject, eventdata, handles)
% hObject    handle to rootButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
configRoot = uigetdir(pwd);
set(handles.rootText,'String', configRoot)
handles.configRoot = configRoot;
handles.savePath = configRoot;
guidata(hObject,handles);


% --- Executes on button press in savePathButton.
function savePathButton_Callback(hObject, eventdata, handles)
% hObject    handle to savePathButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savePath = uigetdir(pwd);
set(handles.savePathText,'String', savePath)
handles.savePath = savePath;
guidata(hObject,handles);



function editAnimal_Callback(hObject, eventdata, handles)
% hObject    handle to editAnimal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAnimal as text
%        str2double(get(hObject,'String')) returns contents of editAnimal as a double


% --- Executes during object creation, after setting all properties.
function editAnimal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAnimal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in preprocessRunButton.
function preprocessRunButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessRunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% NOTE: CURRENTLY EDITING PARAM DOES NOT CHANGE PARAM THAT RUNS THE
% PROGRAMS, NEED CHANGE IN THE FUTURE.
disp(handles.param)
disp(handles.configRoot)
disp(handles.savePath)
allStr = get(handles.preprocessMenu,'String');
currStr = allStr{get(handles.preprocessMenu,'Value')};
switch currStr
case 'Tuning Curve'
    getTuning(handles.param{:});
end

% --- Executes on button press in preprocessParamEditButton.
function preprocessParamEditButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessParamEditButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of preprocessParamEditButton
if hObject.Value==1
    set(handles.uitable2,'ColumnEditable',true(1,length(get(handles.uitable2,'ColumnName'))));
else
    set(handles.uitable2,'ColumnEditable',false(1,length(get(handles.uitable2,'ColumnName'))));
end
guidata(hObject,handles);
