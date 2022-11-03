
function varargout = preprocessMasterGUI(varargin)
% PREPROCESSMASTERGUI MATLAB code for preprocessMasterGUI.fig
%      PREPROCESSMASTERGUI, by itself, creates a new PREPROCESSMASTERGUI or raises the existing
%      singleton*.
%
%      H = PREPROCESSMASTERGUI returns the handle to a new PREPROCESSMASTERGUI or the handle to
%      the existing singleton*.
%
%      PREPROCESSMASTERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESSMASTERGUI.M with the given input arguments.
%
%      PREPROCESSMASTERGUI('Property','Value',...) creates a new PREPROCESSMASTERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preprocessMasterGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preprocessMasterGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preprocessMasterGUI

% Last Modified by GUIDE v2.5 28-Oct-2021 22:57:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preprocessMasterGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @preprocessMasterGUI_OutputFcn, ...
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


% --- Executes just before preprocessMasterGUI is made visible.
function preprocessMasterGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preprocessMasterGUI (see VARARGIN)

% Choose default command line output for preprocessMasterGUI
handles.output = hObject;
% Initialize other parameters
handles.fn1 = '';handles.fn2 = '';handles.fn3 = '';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes preprocessMasterGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = preprocessMasterGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% Set function selection menu; detect main functions from preprocess/main folder
guiDir = which('preprocessMasterGUI');guiDirSplit = strsplit(guiDir,filesep);
paramFnDir = [strjoin(guiDirSplit(1:end-2),filesep) filesep 'param_function'];
mainDir = [strjoin(guiDirSplit(1:end-2),filesep) filesep 'main']; mFile = dir([mainDir filesep '*.m']);
if isempty(mFile); msgbox('WARNING -- NO MAIN FUNCTION DETECTED!'); end
mFile = {mFile.name}'; mFile = cellfun(@(x)(x(1:end-2)),mFile,'UniformOutput',false);mFile=['none';mFile];
handles.fnSelectionMenu.String = mFile;
handles.mainDir = mainDir; handles.paramFnDir = paramFnDir;
guidata(hObject,handles);

% --- Executes on button press in preprocessStatusSetButton.
function preprocessStatusSetButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessStatusSetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
targetVariable = get(handles.preprocessStatusMenu,'String');
targetVariable = targetVariable{get(handles.preprocessStatusMenu,'Value')};
if ~strcmp(targetVariable, 'None')
    configTable = get(handles.mouseInfoTable,'Data');
    configTableColumnName = get(handles.mouseInfoTable,'ColumnName');
    configTable = cell2table(configTable);
    configTable.Properties.VariableNames = configTableColumnName;

    configTable.(targetVariable) = ones(size(configTable,1),1);
    set(handles.mouseInfoTable,'Data', table2cell(configTable));
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
handles.sessionConfigName = [handles.configRoot filesep 'mouse' filesep masterInfo.sessionConfig];
handles.imagingConfigName = [handles.configRoot filesep 'imaging' filesep masterInfo.imagingConfig];
handles.behaviorConfigName = [handles.configRoot filesep 'behavior' filesep masterInfo.behaviorConfig];
configTable = readtable(handles.sessionConfigName);
set(handles.mouseInfoTable,'Data',table2cell(configTable));
set(handles.mouseInfoTable,'ColumnName', configTable.Properties.VariableNames);
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


% --- Executes on selection change in fileListBox.
function fileListBox_Callback(hObject, eventdata, handles)
% hObject    handle to fileListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.listboxSelectedIndex = eventdata.Source.Value;
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns fileListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fileListBox


% --- Executes during object creation, after setting all properties.
function fileListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fnSelectionMenu.
function fnSelectionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fnSelectionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fnSelectionMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fnSelectionMenu
allOpt = get(hObject,'String');
selectedOpt = allOpt{eventdata.Source.Value};
configTable = getConfigTable(handles);
% Load filenames for file selection 
filenames = fn_fileListAutoSelect(configTable,selectedOpt);

set(handles.fileListBox,'String',filenames);

% --- Executes during object creation, after setting all properties.
function fnSelectionMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fnSelectionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function configTable = getConfigTable(handles)
configTable = get(handles.mouseInfoTable,'Data');
configTableColumnName = get(handles.mouseInfoTable,'ColumnName');
configTable = cell2table(configTable);
configTable.Properties.VariableNames = configTableColumnName;


% --- Executes on button press in paramPreviewButton.
%---------GET PARAMETERS SHARED BY ALL FUNCTIONS, IMAGING AND MOUSE CONFIG-----------
function paramPreviewButton_Callback(hObject, eventdata, handles)
% hObject    handle to paramPreviewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% LOAD config table
configTable = getConfigTable(handles);
configTable = func_fillMissingTable(configTable);
configTable = fn_correctConfigTableMisc(configTable);
% READ and check filenames
filenames = get(handles.fileListBox,'String');if ~iscell(filenames);filenames = {filenames};end
filenameIdx = cellfun(@(x)(find(strcmp(configTable.ImagingFile,x))),filenames);
fn_paramPreviewFileNumCheck(handles,filenames);

% add the imaging configuration into parameters
paramValue = {};

tempFieldName = fieldnames(handles.masterInfo);
for i = 1:length(fieldnames(handles.masterInfo)); paramValue = [paramValue;tempFieldName{i};handles.masterInfo.(tempFieldName{i})];end %#ok<*AGROW>

imagingConfig = func_loadImagingConfig(handles.imagingConfigName,'Split', true);
tempFieldName = fieldnames(imagingConfig);
for i = 1:length(tempFieldName); paramValue = [paramValue;tempFieldName{i};{imagingConfig.(tempFieldName{i})}]; end

% save the root and savepath parameters
paramValue = [paramValue;'root';handles.configRoot;'savepath'; handles.savePath];
paramValueDisplay = cellfun(@(x)(fn_cellStrJoin(x)),paramValue,'UniformOutput',false);

% get the filename selected, and the table of information
paramValue = [paramValue;'filename';strjoin(filenames);'filenameIdx';filenameIdx;'nFiles';length(filenames)];
[tempFilenameFlag,filenameFrames,fileMultiLoadFlag] = fn_paramPreviewFileSameFlag(configTable,filenames,imagingConfig,'tcFile');
paramValue = [paramValue;'filenameTCIdx';tempFilenameFlag;'fileMultiLoadTCFlag';fileMultiLoadFlag;'nFrames_oneplane_TC';{filenameFrames}];
%paramValue = [paramValue;'filenameTC';strjoin(filenameSel);'filenameTCFlag';tempFilenameFlag;'nFrames_oneplane_TC';filenameFrames];
[tempFilenameFlag,filenameFrames,fileMultiLoadFlag] = fn_paramPreviewFileSameFlag(configTable,filenames,imagingConfig,'suite2ppath');
paramValue = [paramValue;'filenameBinIdx';tempFilenameFlag;'fileMultiLoadTCFlag';fileMultiLoadFlag;'nFrames_oneplane_bin';{filenameFrames}];
%paramValue = [paramValue;'filenameBin';strjoin(filenameSel);'filenameBinFlag';tempFilenameFlag;'nFrames_oneplane_bin';filenameFrames];

% get all the path and files needed
paramValue = fn_paramPreviewAddTableParam(configTable,filenameIdx,paramValue,imagingConfig);

% display parameters and save it to the handle
paramTableDisplay = table(paramValueDisplay(1:2:end),paramValueDisplay(2:2:end),'VariableNames',{'name','value'});
set(handles.paramPreviewTable,'Data',table2cell(paramTableDisplay));
set(handles.paramPreviewTable,'ColumnName', paramTableDisplay.Properties.VariableNames);
handles.param = paramValue; 
guidata(hObject,handles);

% --- Executes on button press in fileAddButton.
function fileAddButton_Callback(hObject, eventdata, handles)
% hObject    handle to fileAddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listBoxStr = get(handles.fileListBox,'String');
if ischar(listBoxStr)
    listBoxStr = {listBoxStr};
end
configTable = get(handles.mouseInfoTable,'Data');
for i = 1:size(handles.uitableSelectedIndex,1)
    tempFilename = configTable{handles.uitableSelectedIndex(i,1), handles.uitableSelectedIndex(i,2)};
    if ~any(strcmp(listBoxStr, tempFilename))
        listBoxStr = [listBoxStr; tempFilename];
    end
end
%if strcmp(listBoxStr,'SelectFile') || strcmp(listBoxStr{1},'SelectFile')
listBoxStr(strcmp(listBoxStr,'SelectFile')) = [];
%end
set(handles.fileListBox,'String',listBoxStr);

% --- Executes on button press in fileDropButton.
function fileDropButton_Callback(hObject, eventdata, handles)
% hObject    handle to fileDropButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listBoxStr = get(handles.fileListBox,'String');
if ischar(listBoxStr)
   listBoxStr = {listBoxStr}; 
end
selectBoxStr = listBoxStr(get(handles.fileListBox,'Value'));

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

set(handles.fileListBox,'String',listBoxStr);

% --- Executes on button press in preprocessStatusClearButton.
function preprocessStatusClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessStatusClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
targetVariable = get(handles.preprocessStatusMenu,'String');
targetVariable = targetVariable{get(handles.preprocessStatusMenu,'Value')};
if ~strcmp(targetVariable, 'None')
    configTable = get(handles.mouseInfoTable,'Data');
    configTableColumnName = get(handles.mouseInfoTable,'ColumnName');
    configTable = cell2table(configTable);
    configTable.Properties.VariableNames = configTableColumnName;

    configTable.(targetVariable) = zeros(size(configTable,1),1);
    set(handles.mouseInfoTable,'Data', table2cell(configTable));
end

% --- Executes when selected cell(s) is changed in mouseInfoTable.
function mouseInfoTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to mouseInfoTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
handles.uitableSelectedIndex = eventdata.Indices;
guidata(hObject,handles);



% --- Executes when entered data in editable cell(s) in mouseInfoTable.
function mouseInfoTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to mouseInfoTable (see GCBO)
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
    set(handles.mouseInfoTable,'ColumnEditable',true(1,length(get(handles.mouseInfoTable,'ColumnName'))));
else
    set(handles.mouseInfoTable,'ColumnEditable',false(1,length(get(handles.mouseInfoTable,'ColumnName'))));
end
guidata(hObject,handles);

% --- Executes on button press in preprocessStatusSaveButton.
function preprocessStatusSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessStatusSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
configTable = get(handles.mouseInfoTable,'Data');
configTableColumnName = get(handles.mouseInfoTable,'ColumnName');
configTable = cell2table(configTable);
configTable.Properties.VariableNames = configTableColumnName;
writetable(configTable,handles.sessionConfigName);

% --- Executes on key press with focus on fileListBox and none of its controls.
function fileListBox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fileListBox (see GCBO)
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
currStr = getSelFnName(handles);
param = handles.param;
if ~isempty(handles.fn1); param = [param;'fn1';handles.fn1]; end
if ~isempty(handles.fn2); param = [param;'fn2';handles.fn2]; end
if ~isempty(handles.fn3); param = [param;'fn3';handles.fn3]; end
feval(currStr,param{:});

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over fnSelectionMenu.
function fnSelectionMenu_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to fnSelectionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function x = fn_cellStrJoin(x)
if iscell(x); x = strjoin(x);end

% --- Executes on button press in fnButton1.
function fnButton1_Callback(hObject, eventdata, handles)
% hObject    handle to fnButton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fnButton1
paramValue = handles.param;
if hObject.Value==1
    [funcParam, funcParamPath] = uigetfile(handles.paramFnDir);
    addpath(funcParamPath);tempSplit = strsplit(funcParam,'.');
    set(handles.fnButton1,'String', tempSplit{1}); handles.fn1 = tempSplit{1};
    if length(tempSplit{1})>10; set(handles.fnButton1,'FontSize', 6.5); end
else
    set(handles.fnButton1,'String', 'fn1');set(handles.fnButton1,'FontSize', 8.0);handles.fn1 = '';
end
guidata(hObject,handles);


% --- Executes on button press in fnButton2.
function fnButton2_Callback(hObject, eventdata, handles)
% hObject    handle to fnButton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fnButton2
paramValue = handles.param;
if hObject.Value==1
    [funcParam, funcParamPath] = uigetfile(handles.paramFnDir);
    addpath(funcParamPath);tempSplit = strsplit(funcParam,'.');
    set(handles.fnButton2,'String', tempSplit{1}); handles.fn2 = tempSplit{1};
     if length(tempSplit{1})>10; set(handles.fnButton2,'FontSize', 6.5); end
else
    set(handles.fnButton2,'String', 'fn2');set(handles.fnButton2,'FontSize', 8.0);handles.fn2 = '';
end
guidata(hObject,handles);


% --- Executes on button press in fnButton3.
function fnButton3_Callback(hObject, eventdata, handles)
% hObject    handle to fnButton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fnButton3
paramValue = handles.param;
if hObject.Value==1
    [funcParam, funcParamPath] = uigetfile(handles.paramFnDir);
    addpath(funcParamPath);tempSplit = strsplit(funcParam,'.');
    set(handles.fnButton3,'String', tempSplit{1}); handles.fn3 = tempSplit{1};
    if length(tempSplit{1})>10; set(handles.fnButton3,'FontSize', 6.5); end
else
    set(handles.fnButton3,'String', 'fn3'); set(handles.fnButton3,'FontSize', 8.0); handles.fn3 = '';
end
guidata(hObject,handles);


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of helpButton
msg = help(getSelFnName(handles)); mbox = msgbox(msg);    
msgboxFontSize(mbox,11,'FontName','Consolas');     


function currStr = getSelFnName(handles)
allStr = get(handles.fnSelectionMenu,'String');
currStr = allStr{get(handles.fnSelectionMenu,'Value')};


