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

% Last Modified by GUIDE v2.5 04-Jul-2020 13:39:50

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
handles.sep = '\';
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
    case 'extractTC'
        filenames = configTable.('ImagingFile');
    case 'extractRedraw'
        filenames = configTable.('ImagingFile');
    case 'deconvolve'
        filenames = configTable.('ImagingFile');
    case 'roiTracking'
        filenames = configTable.('ImagingFile');
    case 'roiRevision'
        filenames = configTable.('ImagingFile');
    case 'roiRevisionPlot'
        filenames = configTable.('ImagingFile');
    case 'getTuning'
        filenames = 'SelectFile';
    case 'none'
        filenames = 'SelectFile';
    case 'saveMeanImg'
        filenames = configTable.('ImagingFile');
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
%---------GET PARAMETERS SHARED BY ALL FUNCTIONS, IMAGING AND MOUSE CONFIG-----------
configTable = getConfigTable(handles);
configTable = func_fillMissingTable(configTable);

% add the imaging configuration into parameters
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

allStr = get(handles.preprocessMenu,'String');
currStr = allStr{get(handles.preprocessMenu,'Value')};


% get the filename selected, and the table of information
filenames = get(handles.listbox1,'String');if ~iscell(filenames);filenames = {filenames};end
filenameFlag = strcmp(configTable.ImagingFile,filenames);


binSelected = configTable.suite2ppath(filenameFlag);
%if length(tcFileSelected)>1 && ~isequal(tcFileSelected{:}); disp('ERROR: Selected files does not have same tcFile!'); end
binSelected_removeNone = binSelected(~strcmp(binSelected,'None'));
[~,~,j]=unique(binSelected_removeNone);
binCommon = binSelected_removeNone{mode(j)}; % get all tcfile name of selected files
% get all the imaging files with same tcFile name with selected filenames
filename_binFlagCommon = strcmp(configTable.suite2ppath,binCommon);
% nFrames_oneplane_load is used to load TC file corresponding to selected files 
nFrames_oneplane_binStr = configTable.nFrames_oneplane(filename_binFlagCommon);
% filename_loadFlag is used to identify where selected files 
filename_binFlag = filename_binFlagCommon; filename_binFlag(~filenameFlag) = 0; filename_binFlag = filename_binFlag(filename_binFlagCommon);

if ~iscell(nFrames_oneplane_binStr); nFrames_oneplane_binStr = {nFrames_oneplane_binStr}; end
for i = 1:length(nFrames_oneplane_binStr)
    nFrames_oneplane_binStrSplit = strsplit(nFrames_oneplane_binStr{i});
    for j = 1:str2double(imagingConfig.nPlanes)
        nFrames_oneplane_bin(i,j) = str2double(nFrames_oneplane_binStrSplit{j});
    end
end

filename_bin = configTable.ImagingFile(filename_binFlag);
paramValue = [paramValue;'filenameBin';strjoin(filename_bin);'filenameBinFlag';filename_binFlag;'nFrames_oneplane_bin';nFrames_oneplane_bin];


tcFileSelected = configTable.tcFile(filenameFlag);
tcFileSelected_removeNone = tcFileSelected(~strcmp(tcFileSelected,'None'));
[~,~,j]=unique(tcFileSelected_removeNone);
tcFilenameCommon = tcFileSelected{mode(j)}; % get all tcfile name of selected files
% get all the imaging files with same tcFile name with selected filenames
filename_TCFlag = strcmp(configTable.tcFile,tcFilenameCommon);
% nFrames_oneplane_load is used to load TC file corresponding to selected files 
nFrames_oneplane_loadStr = configTable.nFrames_oneplane(filename_TCFlag);
% filename_loadFlag is used to identify where selected files 
filename_loadFlag = filename_TCFlag; filename_loadFlag(~filenameFlag) = 0; filename_loadFlag = filename_loadFlag(filename_TCFlag);

if ~iscell(nFrames_oneplane_loadStr); nFrames_oneplane_loadStr = {nFrames_oneplane_loadStr}; end
for i = 1:length(nFrames_oneplane_loadStr)
    nFrames_oneplane_loadStrSplit = strsplit(nFrames_oneplane_loadStr{i});
    for j = 1:str2double(imagingConfig.nPlanes)
        nFrames_oneplane_load(i,j) = str2double(nFrames_oneplane_loadStrSplit{j});
    end
end

filename_TC = configTable.ImagingFile(filename_TCFlag);
paramValue = [paramValue;'filenameTC';strjoin(filename_TC);'filenameTCFlag';filename_loadFlag;'nFrames_oneplane_TC';nFrames_oneplane_load];

% save the root and savepath parameters
param = [param;'root';'savepath'];paramValueDisplay = [paramValueDisplay;handles.configRoot; handles.savePath];
paramValue = [paramValue;'filename';strjoin(filenames);'root';handles.configRoot;'savepath'; handles.savePath];
% get all the path and files needed
varNames = configTable.Properties.VariableNames;
for i = 1:length(varNames)
    tempParam = configTable.(varNames{i})(filenameFlag);
    
    if strcmp(varNames{i},'nFrames_oneplane')
        saveParam = [];
        if ~iscell(tempParam); tempParam = {tempParam}; end
        for k = 1:length(tempParam)
            nFrames_oneplaneSplit = strsplit(tempParam{k});
            for j = 1:str2double(imagingConfig.nPlanes)
                saveParam(k,j) = str2double(nFrames_oneplaneSplit{j});
            end
        end  
    else
        saveParam = tempParam;
    end
    
    if iscell(saveParam) && length(saveParam) > 1
        
        if isequal(saveParam{:}); saveParam = saveParam{1};
        else;saveParam = strjoin(saveParam);
            if ~strcmp(varNames{i},'BehavFile') && ~strcmp(varNames{i},'ImagingFile') && ~strcmp(varNames{i},'BehavType') 
                disp(['WARNING: Not all files have same ' varNames{i} '!'])
            end
        end
    elseif isnumeric(saveParam)
        saveParam = {saveParam};
    end
    paramValue = [paramValue;varNames{i};saveParam];

    displayString = {'sbxpath';'h5path';'suite2ppath';'root';'savepath';'datapath';'behavpath'};
    if any(contains(displayString,varNames{i}))
        param = [param;varNames{i}];paramValueDisplay = [paramValueDisplay;saveParam];
    end
end

switch currStr
    case 'getTuning'
        if length(filenames) ~=1 
            disp('ERROR: Please select only one file');
        end
        filename = filenames{1};
    case 'extractTC'   
    case 'extractRedraw'    
    case 'roiTracking'
    case 'roiRevision'        
    case 'saveMeanImg'      
end
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
allStr = get(handles.preprocessMenu,'String');
currStr = allStr{get(handles.preprocessMenu,'Value')};
switch currStr
case 'getTuning'
    getTuning(handles.param{:});
case 'extractTC'
    extractTC(handles.param{:});
case 'extractRedraw'
    extractRedraw(handles.param{:});
case 'saveMeanImg'
    getMeanImg(handles.param{:});
case 'roiTracking'
    roiTracking(handles.param{:});
case 'roiRevision'
    roiRevision(handles.param{:});
case 'roiRevisionPlot'
    roiRevisionPlot(handles.param{:});
end

% --- Executes on button press in funcParam.
function funcParam_Callback(hObject, eventdata, handles)
% hObject    handle to funcParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of funcParam
%if hObject.Value==1
%    set(handles.uitable2,'ColumnEditable',true(1,length(get(handles.uitable2,'ColumnName'))));
%else
%    set(handles.uitable2,'ColumnEditable',false(1,length(get(handles.uitable2,'ColumnName'))));
%end
%guidata(hObject,handles);
paramValue = handles.param;
if hObject.Value==1
    [funcParam, funcParamPath] = uigetfile(handles.configRoot);
    set(handles.funcParamText,'String', funcParam); addpath(funcParamPath);
    tempSplit = strsplit(funcParam,'.');
    funcParamIndex = find(strcmp(paramValue, 'funcParam'));
    if isempty(funcParamIndex); paramValue = [paramValue; 'funcParam'; tempSplit{1}];
    else; paramValue{funcParamIndex+1} = tempSplit{1};end
else
    funcParamIndex = find(strcmp(paramValue, 'funcParam'));
    set(handles.funcParamText,'String', 'None');
    if ~isempty(funcParamIndex);paramValue{funcParamIndex:funcParamIndex+1} = [];end
end
handles.param = paramValue; 
guidata(hObject,handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over preprocessMenu.
function preprocessMenu_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to preprocessMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over funcParam.
function funcParam_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to funcParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
