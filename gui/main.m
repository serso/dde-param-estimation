function varargout = main(varargin)
% MAIN M-file for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 17-Jun-2012 02:35:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @main_OpeningFcn, ...
    'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

addpath('../');
% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fEdit_Callback(hObject, eventdata, handles)
% hObject    handle to fEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fEdit as text
%        str2double(get(hObject,'String')) returns contents of fEdit as a double


% --- Executes during object creation, after setting all properties.
function fEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delayCountEdit_Callback(hObject, eventdata, handles)
% hObject    handle to delayCountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delayCountEdit as text
%        str2double(get(hObject,'String')) returns contents of delayCountEdit as a double
[delayCount, error] = checkCount(hObject, handles.prevDelayCount, handles, 0);

if ( ~error )
    handles.prevDelayCount = num2str(delayCount);
    guidata(hObject, handles);
    adjustTable(handles.delaysTable, delayCount);
end

% --- Executes during object creation, after setting all properties.
function delayCountEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delayCountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.prevDelayCount = '1';
set(hObject, 'String', handles.prevDelayCount);
guidata(hObject, handles);




function pCountEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pCountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pCountEdit as text
%        str2double(get(hObject,'String')) returns contents of pCountEdit as a double
[pCount, error] = checkCount(hObject, handles.prevPCount, handles, 1);
if ( ~error )
    handles.prevPCount = num2str(pCount);
    guidata(hObject, handles);
    adjustTable(handles.paramsTable, pCount);
    adjustTable(handles.paramsResultTable, pCount);
end


function [count, error] = checkCount(hObject, prevCount, handles, minCount)
count = str2double(get(hObject,'String'));

error = false;

if ( count ~= floor(count) )
    set(hObject, 'String', prevCount);
    guidata(hObject, handles);
    msgbox('Number must be an integer value!', 'Incorrect number','error');
    error = true;
end

if ( count < minCount )
    set(hObject, 'String', prevCount);
    guidata(hObject, handles);
    msgbox(sprintf('Number must be not less than %i!', minCount), 'Incorrect number','error');
    error = true;
end


function adjustTable(table, newRows)
data = get(table,'Data');
[dataRows, dataCols] = size(data);
if ( newRows > dataRows )
    data = cat(1, data, cell(newRows - dataRows, dataCols));
elseif ( newRows < dataRows )
    data = data(1:newRows,:);
end
set(table, 'Data', data);


% --- Executes during object creation, after setting all properties.
function pCountEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pCountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.prevPCount = '1';
set(hObject, 'String', handles.prevPCount);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function delaysTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delaysTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', {2});


% --- Executes during object creation, after setting all properties.
function paramsTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paramsTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(1, 1));



function dataCountEdit_Callback(hObject, eventdata, handles)
% hObject    handle to dataCountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataCountEdit as text
%        str2double(get(hObject,'String')) returns contents of dataCountEdit as a double
[dataCount, error] = checkCount(hObject, handles.prevDataCount, handles, 10);

if ( ~error )
    handles.prevDataCount = num2str(dataCount);
    guidata(hObject, handles);
    adjustTable(handles.dataTable, dataCount);
end

% --- Executes during object creation, after setting all properties.
function dataCountEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataCountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.prevDataCount = '100';
set(hObject, 'String', handles.prevDataCount);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function dataTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(100, 2));



function constCountEdit_Callback(hObject, eventdata, handles)
% hObject    handle to constCountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of constCountEdit as text
%        str2double(get(hObject,'String')) returns contents of constCountEdit as a double
[constCount, error] = checkCount(hObject, handles.prevConstCount, handles, 0);

if ( ~error )
    handles.prevConstCount = num2str(constCount);
    guidata(hObject, handles);
    adjustTable(handles.constantsTable, constCount);
end

% --- Executes during object creation, after setting all properties.
function constCountEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to constCountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.prevConstCount = '1';
set(hObject, 'String', handles.prevConstCount);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function loadFromTableRadio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadFromTableRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function dataSource_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dataSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in dataSource.
function dataSource_SelectionChangeFcn(hObject, eventdata, handles)
switch get(eventdata.NewValue, 'Tag') % Get Tag of selected object.
    case 'loadFromTableRadio'
        toggleDataSource(hObject, handles, 'fromTable');
    case 'loadFromFileRadio'
        toggleDataSource(hObject, handles, 'fromFile');
        
        % Code for when radiobutton2 is selected.
    otherwise
        % Code for when there is no match.
end
% hObject    handle to the selected object in dataSource
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
function toggleDataSource(hObject, handles, source)
switch source
    case 'fromTable'
        set(handles.chooseDataFileButton, 'Enable', 'off');
        set(handles.dataTable, 'Enable', 'on');
        set(handles.dataCountEdit, 'Enable', 'on');
        guidata(hObject, handles);
    case 'fromFile'
        set(handles.chooseDataFileButton, 'Enable', 'on');
        set(handles.dataTable, 'Enable', 'off');
        set(handles.dataCountEdit, 'Enable', 'off');
        guidata(hObject, handles);
end

% --- Executes on button press in chooseDataFileButton.
function chooseDataFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to chooseDataFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[dataFileName, dataFilePath] = uigetfile('*.csv', 'Choose cvs file with data');

if ( dataFileName )
    handles.dataFileName = dataFileName;
    handles.dataFilePath = dataFilePath;
    set(handles.dataFileNameText, 'String', handles.dataFileName);
    set(handles.dataFilePathText, 'String', handles.dataFilePath);
    guidata(hObject, handles);
    
    checkProblemData(handles);
end


% --- Executes during object creation, after setting all properties.
function dataSource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function afterAllText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to afterAllText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function chooseDataFileButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chooseDataFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over chooseDataFileButton.
function chooseDataFileButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to chooseDataFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function startFromRowEdit_Callback(hObject, eventdata, handles)
% hObject    handle to startFromRowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFromRowEdit as text
%        str2double(get(hObject,'String')) returns contents of startFromRowEdit as a double
[startFromRow, error] = checkCount(hObject, handles.startFromRow, handles, 1);
if ( ~error )
    [error] = checkProblemData(handles);
    if ( ~error )
        handles.startFromRow = num2str(startFromRow);
        guidata(hObject, handles);
    end
end


% --- Executes during object creation, after setting all properties.
function startFromRowEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFromRowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.startFromRow = '1';
set(hObject, 'String', handles.startFromRow);
guidata(hObject, handles);


function startFromColumnEdit_Callback(hObject, eventdata, handles)
% hObject    handle to startFromColumnEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFromColumnEdit as text
%        str2double(get(hObject,'String')) returns contents of startFromColumnEdit as a double
[startFromColumn, error] = checkCount(hObject, handles.startFromColumn, handles, 1);
if ( ~error )
    [error] = checkProblemData(handles);
    if ( ~error )
        handles.startFromColumn = num2str(startFromColumn);
        guidata(hObject, handles);
    end
end


% --- Executes during object creation, after setting all properties.
function startFromColumnEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFromColumnEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.startFromColumn = '1';
set(hObject, 'String', handles.startFromColumn);
guidata(hObject, handles);

function [error, fStr, np, ndelays, userConstants, xHistoryStr, tData, xData, delays, p0] = checkProblemDef (handles)

error = true;
fStr = '';
np = 0;
ndelays = 0;
userConstants = {};

try
    
    fStr = get(handles.fEdit, 'String');
    xHistoryStr = get(handles.xHistoryEdit, 'String');
    np = str2double(get(handles.pCountEdit, 'String'));
    ndelays = str2double(get(handles.delayCountEdit, 'String'));
    userConstantsCount = str2double(get(handles.constCountEdit, 'String'));
    
    checkTable(handles.delaysTable, 'delays table');
    if ( ndelays > 0 )
        delays = cell2mat(get(handles.delaysTable, 'Data'))';
    else
        delays = [];
    end

    checkTable(handles.constantsTable, 'constants table');
    
    try
        checkTable(handles.paramsTable, 'parameters table');
        p0 = cell2mat(get(handles.paramsTable, 'Data'));
    catch e
        if ( strcmp(e.identifier, 'IllegalArgumentException:NotEnoughTableData'))
            p0 = [];
        else
            rethrow(e);
        end
    end
    
    userConstantsM = get(handles.constantsTable, 'Data');
    
    [rows, cols] = size(userConstantsM);
    for row = 1:min([rows, userConstantsCount])
        for col = 1:cols
            userConstants{row}{col} = userConstantsM(row, col);
        end
    end
    
    ddeParamEst_checkUserInput( fStr, np, ndelays, userConstants, xHistoryStr );
    
    [error, tData, xData] = checkProblemData(handles);    
catch e
    msgbox(e.message, e.identifier,'error');
    if ( ~strcmp(e.identifier, 'IllegalArgumentException:NotEnoughTableData') && ...
            ~strcmp(e.identifier, 'ArgumentCheck:IllegalArgument'))
        rethrow(e);
    end
end

function [error, tData, xData] = checkProblemData (handles)

error = true;
data = [];
tData = [];
xData = [];

try
    if ( get(handles.loadFromFileRadio, 'Value') )
        
        dataFileName = get(handles.dataFileNameText, 'String');
        if ( isempty(dataFileName) )
            throw (MException ('ArgumentCheck:IllegalArgument', 'Choose the file!'))
        end
        
        dataFilePath = get(handles.dataFilePathText, 'String');
        if ( isempty(dataFilePath) )
            throw (MException ('ArgumentCheck:IllegalArgument', 'Choose the file!'))
        end
        
        fileTaskName = strcat(dataFilePath, dataFileName);
        
        startRow = str2double(get(handles.startFromRowEdit, 'String')) - 1;
        startColumn = str2double(get(handles.startFromColumnEdit, 'String')) - 1;
        
        data = csvread(fileTaskName, startRow, startColumn);
        
    elseif (get(handles.loadFromTableRadio, 'Value'))
        data = cell2mat(get(handles.dataTable, 'Data'));
    else
        throw (MException ('ArgumentCheck:IllegalArgument', 'Choose either "Load from file" or "Load from table"!'))
    end
    
    tData = data(1:end,1);
    xData = data(1:end,2);
    
    if ( length(tData) ~= length(xData) )
        throw (MException ('ArgumentCheck:IllegalArgument', 't and x must be the same length!'))
    end

    plot(handles.resultFigure, tData, xData, 'xr');
    set(handles.resultFigure,'XGrid','on');
    set(handles.resultFigure,'YGrid','on');
    %set(handles.resultFigure, 'XData', tData);
    %set(handles.resultFigure, 'YData', xData);
    
    %guidata(hObject, handles);
    
    error = false;
catch e
    msgbox(e.message, e.identifier,'error');
    if ( ~strcmp(e.identifier, 'ArgumentCheck:IllegalArgument'))
        rethrow(e);
    end
end



% --- Executes on button press in solveButton.
function solveButton_Callback(hObject, eventdata, handles)
% hObject    handle to solveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[error, fStr, np, ndelays, userConstants, xHistoryStr, tData, xData, delays, p0] = checkProblemDef(handles);
if ( ~error )
    [ f, fg, fh, ~, ~, ~, xHistory, ~] = ddeParamEst_checkUserInput( fStr, np, ndelays, userConstants, xHistoryStr);
    
    options = [];
    options.optOptions = [];
    options.sqpOptions = [];
    options.hessian_method = 'gauss-newton';
    options.sqpOptions.stepMethodIterative = false;
    options.sqpOptions.stepMethod = 'ldls';
    options.pTol = 0.1;
    options.plotResult = true;
    options.showResult = true;
    options.taskName = fStr;
    options.plotExtResult = true;
    options.extTMax = 2300;
    options.extFigureHandle = handles.resultFigure;
    
    
    [x, p, info] = ddeParamEst(...
        tData, xData, ...
        f, fg, fh, ...
        delays, xHistory, max(delays), ...
        options, np, ...
        [], [], p0);
    
    paramsData = get(handles.paramsResultTable,'Data');
    for i = length(p)
        paramsData{i} = p(i);
        set(handles.paramsResultTable,'Data', paramsData);
        guidata(hObject, handles);
    end
    
end

function checkTable(table, tableName)

m = get(table, 'Data');

[rows, cols] = size(m);

for row = 1:rows
    for col = 1:cols
        tmp = m(row, col);
        if ( isempty(tmp{1}) )
            throw(MException('IllegalArgumentException:NotEnoughTableData', sprintf('Not enough data is entered to table %s', tableName)));
        end
    end
end

% --- Executes on button press in computeSymbolicData.
function computeSymbolicData_Callback(hObject, eventdata, handles)
% hObject    handle to computeSymbolicData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[error, fStr, np, ndelays, userConstants, xHistoryStr] = checkProblemDef(handles);
if ( ~error )
    [ ~, ~, ~, fSym, fgSym, fhSym] = ddeParamEst_checkUserInput( fStr, np, ndelays, userConstants, xHistoryStr);
    
    fprintf('f = \n');
    pretty(fSym);
    fprintf('fg = \n');
    pretty(fgSym);
    fprintf('fh = \n');
    pretty(fhSym);
    
    set(handles.fgSymText, 'String', utils.symToStr(fgSym));
    set(handles.fhSymText, 'String', utils.symToStr(fhSym));
    guidata(hObject, handles);
    
end



function xHistoryEdit_Callback(hObject, eventdata, handles)
% hObject    handle to xHistoryEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xHistoryEdit as text
%        str2double(get(hObject,'String')) returns contents of xHistoryEdit as a double


% --- Executes during object creation, after setting all properties.
function xHistoryEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xHistoryEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function paramsResultTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paramsResultTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(1, 2));


% --- Executes when entered data in editable cell(s) in dataTable.
function dataTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to dataTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%checkProblemData(handles);



function dataDelimiterEdit_Callback(hObject, eventdata, handles)
% hObject    handle to dataDelimiterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataDelimiterEdit as text
%        str2double(get(hObject,'String')) returns contents of dataDelimiterEdit as a double


% --- Executes during object creation, after setting all properties.
function dataDelimiterEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataDelimiterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
