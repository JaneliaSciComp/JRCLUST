function varargout = jrc_gui(varargin)
% JRC_GUI MATLAB code for jrc_gui.fig
%      JRC_GUI, by itself, creates a new JRC_GUI or raises the existing
%      singleton*.
%
%      H = JRC_GUI returns the handle to a new JRC_GUI or the handle to
%      the existing singleton*.
%
%      JRC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JRC_GUI.M with the given input arguments.
%
%      JRC_GUI('Property','Value',...) creates a new JRC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before jrc_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to jrc_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help jrc_gui

% Last Modified by GUIDE v2.5 29-Sep-2017 14:42:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @jrc_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @jrc_gui_OutputFcn, ...
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


% --- Executes just before jrc_gui is made visible.
function jrc_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to jrc_gui (see VARARGIN)

% Background image
% image(handles.axes1, imread('./img/gui_background.png'));
% set(handles.axes1, 'XTick', [], 'YTick', [], 'Visible', 'off');
% axis(handles.axes1, 'equal');

% Disable buttons and menus
enable_(handles, {'Spikesort', 'Detect', 'Sort', 'ClearPrm', 'EditPrm', ...
    'Auto', 'Probe', 'Preview', 'Traces', 'Manual', 'Describe', 'Drift'}, 0);

% Choose default command line output for jrc_gui
handles.output = hObject;
if ~isempty(varargin{1})    
    handles.vcFile_prm = load_prm_(hObject, varargin{1});
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes jrc_gui wait for user response (see UIRESUME)
% uiwait(handles.Fig_gui);


% --- Outputs from this function are returned to the command line.
function varargout = jrc_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ui_file_Callback(hObject, eventdata, handles)
% hObject    handle to ui_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_18_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuWiki_Callback(hObject, eventdata, handles)
% hObject    handle to menuWiki (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuGitPull_Callback(hObject, eventdata, handles)
% hObject    handle to menuGitPull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to menuUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuAbout_Callback(hObject, eventdata, handles)
% hObject    handle to menuAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuIssue_Callback(hObject, eventdata, handles)
% hObject    handle to menuIssue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuTraces_Callback(hObject, eventdata, handles)
% hObject    handle to menuTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuPreview_Callback(hObject, eventdata, handles)
% hObject    handle to menuPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuManual_Callback(hObject, eventdata, handles)
% hObject    handle to menuManual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuSort_Callback(hObject, eventdata, handles)
% hObject    handle to menuSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuAuto_Callback(hObject, eventdata, handles)
% hObject    handle to menuAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuSpikesort_Callback(hObject, eventdata, handles)
% hObject    handle to menuSpikesort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuLoadPrm_Callback(hObject, eventdata, handles)
% hObject    handle to menuLoadPrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_prm_(hObject);


function vcFile_prm = load_prm_(hObject, vcFile_prm)
if nargin<2, vcFile_prm = ''; end

if isempty(vcFile_prm)
    [vcFile_prm, PathName, FilterIndex] = uigetfile('.prm');
    if FilterIndex==0, return; end
    vcFile_prm = fullfile(PathName, vcFile_prm);
end
S_gui = guidata(hObject);
set(S_gui.Fig_gui, 'Name', vcFile_prm);
enable_(S_gui, {'EditPrm', 'Spikesort', 'Probe', 'Preview', 'Traces', 'Detect', 'ClearPrm'}, 1);
if ~isempty(dir(strrep(vcFile_prm, '.prm', '_jrc.mat')))    
    enable_(S_gui, {'Manual', 'Describe', 'Drift', 'Sort', 'Auto'}, 1);
end
edit(vcFile_prm);
S_gui.vcFile_prm = vcFile_prm;
guidata(hObject, S_gui)


% --------------------------------------------------------------------
function menuMakePrm_Callback(hObject, eventdata, handles)
% hObject    handle to menuMakePrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuEditPrm_Callback(hObject, eventdata, handles)
% hObject    handle to menuEditPrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_file_openbatch_Callback(hObject, eventdata, handles)
% hObject    handle to gui_file_openbatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_view_validate_Callback(hObject, eventdata, handles)
% hObject    handle to gui_view_validate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_file_makebatch_Callback(hObject, eventdata, handles)
% hObject    handle to gui_file_makebatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_help_compile_Callback(hObject, eventdata, handles)
% hObject    handle to gui_help_compile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnLoadPrm.
function btnLoadPrm_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadPrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_prm_(hObject);


% --- Executes on button press in btnDetect.
function btnDetect_Callback(hObject, eventdata, handles)
% hObject    handle to btnDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('detect', S_gui.vcFile_prm);
    enable_(S_gui, {'Describe', 'Sort'}, 1);
    enable_(S_gui, {'Drift', 'Auto'}, 0);
catch
    disperr_();
end


% --- Executes on button press in btnSpikesort.
function btnSpikesort_Callback(hObject, eventdata, handles)
% hObject    handle to btnSpikesort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('spikesort', S_gui.vcFile_prm);
    enable_(S_gui, {'Drift', 'Describe', 'Sort', 'Auto', 'Manual'}, 1);
catch
    disperr_();
end


% --- Executes on button press in btnSort.
function btnSort_Callback(hObject, eventdata, handles)
% hObject    handle to btnSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('sort', S_gui.vcFile_prm);
    enable_(S_gui, {'Drift', 'Describe', 'Auto', 'Manual'}, 1);
catch
    disperr_();
end


% --- Executes on button press in btnWiki.
function btnWiki_Callback(hObject, eventdata, handles)
% hObject    handle to btnWiki (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    jrc('wiki');
catch
    disperr_();
end


% --- Executes on button press in btnEditPrm.
function btnEditPrm_Callback(hObject, eventdata, handles)
% hObject    handle to btnEditPrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    edit(S_gui.vcFile_prm);
catch
    ;
end


% --- Executes on button press in btnTraces.
function btnTraces_Callback(hObject, eventdata, handles)
% hObject    handle to btnTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('traces', S_gui.vcFile_prm);
catch
    disperr_();
end


% --------------------------------------------------------------------
function menuProbe_Callback(hObject, eventdata, handles)
% hObject    handle to menuProbe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnPreview.
function btnPreview_Callback(hObject, eventdata, handles)
% hObject    handle to btnPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('preview', S_gui.vcFile_prm);
catch
    disperr_();
end


% --- Executes on button press in btnManual.
function btnManual_Callback(hObject, eventdata, handles)
% hObject    handle to btnManual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('manual', S_gui.vcFile_prm);
catch
    disperr_();
end


% --- Executes on button press in btnAuto.
function btnAuto_Callback(hObject, eventdata, handles)
% hObject    handle to btnAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('auto', S_gui.vcFile_prm);
    enable_(S_gui, {'Drift', 'Describe', 'Manual'}, 1);
catch
    disperr_();
end


% --- Executes on button press in btnClear.
function btnClear_Callback(hObject, eventdata, handles)
% hObject    handle to btnClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
jrc('clear');

% --- Executes on button press in btnProbe.
function btnProbe_Callback(hObject, eventdata, handles)
% hObject    handle to btnProbe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('probe', S_gui.vcFile_prm);
catch
    disperr_();
end


% --- Executes on button press in btnAbout.
function btnAbout_Callback(hObject, eventdata, handles)
% hObject    handle to btnAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    S_out = jrc('test', 'about_', {}, 1, 0);
    msgbox(S_out.out1);
    fprintf('%s\n', S_out.out1{:});
catch
    disperr_();
end

% --- Executes on button press in btnIssue.
function btnIssue_Callback(hObject, eventdata, handles)
% hObject    handle to btnIssue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
jrc('issue');


%--------------------------------------------------------------------------
% function disperr_()
% Display error message and the error stack
function disperr_(vcMsg)
% ask user to email jrclust@vidriotech.com ? for the error ticket?
dbstack('-completenames'); % display an error stack
vcErr = lasterr();
if nargin==0
    fprintf(2, '%s\n', vcErr);
else
    fprintf(2, '%s:\n\t%s\n', vcMsg, vcErr);
end


% --- Executes on button press in btnDescribe.
function btnDescribe_Callback(hObject, eventdata, handles)
% hObject    handle to btnDescribe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    S_out = jrc('test', 'describe_', {}, 1, 0);
    msgbox(S_out.out1);
    fprintf('%s\n', S_out.out1{:});
catch
    disperr_();
end

% --------------------------------------------------------------------
function menuClearPrm_Callback(hObject, eventdata, handles)
% hObject    handle to menuClearPrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuClear_Callback(hObject, eventdata, handles)
% hObject    handle to menuClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnDrift.
function btnDrift_Callback(hObject, eventdata, handles)
% hObject    handle to btnDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('drift', S_gui.vcFile_prm);
catch
    disperr_();
end

% --------------------------------------------------------------------
function menuDrift_Callback(hObject, eventdata, handles)
% hObject    handle to menuDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnClearPrm.
function btnClearPrm_Callback(hObject, eventdata, handles)
% hObject    handle to btnClearPrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S_gui = guidata(hObject);
try
    jrc('clear', S_gui.vcFile_prm);
    enable_(S_gui, {'Sort', 'Auto', 'Manual', 'Describe', 'Drift'}, 0);
catch
    disperr_();
end

% --- Executes on button press in btnMakePrm.
function btnMakePrm_Callback(hObject, eventdata, handles)
% hObject    handle to btnMakePrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
csAns = inputdlg({'raw recording file', 'Probe file', 'Template file'}, 'Recording format', 1, {'', '', 'default.prm'});
if isempty(csAns), return; end
% S = struct('vcFile', csAns{1}, 'probeFile', csAns{2}, 'template_file', csAns{3});
vcFile_prm = jrc('makeprm', csAns{1}, csAns{2}, csAns{3});
% guidata_set_(hObject, vcFile_prm);
load_prm_(hObject, vcFile_prm);


function S_gui = guidata_set_(hObject, varargin)
S_gui = guidata(hObject);
for i=1:numel(varargin)
    S_gui.(inputname(i+1)) = varargin{i};
end
guidata(hObject, S_gui);


function varargout = guidata_get_(hObject, varargin)
S_gui = guidata(hObject);
for i=1:numel(varargin)
    if isfield(S_gui, varargin{i})
        varargout{i} = S_gui.(varargin{i});
    end
end


% --------------------------------------------------------------------
function menuDetect_Callback(hObject, eventdata, handles)
% hObject    handle to menuDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function enable_(S, csStr, fEnable)
% Enable handle object
if fEnable, vcEnable = 'on'; else vcEnable='off'; end
if ischar(csStr), csStr = {csStr}; end
for i=1:numel(csStr)
    try
        eval(sprintf('set(S.btn%s, ''Enable'', ''%s'');', csStr{i}, vcEnable));
        eval(sprintf('set(S.menu%s, ''Enable'', ''%s'');', csStr{i}, vcEnable));
    catch
        disp(lasterr());
    end
end


% --------------------------------------------------------------------
function menuDescribe_Callback(hObject, eventdata, handles)
% hObject    handle to menuDescribe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
