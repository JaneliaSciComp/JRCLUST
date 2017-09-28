function varargout = jrc3_gui(varargin)
% JRC3_GUI MATLAB code for jrc3_gui.fig
%      JRC3_GUI, by itself, creates a new JRC3_GUI or raises the existing
%      singleton*.
%
%      H = JRC3_GUI returns the handle to a new JRC3_GUI or the handle to
%      the existing singleton*.
%
%      JRC3_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JRC3_GUI.M with the given input arguments.
%
%      JRC3_GUI('Property','Value',...) creates a new JRC3_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before jrc3_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to jrc3_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help jrc3_gui

% Last Modified by GUIDE v2.5 28-Sep-2017 07:16:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @jrc3_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @jrc3_gui_OutputFcn, ...
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


% --- Executes just before jrc3_gui is made visible.
function jrc3_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to jrc3_gui (see VARARGIN)
image(handles.axes1, imread('./img/gui_background.png'));
set(handles.axes1, 'XTick', [], 'YTick', [], 'Visible', 'off');
axis(handles.axes1, 'equal');
% Choose default command line output for jrc3_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes jrc3_gui wait for user response (see UIRESUME)
% uiwait(handles.Fig_gui);


% --- Outputs from this function are returned to the command line.
function varargout = jrc3_gui_OutputFcn(hObject, eventdata, handles) 
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
function gui_help_wiki_Callback(hObject, eventdata, handles)
% hObject    handle to gui_help_wiki (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_help_gitpull_Callback(hObject, eventdata, handles)
% hObject    handle to gui_help_gitpull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_help_update_Callback(hObject, eventdata, handles)
% hObject    handle to gui_help_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_help_about_Callback(hObject, eventdata, handles)
% hObject    handle to gui_help_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_help_issue_Callback(hObject, eventdata, handles)
% hObject    handle to gui_help_issue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_view_traces_Callback(hObject, eventdata, handles)
% hObject    handle to gui_view_traces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_view_preview_Callback(hObject, eventdata, handles)
% hObject    handle to gui_view_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_view_manual_Callback(hObject, eventdata, handles)
% hObject    handle to gui_view_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_edit_sort_Callback(hObject, eventdata, handles)
% hObject    handle to gui_edit_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_edit_auto_Callback(hObject, eventdata, handles)
% hObject    handle to gui_edit_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_edit_spikesort_Callback(hObject, eventdata, handles)
% hObject    handle to gui_edit_spikesort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_file_openprm_Callback(hObject, eventdata, handles)
% hObject    handle to gui_file_openprm (see GCBO)
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
vhBtn = [S_gui.btnEditPrm, S_gui.btnSpikesort, S_gui.btnProbe, S_gui.btnPreview, S_gui.btnTraces, S_gui.btnDetect, S_gui.btnClearPrm];
set(vhBtn, 'Enable', 'on');
if ~isempty(dir(strrep(vcFile_prm, '.prm', '_jrc.mat')))
    set([S_gui.btnManual, S_gui.btnDescribe, S_gui.btnDrift, S_gui.btnSort, S_gui.btnAuto], 'Enable', 'on');
end
edit(vcFile_prm);
S_gui.vcFile_prm = vcFile_prm;
guidata(hObject, S_gui)


% --------------------------------------------------------------------
function gui_file_makeprm_Callback(hObject, eventdata, handles)
% hObject    handle to gui_file_makeprm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_edit_prm_Callback(hObject, eventdata, handles)
% hObject    handle to gui_edit_prm (see GCBO)
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
    set([S_gui.btnDescribe, S_gui.btnSort], 'Enable', 'on');
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
    set([S_gui.btnDrift, S_gui.btnDescribe], 'Enable', 'on');
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
    set([S_gui.btnDrift, S_gui.btnDescribe, S_gui.btnAuto], 'Enable', 'on');
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
function gui_view_probe_Callback(hObject, eventdata, handles)
% hObject    handle to gui_view_probe (see GCBO)
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
    set([S_gui.btnDrift, S_gui.btnDescribe], 'Enable', 'on');
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
    S_out = jrc3('test', 'about_', {}, 1, 0);
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
    S_out = jrc3('test', 'describe_', {}, 1, 0);
    msgbox(S_out.out1);
    fprintf('%s\n', S_out.out1{:});
catch
    disperr_();
end

% --------------------------------------------------------------------
function gui_edit_clearprm_Callback(hObject, eventdata, handles)
% hObject    handle to gui_edit_clearprm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gui_edit_clear_Callback(hObject, eventdata, handles)
% hObject    handle to gui_edit_clear (see GCBO)
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
function gui_view_drift_Callback(hObject, eventdata, handles)
% hObject    handle to gui_view_drift (see GCBO)
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
    set([S_gui.btnSort, S_gui.btnAuto, S_gui.btnManual, S_gui.btnDescribe, S_gui.btnDrift], 'Enable', 'off');
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
% S = struct('vcFile', csAns{1}, 'probe_file', csAns{2}, 'template_file', csAns{3});
vcFile_prm = jrc3('makeprm', csAns{1}, csAns{2}, csAns{3});
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
