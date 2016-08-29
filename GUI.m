function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 29-Aug-2016 06:15:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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

% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Load preferences
load preferences.mat % This loads a matrix pref
handles.pref = pref; 
%runDTButton_Callback(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function dtEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to dtEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function NEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to NEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function tfEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to tfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function LEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to LEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function aEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to aEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in LRadioButton.
function LRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to LRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LRadioButton
set(handles.LEdit, 'Enable', 'on')
set(handles.aLMultEdit, 'Enable', 'off')


% --- Executes on button press in aRadioButton.
function aRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to aRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aRadioButton
set(handles.LEdit, 'Enable', 'off')
set(handles.aLMultEdit, 'Enable', 'on')

% --- Executes during object creation, after setting all properties.
function potEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to potEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function parabolicEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to parabolicEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cubicRadioButton.
function cubicRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to cubicRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cubicRadioButton
set(handles.potEdit, 'Enable', 'on')
set(handles.parabolicEdit, 'Enable', 'off')
set(handles.vOtherEdit, 'Enable', 'off')

% --- Executes on button press in parabolicRadioButton.
function parabolicRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to parabolicRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parabolicRadioButton
set(handles.potEdit, 'Enable', 'off')
set(handles.vOtherEdit, 'Enable', 'off')
set(handles.parabolicEdit, 'Enable', 'on')

% --- Executes during object creation, after setting all properties.
function vOtherEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to vOtherEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in vOtherRadioButton.
function vOtherRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to vOtherRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vOtherRadioButton
set(handles.potEdit, 'Enable', 'off')
set(handles.vOtherEdit, 'Enable', 'on')
set(handles.parabolicEdit, 'Enable', 'off')

% --- Executes during object creation, after setting all properties.
function iLabel1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to handles.iLabel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function iParamEdit1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to handles.iParamEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function psi0Listbox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to psi0Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% CALLBACK OF INITIAL CONDITION LISTBOX. 
% This callback handles hiding or enabling entry boxes depending on
% what type of IC is enabled.
function psi0Listbox_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to psi0Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns psi0Listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from psi0Listbox
selection=hObject.Value;

if (selection ~= 1)
    handles.xjLabel.Visible = 'off';
    handles.tjLabel.Visible = 'off';
    handles.xjEdit.Visible = 'off';
    handles.tjEdit.Visible = 'off';
    handles.intensityEdit.Visible = 'off';
    handles.intensityCheck.Visible = 'off';
    handles.iCheck1.Visible = 'off';
    handles.decimalsEdit.Visible = 'off';
    handles.ratioEdit.Visible = 'off';
    handles.ratioLabel.Visible = 'off';
    handles.ratioInfo.Visible = 'off';
end
if (selection ~= 2)
    handles.iParam1.Visible = 'on';
    handles.iParamEdit1.Visible = 'on';
    handles.xjEdit.Visible = 'off';
    handles.tjEdit.Visible = 'off';
    handles.intensityEdit.Visible = 'off';
    handles.intensityCheck.Visible = 'off';
end


if (selection == 1)
    handles.iLabel1.String = 'psi0=A0+\Sum_j 2*Aj*cos(j*w*x)';
    handles.iLabel2.Visible = 'off';
    handles.iParam1.String = 'Order';
    handles.iParamEdit1.String = '2';
    handles.iParam2.Visible = 'on';
    handles.iParamEdit2.Visible = 'on';
    handles.iParam2.String = 'T';
    handles.iParamEdit2.String = '-10';
    handles.intensityEdit.Visible = 'on';
    handles.intensityCheck.Visible = 'on';
    handles.iCheck1.Visible = 'on';
    handles.decimalsEdit.Visible = 'on';
    handles.ratioEdit.Visible = 'on';
    handles.ratioLabel.Visible = 'on';
    handles.ratioInfo.Visible = 'on';
    handles.xjLabel.Visible = 'on';
    handles.tjLabel.Visible = 'on';
    handles.xjEdit.Visible = 'on';
    handles.tjEdit.Visible = 'on';
    handles.aEdit.Visible = 'on';
    handles.aLabel.Visible = 'on';
    handles.gCNLabel.Visible = 'on';
    handles.gCNEdit.Visible = 'on';
    handles.seedLabel.Visible = 'on';
    handles.seedMenu.Visible = 'on';

elseif (selection == 2)
    handles.iLabel1.String = 'psi0=A0+2*A1*cos(w*x)';
    handles.iParam1.Visible = 'off';
    handles.iParamEdit1.Visible = 'off';
    handles.iParam2.Visible = 'off';
    handles.iParamEdit2.Visible = 'off';
    handles.iCheck1.Visible = 'off';
    handles.decimalsEdit.Visible = 'off';
    handles.aEdit.Visible = 'on';
    handles.aLabel.Visible = 'on';
    handles.gCNLabel.Visible = 'off';
    handles.gCNEdit.Visible = 'off';
    handles.seedLabel.Visible = 'off';
    handles.seedMenu.Visible = 'off';
elseif (selection == 3)
    handles.iLabel1.String = 'psi0=exp(-x^2/2/s^2)H_n(x/s)';
    handles.iLabel2.String = 'H_n = nth order Hermite polynomial';
    handles.iLabel3.Visible = 'off';
    handles.iParam1.String = 's';
    handles.iParam2.String = 'n';
    handles.iParamEdit1.String = '5';
    handles.iParamEdit2.String = '0';
    handles.iParam2.Visible = 'on';
    handles.iParamEdit2.Visible = 'on';
    handles.iCheck1.Visible = 'off';
    handles.decimalsEdit.Visible = 'off';
    handles.iParam1.Visible = 'on';
    handles.iParamEdit1.Visible = 'on';
    handles.aEdit.Visible = 'off';
    handles.aLabel.Visible = 'off';
    handles.gCNLabel.Visible = 'off';
    handles.gCNEdit.Visible = 'off';
    handles.seedLabel.Visible = 'off';
    handles.seedMenu.Visible = 'off';
elseif (selection == 4)
    handles.iLabel1.String = 'psi0=exp(-x^2/2/s^2)J_n(x)';
    handles.iLabel2.String = 'J_n = nth order Bessel polynomial';
    handles.iLabel3.Visible = 'off';
    handles.iParam1.String = 's';
    handles.iParam2.String = 'n';
    handles.iParamEdit1.String = '10';
    handles.iParamEdit2.String = '0';
    handles.iParam2.Visible = 'on';
    handles.iParamEdit2.Visible = 'on';
    handles.iCheck1.Visible = 'off';
    handles.decimalsEdit.Visible = 'off';
    handles.iParam1.Visible = 'on';
    handles.iParamEdit1.Visible = 'on';
    handles.aEdit.Visible = 'off';
    handles.aLabel.Visible = 'off';
    handles.gCNLabel.Visible = 'off';
    handles.gCNEdit.Visible = 'off';
    handles.seedLabel.Visible = 'off';
    handles.seedMenu.Visible = 'off';
elseif (selection == 5) 
    handles.iLabel1.String = 'psi0=exp(a*x)*Ai(x)';
    handles.iLabel2.String = 'Ai(x) = Airy function';
    handles.iLabel3.Visible = 'off';
    handles.iParam1.String = 'a';
    handles.iParam2.String = 'n';
    handles.iParamEdit1.String = '0.1';
    handles.iParamEdit2.String = '0';
    handles.iParam2.Visible = 'off';
    handles.iParamEdit2.Visible = 'off';
    handles.iCheck1.Visible = 'off';
    handles.decimalsEdit.Visible = 'off';
    handles.iParam1.Visible = 'on';
    handles.iParamEdit1.Visible = 'on';
    handles.aEdit.Visible = 'off';
    handles.aLabel.Visible = 'off';
    handles.gCNLabel.Visible = 'off';
    handles.gCNEdit.Visible = 'off';
    handles.seedLabel.Visible = 'off';
    handles.seedMenu.Visible = 'off';
end

% --- Executes during object creation, after setting all properties.
function iParamEdit2_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to handles.iParamEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function orderBox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to orderBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function aLMultEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to aLMultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function aLMultEdit_Callback(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to aLMultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in iCheck1.
function iCheck1_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to iCheck1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value == 1
    handles.decimalsEdit.Enable = 'on';
else
    handles.decimalsEdit.Enable = 'off';
end
 
% --- Executes on button press in analSwitch.
function analSwitch_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to analSwitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pref = handles.pref;
PSI = handles.PSI_anal;
x = handles.x_anal;
t = handles.t_anal;
if strcmp(handles.analPlotMode, 'density'); % Axes are in density plot mode
    nlsePlot(abs(PSI).^pref.pwr, x, t, 1, 1, handles.analRealAxes, '3D');
    handles.analPlotMode = '3D';
else
    nlsePlot(abs(PSI).^pref.pwr, x, t, 1, 1, handles.analRealAxes, 'density');
    handles.analPlotMode = 'density';
end
title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^pref.pwr))));
guidata(hObject, handles);

% --- Executes on button press in numSwitch.
function numSwitch_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to numSwitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pref = handles.pref;
PSI = handles.PSI_num;
x = handles.x_num;
t = handles.t_num;
Nx = length(x);
Nt = length(t);
if strcmp(handles.numPlotMode, 'density'); % Axes are in density plot mode
    nlsePlot(abs(PSI).^pref.pwr, x, t, ceil(Nt/1000), ceil(Nx/256), handles.numRealAxes, '3D');
    handles.numPlotMode = '3D';
else
    nlsePlot(abs(PSI).^pref.pwr, x, t, ceil(Nt/1000), ceil(Nx/256), handles.numRealAxes, 'density');
    handles.numPlotMode = 'density';
end
title(sprintf('Numerical. Max = %.3f', max(max(abs(PSI).^pref.pwr))));
guidata(hObject, handles);

% --- Executes on button press in undockDT1.
function undockDT1_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to undockDT1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
h = allchild(handles.analRealAxes);
copyobj(handles.analRealAxes,f); colormap('jet'); colorbar;
%title(' ');
if (size(h) == 1); % Density Plot mode
    colorbar
    set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);
else
    colorbar off;
    set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);
end

% --- Executes on button press in undockDT2.
function undockDT2_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to undockDT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
leg = legend(handles.analInverseAxes);
if ~isempty(leg)
    copyobj([leg, handles.analInverseAxes],f);
else
    copyobj(handles.analInverseAxes,f);
    colorbar; colormap('jet');
end
set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);

% --- Executes on button press in undockNum1.
function undockNum1_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to undockNum1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
h = allchild(handles.numRealAxes);
copyobj(handles.numRealAxes,f); colormap('jet'); colorbar;
%title(' ');
if (size(h) == 1); % Density Plot mode
    colorbar
    set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);
else
    colorbar off;
    set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);
end

% --- Executes on button press in undockNum2.
function undockNum2_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to undockNum2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
leg = legend(handles.numInverseAxes);
if ~isempty(leg)
    copyobj([leg, handles.numInverseAxes],f);
else
    copyobj(handles.analInverseAxes,f);
    colorbar; colormap('jet');
end
set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);

% --- Executes during object creation, after setting all properties.
function xjEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to xjEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function tjEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to tjEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function decimalsEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to decimalsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(6, 1));
set(hObject, 'RowName', {'A0', 'A1', 'A2', 'A3', 'A4', 'A5'});

% --- Executes on button press in pubSizeButton.
function pubSizeButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to pubSizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = handles.numRealAxes;
axes(ax);
maxY = ax.YLim(2);
ylim([maxY/2-2, maxY/2+2]);
maxX = ax.XLim(2);
if maxX > 4;
    xlim([-4, 4]);
end
title('');

% --------------------------------------------------------------------
function menuPreferences_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to menuPreferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Preferences;
load preferences.mat;
handles.pref = pref;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function intensityEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to intensityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function potMenu_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to potMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns potMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from potMenu
str = get(handles.potMenu, 'String');
val = get(handles.potMenu,'Value');
% Set current data to the selected data set.
switch str{val};
case 'g*|psi|^2 (cubic)' % User selects peaks.
    handles.potEdit.String = num2str(-1);
case 'alpha*x^2 (parabolic)' % User selects membrane.
    handles.potEdit.String = num2str(0.3);
end

% --- Executes during object creation, after setting all properties.
function potMenu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to potMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = handles.numInverseAxes;
axes(ax);
h = allchild(ax);
if handles.SpectrumPubMode
CO =     [         0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
    handles.SpectrumPubMode = 0;
    title('Numerical');
    width = 1;
    maxX = ax.XLim(2);
    minX = ax.XLim(1);
    xlim([0, minX + maxX]);
    ylim ([-10, 2]);
    ax.XTick = 0:(minX+maxX)/4:minX+maxX;
    ax.YTick = -10:2:2;
else
CO=[0.0314    0.1882    0.4196; 
    0         0.5       0.1412;  
    0.8941    0.1020    0.1098; 
    1.0000    0.4980    0     ; 
    0.3059    0.7020    0.8275;
    0.6824    0.0039    0.4941;
    0.4980         0         0];
    handles.SpectrumPubMode = 1;
    title('');
    width = 1.5;
    ylim([-8, 2]);
    maxX = ax.XLim(2);
    xlim([maxX/2-2, maxX/2+2]);
    ax.XTick = maxX/2-2:1:maxX/2+2;
    ax.YTick = -8:2:2;
end

nLines = length(h);
for i = 1:nLines
	h(i).Color = CO(nLines-i+1, :);
    h(i).LineWidth = width;
end

guidata(hObject,handles)

% --------------------------------------------------------------------
function vidSettings_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to vidSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
videoSettings

% --- Executes on button press in numericalMaxima.
function numericalMaxima_Callback(hObject, ~, handles)
% hObject    handle to numericalMaxima (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maximaType = 1;     % This indicates function should use numerical data to compute maxima
guidata(hObject, handles);
maximaDisplay;

% --- Executes on button press in analyticalMaxima.
function analyticalMaxima_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to analyticalMaxima (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maximaType = 2;         % This indicates function should use DT data to compute maxima
guidata(hObject, handles);
maximaDisplay;

% --------------------------------------------------------------------
function menuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to menuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function viewMenu_Callback(hObject, eventdata, handles)
% hObject    handle to viewMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function analParamView_Callback(hObject, eventdata, handles)
% hObject    handle to analParamView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
analProp;


% --- Executes during object creation, after setting all properties.
function ratioEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ratioEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in numEnergy.
function numEnergy_Callback(hObject, eventdata, handles)
% hObject    handle to numEnergy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g = str2double(handles.potEdit.String);
Nx = str2double(handles.NEdit.String);
dt = str2double(handles.dtEdit.String);
[dE, E, ke, pe] = energy(handles.PSI_num, handles.t_num, handles.k_num.^2, Nx, g, dt);
figure;
plot(handles.t_num, dE, 'LineWidth', 1.5); title('dE'); title('Integrated Energy Error');
xlabel('t'); ylabel('dE'); grid on;

figure;
plot(handles.t_num, E, 'LineWidth', 1.5); title('Energy');
hold on;
plot(handles.t_num, ke, 'LineWidth', 1.5);  
hold on;
plot(handles.t_num, pe, 'LineWidth', 1.5);
legend('E', 'T', 'V', 'Location', 'Best');
xlabel('t'); ylabel('E'); grid on;


% --------------------------------------------------------------------
function pubModeSettings_Callback(hObject, eventdata, handles)
% hObject    handle to pubModeSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pubModeSettings;

% --- Executes on selection change in seedMenu.
function seedMenu_Callback(hObject, eventdata, handles)
% hObject    handle to seedMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = hObject.String;
val = hObject.Value;
seed = str{val};
switch lower(seed)
    case {'breather', 'soliton'}
        handles.gCNEdit.Enable = 'off';
    case {'cn', 'dn'}
        handles.gCNEdit.Enable = 'on';
    otherwise
        error('Unknown Seed.');
end

% --- Executes during object creation, after setting all properties.
function seedMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seedMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function gCNEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gCNEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runDTButton.
function runDTButton_Callback(hObject, eventdata, handles)
% hObject    handle to runDTButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
order =  str2double(handles.iParamEdit1.String);
a = eval(handles.aEdit.String); 
L = pi/sqrt(1-2*a(1));                        %#ok<NASGU> % Read length L (will be used in shifts)    
lambda = sqrt(8*a(1)*(1-2*a(1)));             %#ok<NASGU> % Growth factor (used in shifts)
Omega = 2*sqrt(1-2*a(1));                     %#ok<NASGU> % Fundamental wavenumber (used in shifts)  
xj = eval(handles.xjEdit.String);
tj = eval(handles.tjEdit.String);
T = str2double(handles.iParamEdit2.String);
R = eval(handles.ratioEdit.String);
g = eval(handles.gCNEdit.String);
pref = handles.pref;
str = get(handles.seedMenu, 'String');
val = get(handles.seedMenu,'Value');
seed = str{val};
[PSI, x, t] = calcDarboux(order, a, R, T, g, seed, xj, tj, pref.Nx, pref.Nt, pref.Lmode, pref.L, pref.aLMult);  % FIX PLEASE

peak = peakPredict(a, order, R, seed, g);
handles.peakLabel.String = num2str(abs(peak).^handles.pref.pwr);

nlsePlot(abs(PSI).^pref.pwr, x, t, 1, 1, handles.analRealAxes, pref.plotType);
%k = 2*(-pref.Nx/2:1:pref.Nx/2-1)'*pi/2/max(x);          % Wave number
k = (-pref.Nx/2:pref.Nx/2-1);
fourierPlot(PSI, t, k, pref.fourierLines, 1, handles.analInverseAxes, pref.fourierPlotType);

axes(handles.analRealAxes); title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^pref.pwr))));
axes(handles.analInverseAxes); title('Analytical');

handles.analPlotMode = pref.plotType;
handles.analFourierPlotMode = pref.fourierPlotType;
handles.PSI_anal = PSI;                                 % Save result to handles struct
handles.x_anal = x;
handles.t_anal = t;
handles.order = order;
handles.k_anal = k;
guidata(hObject, handles);

% --- Executes on key press with focus on main_GUI and none of its controls.
function main_GUI_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to main_GUI (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key
  case 'return'
    runDTButton_Callback(hObject, eventdata, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN BUTTON CALL BACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in runButton.
function runButton_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    han    5dle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    dt = str2double(handles.dtEdit.String);                    % Read dt from edit
    Nx =  str2double(handles.NEdit.String);                    % read Nx from edit
    Tmax = str2double(handles.tfEdit.String);                  % Read Tmax from edit
    Nt = Tmax/dt;                                              % Calculate Nt
    mult = str2double(handles.aLMultEdit.String);              % Read mult from edit box
    R = eval(handles.ratioEdit.String);
    a = eval(handles.aEdit.String); 
    
    if strcmp(handles.LEdit.Enable, 'on')
        L = str2double(handles.LEdit.String);
    else
        [~, D] = rat(R);
        if ~imag(a) % If parameter 'a' is used.
            aa = a(1); 
        else % Otherwise, complex eigenvalue form is assumed
            aa = imag(a(1))^2/2; % Calculate parameter 'a' from 'l'
        end  
        disp(aa);
        L = D*mult*pi/sqrt(1-2*aa);     % Periodic length
        Omega = 2*sqrt(1-2*aa);      % Principal wave number
    end
        
    xj = eval(handles.xjEdit.String);
    tj = eval(handles.tjEdit.String);
    
    dx = L/Nx;                             % Spatial step size
    x = (-Nx/2:1:Nx/2-1)'*dx;
    
    % Determine the selected data set.
    str = get(handles.potMenu, 'String');
    val = get(handles.potMenu,'Value');
    % Set current data to the selected data set.
    switch str{val};
    case 'g*|psi|^2 (cubic)' % Cubic NLSE
        g = str2double(handles.potEdit.String);
        V = @(psi, x) (g*abs(psi).^2);
    case 'alpha*x^2 (parabolic)' % Parabolic SE
        alpha = str2double(handles.potEdit.String);
        V = @(psi, x) 1/2*alpha*x.^2;
    end
    
    psi0_selection = handles.psi0Listbox.Value;
    
    if (psi0_selection == 1) % THIS SHOULD BE MOVED TO GEN COEFF!!!
        order =  str2double(handles.iParamEdit1.String);
        T = str2double(handles.iParamEdit2.String);
        R = eval(handles.ratioEdit.String);
        g = eval(handles.gCNEdit.String);
        str = get(handles.seedMenu, 'String');
        val = get(handles.seedMenu,'Value');
        seed = str{val};
        [psi_dt, x_dt, ~] = calcDarboux(order, a, R, T, g, seed, xj, tj, Nx, 0, 'periodic', 0, mult);
        if length(R) == 1 && R == 2
            R = 1:order;
        else
            R = [1, R];
        end
        [A0, A] = genCoeff(x_dt, psi_dt, a, order, R, 1, handles.uipanel2, handles.iCheck1.Value);
             
        psi_0 = A0;
        for i=1:order
            psi_0 = psi_0 + 2*A(i)*cos(R(i)*Omega*x);
        end
%         norm = dx*(psi_0')*(psi_0)/L;
%         disp(['Norm is: ', num2str(norm)]);
        handles.uitable1.Data = [A0, A].';
        
    elseif (psi0_selection == 2)
        
        A = handles.uitable1.Data;
        if iscell(A)
            for i = 1:length(A)
                if isempty(A{i})
                    B(i) = 0;
                else
                    B(i) = A{i};
                end
            end
            A = B.';
        end
        A = A(2:end);
        
        A02 = 1;
        for i = 1:length(A)
            A02 = A02 - 2*abs(A(i))^2;
        end
        A0 = sqrt(A02);
        
        psi_0 = A0;
        for i=1:length(A)
            psi_0 = psi_0 + 2*A(i)*cos(R(i)*Omega*x);
        end
        handles.uitable1.Data = [A0; A];
    elseif (psi0_selection == 3)
        s = str2double(handles.iParamEdit1.String);
        n = str2double(handles.iParamEdit2.String);
        psi_0 = exp(-x.^2/2/s^2).*hermiteH(n, x/s);
    elseif (psi0_selection == 4)
        s = str2double(handles.iParamEdit1.String);
        n = str2double(handles.iParamEdit2.String);
        psi_0 = exp(-x.^2/2/s^2).*besselj(n, x);
    elseif (psi0_selection == 5)
        a = str2double(handles.iParamEdit1.String);
        psi_0 = exp(a*x).*airy(x);
    end
    
    orderSelection = handles.orderBox.Value;
    if orderSelection == 1
        method = 'T1';
    elseif orderSelection == 2
        method = 'T2';
    elseif orderSelection == 3
        method = 'T4S';
    elseif orderSelection == 4
        method = 'T4M';
    elseif orderSelection == 5
        method = 'T6S';
    elseif orderSelection == 6
        method = 'TM';
    elseif orderSelection == 7
        method = 'T8S';    
    elseif orderSelection == 8
        method = 'T8M';
    end
    [PSI, x, t, k] = solve(dt, Nx, Tmax, L, mult, V, psi_0, method); 
    handles.k_num = k;

%    initialPlot(PSI(1, :), x, handles.axes1);
    pref = handles.pref;
    nlsePlot(abs(PSI).^pref.pwr, x, t, ceil(Nt/1000), ceil(Nx/256), handles.numRealAxes, pref.plotType);
    k = (-Nx/2:Nx/2-1);
    fourierPlot(PSI, t, k, pref.fourierLines, 1, handles.numInverseAxes, pref.fourierPlotType);
  
    axes(handles.numRealAxes); title(sprintf('Numerical. Max = %.3f', max(max(abs(PSI).^pref.pwr))));
    axes(handles.numInverseAxes); title('Numerical'); 
    handles.PSI_num = PSI;           % Save result to handles structure.
    handles.x_num = x;
    handles.t_num = t;
    handles.numPlotMode = pref.plotType;
    handles.numFourierPlotMode = pref.fourierPlotType;
    
    handles.SpectrumPubMode = 0;
    handles.IntensityPubMode = 0;
    guidata(hObject,handles)


% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function analExportMenu_Callback(hObject, eventdata, handles)
k = any(strcmp('PSI_anal',fieldnames(handles))); % Check if numerical data exists in memory
if ~k % If it doesn't, display error
    errordlg('No DT data exists in memory, nothing to export.');
end
PSI = handles.PSI_anal;
x = handles.x_anal;
t = handles.t_anal;
exportData(PSI, x, t, handles);

% --------------------------------------------------------------------
function numExportMenu_Callback(hObject, eventdata, handles)
% hObject    handle to numExportMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load data to be saved
k = any(strcmp('PSI_num',fieldnames(handles))); % Check if numerical data exists in memory
if ~k % If it doesn't, display error
    errordlg('No numerical data exists in memory, nothing to export.');
end
PSI = handles.PSI_num;
x = handles.x_num;
t = handles.t_num;
exportData(PSI, x, t, handles);

% --------------------------------------------------------------------
function exportData(PSI, x, t, handles)
% hObject    handle to exportDataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

defaultName = [datestr(date, 'yyyymmdd'), '.mat'];

extensions = {'*.mat',...
 'MATLAB MAT file (*.mat)';
 '*.csv', 'Comma Separate Values (*.csv)';...
 '*.txt','Space-delimited TXT (*.txt)';...
 '*.*',  'All Files (*.*)'};

[filename, pathname, index] = uiputfile(extensions, 'Save as', defaultName);

if isequal(filename,0) || isequal(pathname,0) % Canceled
   return;
else
   fullfile = [pathname, '\', filename];
   if index == 1 % mat
       save(fullfile, 'PSI', 'x', 't');
   elseif index == 2 %csv
       full = [[NaN, x']; t PSI];
       dlmwrite(fullfile,full,'delimiter',',','precision',handles.pref.dec)
   elseif index == 3 % space-delimited txt
       full = [[NaN, x']; t PSI];
       dlmwrite(fullfile,full,'delimiter',' ','precision',handles.pref.dec)
   else
       error('Cant save in this extensions');
   end
end


% --------------------------------------------------------------------
function exportDataMenu_Callback(hObject, eventdata, handles)
% hObject    handle to exportDataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on main_GUI or any of its controls.
function main_GUI_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to main_GUI (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key
  case 'return'
    runDTButton_Callback(hObject, eventdata, handles);
end


function analFourierSwitch_Callback(hObject, eventdata, handles)
pref = handles.pref;
PSI = handles.PSI_anal;
t = handles.t_anal;
k = handles.k_anal;
if strcmp(handles.analFourierPlotMode, 'lines'); % Axes are in density plot mode
    fourierPlot(PSI, t, k, pref.fourierLines, 1, handles.analInverseAxes, 'density');
    handles.analFourierPlotMode = 'density';
else
    fourierPlot(PSI, t, k, pref.fourierLines, 1, handles.analInverseAxes, 'lines');
    handles.analFourierPlotMode = 'lines';
end
title('Analytical');
guidata(hObject, handles);


function numFourierSwitch_Callback(hObject, eventdata, handles)
pref = handles.pref;
PSI = handles.PSI_num;
t = handles.t_num;
k = handles.k_num;
if strcmp(handles.numFourierPlotMode, 'lines'); % Axes are in density plot mode
    fourierPlot(PSI, t, k, pref.fourierLines, 1, handles.numInverseAxes, 'density');
    handles.numFourierPlotMode = 'density';
else
    fourierPlot(PSI, t, k, pref.fourierLines, 1, handles.numInverseAxes, 'lines');
    handles.numFourierPlotMode = 'lines';
end
title('Numerical');
guidata(hObject, handles);
