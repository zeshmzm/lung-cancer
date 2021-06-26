function varargout = contour_fig(varargin)
% CONTOUR_FIG MATLAB code for contour_fig.fig
%      CONTOUR_FIG, by itself, creates a new CONTOUR_FIG or raises the existing
%      singleton*.
%
%      H = CONTOUR_FIG returns the handle to a new CONTOUR_FIG or the handle to
%      the existing singleton*.
%
%      CONTOUR_FIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTOUR_FIG.M with the given input arguments.
%
%      CONTOUR_FIG('Property','Value',...) creates a new CONTOUR_FIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before contour_fig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to contour_fig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help contour_fig

% Last Modified by GUIDE v2.5 12-May-2019 17:26:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @contour_fig_OpeningFcn, ...
                   'gui_OutputFcn',  @contour_fig_OutputFcn, ...
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


% --- Executes just before contour_fig is made visible.
function contour_fig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to contour_fig (see VARARGIN)

% Choose default command line output for contour_fig
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes contour_fig wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = contour_fig_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run('D:\Study\Final Year Project\FYP 3rd rb active contour\lung_cancer_contour.m')