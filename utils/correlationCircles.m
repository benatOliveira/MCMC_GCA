function h = correlationCircles(data, fig_number,varargin)
%CORRELATIONCIRCLES Represent correlation matrix using colored circles
%
%   correlationCircles(DATA)
%
%   Example
%     % Simple example on iris
%     s = load('fisheriris');
%     correlationCircles(s.meas, 'varNames', {'SepalLength', 'SepalWidth', 'PetalLength', 'PetalWidth'})
%
%     % Another example with more variables
%     load cities
%     correlationCircles(ratings, 'varNames', categories)
%
%   See also
%     corrcoef
%

% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2012-07-16,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2012 INRA - Cepia Software Platform.


%% Parse input arguments

varNames = {};
while length(varargin) >= 2
    argName = varargin{1};
    switch lower(argName)
        case 'varnames'
            varNames = varargin{2};
        otherwise
            error(['Unkown optional argument: ' argName]);
    end
    varargin(1:2) = [];
end

% varNames

%% Initialisations

% number of variables
n = size(data, 2);
nRows = n;
nCols = n;

% compute the correlation matrix of input array
% its values will be used for computing radius and color of circles
cc = corrcoef(data);


%% Prepare display

% initialize the colormap
cmap = jet(256);

% Create the main axis containing all small axes
figure(fig_number);
cax = gca;
axRef = newplot(cax);
fig = ancestor(axRef, 'figure');
holdState = ishold(axRef);
set(axRef, ...
    'Visible', 'off', ...
    'xlim', [-1 1], ...
    'ylim', [-1 1], ...
    'DataAspectRatio', [1 1 1], ...
    'PlotBoxAspectRatio', [1 1 1], ...
    'color', 'none');

pos = get(axRef, 'Position');

% size of sub-plots
width   = pos(3) / (nRows+1);
height  = pos(4) / (nCols+1);

% 2 percent space between axes
space = .02;
pos(1:2) = pos(1:2) + space * [width height];

% grab some info on ref axis
axRefVisibility = get(axRef, 'HandleVisibility');
axRefParent = get(axRef, 'Parent');

% pre-compute data for drawing circles
t = linspace(0, 2*pi, 100)';
cx = cos(t);
cy = sin(t);


%% Plot all circles

% iterate over all cells
for i = nRows:-1:1
    for j = nCols:-1:1
        % compute the position within the main figure
        % (may be necessary to adjust x and y depending on cell number...)
        axPos = [...
            pos(1) + j * width  ...
            pos(2) + (nRows-i) * height ...
            width * (1-space) ...
            height * (1-space)];
        
        % create the axes
        ax(i,j) = axes(...
            'Position', axPos, ...
            'HandleVisibility', axRefVisibility, ...
            'parent', axRefParent);
        
        % color and radius of current correlation circle
        indColor = min(floor((cc(i,j) + 1) * 128) + 1, 256);
        color = cmap(indColor, :);
        r = abs(cc(i, j));
        
        if i>j
        % fill a disc
        hh(i,j) = fill(cx*r, cy*r, color);
        end
        % normalise the display of each axis
        set(ax(i,j), ...
            'xlim', [-1 1], ...
            'ylim', [-1 1], ...
            'DataAspectRatio', [1 1 1], ...
            'xgrid', 'off', ...
            'ygrid', 'off', ...
            'Visible', 'off');
        
        if i<=j
%         text('Parent',ax(i,j),'String',num2str(cc(i, j)),'Position',[i j 0],'Visible','on');
        text('Parent',ax(i,j),'String',num2str(cc(i, j)),'Visible','on');
        end
        
    end
end


%% Display fixup

% remove all labels
set(ax(:), 'xticklabel', '')
set(ax(:), 'yticklabel', '')

% setup tags
set(axRef, ...
    'userdata', ax, ...
    'tag', 'CorrelationCirclesRefAxis');
set(ax, 'tag', 'CorrelationCirclesAxis');

% make axRef the current axis
set(fig, 'CurrentAx', axRef)
if ~holdState,
    set(fig, 'NextPlot', 'replace')
end

% set Title and X/YLabel visible, but use empty labels
textHandles = [get(axRef,'Title'); get(axRef,'XLabel'); get(axRef,'YLabel')]; 
set(textHandles, 'String', '', 'Visible', 'on')

% display labels on the left and on the top of the circle array
if ~isempty(varNames)
    % eventually convert to cell array
    if ischar(varNames)
        varNames = cellstr(varNames);
    end
    
    % display each variable name 
    for i = 1:n
        set(ax(1,i), 'XAxisLocation', 'Top');
        xlabel(ax(1, i), varNames{i}, ...
            'Visible', 'On', ...
            'FontSize', 12, 'Rotation', 45, ...
            'VerticalAlignment', 'Middle', ...
            'HorizontalAlignment', 'Left');
        ylabel(ax(i, 1), varNames{i}, ...
            'Visible', 'On', ...
            'FontSize', 12, 'Rotation', 0, ...
            'VerticalAlignment', 'Middle', ...
            'HorizontalAlignment', 'Right');
    end
end

% Create colorbar
colormap('jet')
colorbar(axRef,'Position',...
    [0.120481927710843 0.103323699421965 0.0113394755492558 0.815028901734101],'XTickLabel',{'-1','-0.8','-0.6','-0.4','-0.2','0','0.2','0.4','0.6','0.8','1'}, ...
               'XTick', 0:0.1:1);

% return handles if needed
if nargout ~= 0
    h = hh;
end


