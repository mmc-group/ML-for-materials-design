clc; clear all; close all;

% Note: this code only allows visualizing up to TWO stiffness at a time.

% To-do: enter one OR two 9-dimensional stiffness vector with single []
% brackets inside {} brackets

S_list = {
    [0.3068, 0.0560, 0.1079, 0.2979, 0.1073, 0.6205, 0.1163, 0.1183, 0.0498],
    [0.2041, 0.0327, 0.0377, 0.0758, 0.0246, 0.1285, 0.0288, 0.0489, 0.0412]
     };

%------------------------------------------------------------------------
% Code: Do not edit anything below this line

visualize_2d(S_list)

function [] =  visualize_2d(S_list)
for s_iter =1:length(S_list)

    S = S_list{s_iter};


    np = 100;
    cols={'-b','-r','-k'};
    global Rmax
    
    fontsize = 24;
    
    contrast = 1;
    Cvec = get_Cvec(S);
    
    
    S = getCompliance(Cvec);
    if s_iter==1
        Rmax = 1.05*max([1/S(1,1,1,1),1/S(2,2,2,2),1/S(3,3,3,3)]);
    end
    
    subplot(1,3,1)
    if s_iter>1
        hold on
    end
    makePlotXY(Cvec,np,Rmax,cols{s_iter});
    hold on
    set(gca,'FontSize',fontsize)
    quiver(0,0,1.1*Rmax,0,'-k','LineWidth',1)
    quiver(0,0,0,1.1*Rmax,'-k','LineWidth',1)
    text(1.1*Rmax,0,'x','HandleVisibility', 'off','FontSize',fontsize);
    text(0,1.1*Rmax,'y','HandleVisibility', 'off','FontSize',fontsize);
    hold off
    
    subplot(1,3,2)
    if s_iter>1
        hold on
    end
    makePlotYZ(Cvec,np,Rmax,cols{s_iter});
    hold on
    set(gca,'FontSize',fontsize)
    quiver(0,0,1.1*Rmax,0,'-k','LineWidth',1)
    quiver(0,0,0,1.1*Rmax,'-k','LineWidth',1)
    text(1.1*Rmax,0,'y','HandleVisibility', 'off','FontSize',fontsize);
    text(0,1.1*Rmax,'z','HandleVisibility', 'off','FontSize',fontsize);
    hold off
    
    subplot(1,3,3)
    if s_iter>1
        hold on
    end
    makePlotXZ(Cvec,np,Rmax,cols{s_iter});
    hold on
    set(gca,'FontSize',fontsize)
    quiver(0,0,1.1*Rmax,0,'-k','LineWidth',1)
    quiver(0,0,0,1.1*Rmax,'-k','LineWidth',1)
    text(1.1*Rmax,0,'x','HandleVisibility', 'off','FontSize',fontsize);
    text(0,1.1*Rmax,'z','HandleVisibility', 'off','FontSize',fontsize);
    hold off
    
end
set(gcf,'Position',[189 198 1298 398])
end

function [Cvec] =  get_Cvec(S)

C = zeros(6,6);
C(1,1)=S(1);
C(1,2)=S(2);
C(1,3)=S(3);

C(2,1) = C(1,2);
C(2,2) = S(4);
C(2,3) = S(5);

C(3,1) = C(1,3);
C(3,2) = C(2,3);
C(3,3) = S(6);

C(4,4) = S(7);
C(5,5) = S(8);
C(6,6) = S(9);

Cvec = reshape(C,6,6)';
end



function [h] = makePlotXY(Cvec,np,Rmax,str)
S = getCompliance(Cvec);
theta=linspace(0,2*pi,np);
E=zeros(1,np);
Ex=zeros(1,np);
Ey=zeros(1,np);
Ez=zeros(1,np);
for i=1:length(theta)
    [Ex(i),Ey(i),Ez(i),E(i)] = getE(S,pi/2,theta(i));
end
% polar(theta,R*ones(size(theta)),':w'); hold on
h=polar(theta,E,str);
end


function [h] = makePlotYZ(Cvec,np,Rmax,str)
S = getCompliance(Cvec);
phi=linspace(0,2*pi,np);
E=zeros(1,np);
Ex=zeros(1,np);
Ey=zeros(1,np);
Ez=zeros(1,np);
for i=1:length(phi)
    [Ex(i),Ey(i),Ez(i),E(i)] = getE(S,phi(i)+pi/2,pi/2);
end
% polar(phi,R*ones(size(phi)),':w'); hold on
h=polar(phi,E,str);
end


function [h] = makePlotXZ(Cvec,np,Rmax,str)
S = getCompliance(Cvec);
phi=linspace(0,2*pi,np);
E=zeros(1,np);
Ex=zeros(1,np);
Ey=zeros(1,np);
Ez=zeros(1,np);
for i=1:length(phi)
    [Ex(i),Ey(i),Ez(i),E(i)] = getE(S,phi(i)+pi/2,0);
end
% polar(phi,R*ones(size(phi)),':w'),str; hold on
h=polar(phi,E,str);
end




function [S] = getCompliance(Cvec)


map = [...
    1,6,5;...
    6,2,4;...
    5,4,3];


Cv = reshape(Cvec,[6,6]);
Cv = 0.5*(Cv+Cv');
Sv = inv(Cv);

S=zeros([3,3,3,3]);

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                vr = map(i,j);
                vc = map(k,l);
                
                factor=1;
                if(vr>3)
                    factor = 2*factor;
                end
                if(vc>3)
                    factor = 2*factor;
                end
                
                S(i,j,k,l) = 1/factor * Sv(vr,vc);
                    
            end
        end
    end
end

end

function [Ex,Ey,Ez,E] = getE_3d(S,phi,theta)
E=zeros(length(phi),length(theta));
Ex=E;
Ey=E;
Ez=E;
for p = 1:length(phi)
   for t = 1:length(theta)
       [Ex(p,t),Ey(p,t),Ez(p,t),E(p,t)] = getE(S,phi(p),theta(t));
   end
end

end


function [x,y,z,E] = getE(S,phi,theta)
d=[cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)];
E = 1/doubleContraction(S,d);
x = E*d(1);
y = E*d(2);
z = E*d(3);
end

function [x] = doubleContraction(S,d)
x=0;
for i=1:3
       for j=1:3
           for k=1:3
               for l=1:3
                   x = x+S(i,j,k,l)*d(i)*d(j)*d(k)*d(l);
               end
           end
       end
end
end

function hpol = polar(varargin)
    %POLAR  Polar coordinate plot.
    %   POLAR(THETA, RHO) makes a plot using polar coordinates of
    %   the angle THETA, in radians, versus the radius RHO.
    %   POLAR(THETA, RHO, S) uses the linestyle specified in string S.
    %   See PLOT for a description of legal linestyles.
    %
    %   POLAR(AX, ...) plots into AX instead of GCA.
    %
    %   H = POLAR(...) returns a handle to the plotted object in H.
    %
    %   Example:
    %      t = 0 : .01 : 2 * pi;
    %      polar(t, sin(2 * t) .* cos(2 * t), '--r');
    %
    %   See also POLARPLOT, PLOT, LOGLOG, SEMILOGX, SEMILOGY.
    
    %   Copyright 1984-2015 MathWorks, Inc.
    
    % Parse possible Axes input
    
    global Rmax
    
    [cax, args, nargs] = axescheck(varargin{:});
    
    if nargs < 1
        error(message('MATLAB:narginchk:notEnoughInputs'));
    elseif nargs > 3
        error(message('MATLAB:narginchk:tooManyInputs'));
    end
    
    if isa(cax,'matlab.ui.control.UIAxes')
        error(message('MATLAB:ui:uiaxes:general'));
    end
                
    if nargs < 1 || nargs > 3
        error(message('MATLAB:polar:InvalidDataInputs'));
    elseif nargs == 2
        theta = args{1};
        rho = args{2};
        if ischar(rho)
            line_style = rho;
            rho = theta;
            [mr, nr] = size(rho);
            if mr == 1
                theta = 1 : nr;
            else
                th = (1 : mr)';
                theta = th(:, ones(1, nr));
            end
        else
            line_style = 'auto';
        end
    elseif nargs == 1
        theta = args{1};
        line_style = 'auto';
        rho = theta;
        [mr, nr] = size(rho);
        if mr == 1
            theta = 1 : nr;
        else
            th = (1 : mr)';
            theta = th(:, ones(1, nr));
        end
    else % nargs == 3
        [theta, rho, line_style] = deal(args{1 : 3});
    end
    if ischar(theta) || ischar(rho)
        error(message('MATLAB:polar:InvalidInputType'));
    end
    if ~isequal(size(theta), size(rho))
        error(message('MATLAB:polar:InvalidInputDimensions'));
    end
    try
        theta = full(double(theta));
        rho = full(double(rho));
    catch
        error(message('MATLAB:specgraph:private:specgraph:nonNumericInput'));
    end
    
    % get hold state
    cax = newplot(cax);
    
    next = lower(get(cax, 'NextPlot'));
    hold_state = ishold(cax);

    if isa(handle(cax),'matlab.graphics.axis.PolarAxes')
        error(message('MATLAB:polar:PolarAxes'));
    end
    
    % The grid color will be based on the axes background and grid color.
    axColor = cax.Color;
    if strcmp(axColor,'none')
        % If the axes is transparent, fall back to the parent container
        parent = cax.Parent;
        
        if isprop(parent,'BackgroundColor')
            % Panels and Tabs use BackgroundColor
            axColor = parent.BackgroundColor;
        else
            % Figures use Color
            axColor = parent.Color;
        end
        
        if strcmp(axColor,'none')
            % A figure/tab with Color none is black.
            axColor = [0 0 0];
        end
    end
    
    gridColor = cax.GridColor;
    gridAlpha = cax.GridAlpha;
    if strcmp(gridColor,'none')
        % Grid color is none, ignore transparency.
        tc = gridColor;
    else
        % Manually blend the color of the axes with the grid color to mimic
        % the effect of GridAlpha.
        tc = gridColor.*gridAlpha + axColor.*(1-gridAlpha);
    end
    ls = cax.GridLineStyle;
    
    % Hold on to current Text defaults, reset them to the
    % Axes' font attributes so tick marks use them.
    fAngle = get(cax, 'DefaultTextFontAngle');
    fName = get(cax, 'DefaultTextFontName');
    fSize = get(cax, 'DefaultTextFontSize');
    fWeight = get(cax, 'DefaultTextFontWeight');
    fUnits = get(cax, 'DefaultTextUnits');
    set(cax, ...
        'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
        'DefaultTextFontName', get(cax, 'FontName'), ...
        'DefaultTextFontSize', get(cax, 'FontSize'), ...
        'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
        'DefaultTextUnits', 'data');
    
    % only do grids if hold is off
    if ~hold_state
        
        % make a radial grid
        hold(cax, 'on');
        % ensure that Inf values don't enter into the limit calculation.
        arho = abs(rho(:));
        maxrho = max(arho(arho ~= Inf));
        hhh = line([-maxrho, -maxrho, maxrho, maxrho], [-maxrho, maxrho, maxrho, -maxrho], 'Parent', cax);
        set(cax, 'DataAspectRatio', [1, 1, 1], 'PlotBoxAspectRatioMode', 'auto');
        v = [get(cax, 'XLim') get(cax, 'YLim')];
        ticks = sum(get(cax, 'YTick') >= 0);
        delete(hhh);
        % check radial limits and ticks
        rmin = 0;
        rmax = Rmax;
        rticks = max(ticks - 1, 2);
        rticks = 2;
        if rticks > 5   % see if we can reduce the number
            if rem(rticks, 2) == 0
                rticks = rticks / 2;
            elseif rem(rticks, 3) == 0
                rticks = rticks / 3;
            end
        end
        
        % define a circle
        th = 0 : pi / 50 : 2 * pi;
        xunit = cos(th);
        yunit = sin(th);
        % now really force points on x/y axes to lie on them exactly
        inds = 1 : (length(th) - 1) / 4 : length(th);
        xunit(inds(2 : 2 : 4)) = zeros(2, 1);
        yunit(inds(1 : 2 : 5)) = zeros(3, 1);
        % plot background if necessary
        if ~ischar(get(cax, 'Color'))
            patch('XData', xunit * rmax, 'YData', yunit * rmax, ...
                'EdgeColor', tc, 'FaceColor', get(cax, 'Color'), ...
                'HandleVisibility', 'off', 'Parent', cax);
        end
        
        % draw radial circles
        c82 = cos(34 * pi / 180);
        s82 = sin(34 * pi / 180);
        rinc = (rmax - rmin) / rticks;
        for i = (rmin + rinc) : rinc : rmax
            hhh = line(xunit * i, yunit * i, 'LineStyle', ls, 'Color', tc, 'LineWidth', 1, ...
                'HandleVisibility', 'off', 'Parent', cax);
            text((i + rinc / 20) * c82, (i + rinc / 20) * s82, ...
                ['  ' num2str(round(i,2))], 'VerticalAlignment', 'bottom', ...
                'HandleVisibility', 'off', 'Parent', cax,'FontSize',24);
        end
        set(hhh, 'LineStyle', '-'); % Make outer circle solid
        
        % plot spokes
        th = (1 : 6) * 2 * pi / 12;
        cst = cos(th);
        snt = sin(th);
        cs = [-cst; cst];
        sn = [-snt; snt];
        line(rmax * cs, rmax * sn, 'LineStyle', ls, 'Color', tc, 'LineWidth', 1, ...
            'HandleVisibility', 'off', 'Parent', cax);
        
        % annotate spokes in degrees
        rt = 1.1 * rmax;
        for i = 1 : length(th)
%             text(rt * cst(i), rt * snt(i), int2str(i * 30),...
%                 'HorizontalAlignment', 'center', ...
%                 'HandleVisibility', 'off', 'Parent', cax);
            if i == length(th)
                loc = int2str(0);
            else
                loc = int2str(180 + i * 30);
            end
%             text(-rt * cst(i), -rt * snt(i), loc, 'HorizontalAlignment', 'center', ...
%                 'HandleVisibility', 'off', 'Parent', cax);
        end
        
        % set view to 2-D
        view(cax, 2);
        % set axis limits
        axis(cax, rmax * [-1, 1, -1.15, 1.15]);
    end
    
    % Reset defaults.
    set(cax, ...
        'DefaultTextFontAngle', fAngle , ...
        'DefaultTextFontName', fName , ...
        'DefaultTextFontSize', fSize, ...
        'DefaultTextFontWeight', fWeight, ...
        'DefaultTextUnits', fUnits );
    
    % transform data to Cartesian coordinates.
    xx = rho .* cos(theta);
    yy = rho .* sin(theta);
    
    % plot data on top of grid
    if strcmp(line_style, 'auto')
        q = plot(xx, yy, 'Parent', cax,'LineWidth',2);
    else
        q = plot(xx, yy, line_style, 'Parent', cax,'LineWidth',2);
    end
    
    if nargout == 1
        hpol = q;
    end
    
    if ~hold_state
        set(cax, 'DataAspectRatio', [1, 1, 1]), axis(cax, 'off');
        set(cax, 'NextPlot', next);
    end
    set(get(cax, 'XLabel'), 'Visible', 'on');
    set(get(cax, 'YLabel'), 'Visible', 'on');
    
    % Disable pan and zoom
    p = hggetbehavior(cax, 'Pan');
    p.Enable = false;
    z = hggetbehavior(cax, 'Zoom');
    z.Enable = false;
    
    if ~isempty(q) && ~isdeployed
        makemcode('RegisterHandle', cax, 'IgnoreHandle', q, 'FunctionName', 'polar');
    end
end