classdef ShowTimeseries < matlab.apps.AppBase
  
  % Properties
  properties (Access = public)
    
    UIFigure                    matlab.ui.Figure
    TabGroup                    matlab.ui.container.TabGroup
    OverviewTab                 matlab.ui.container.Tab
    ElementsDropDownLabel       matlab.ui.control.Label
    ElementsDropDown            matlab.ui.control.DropDown
    TimesyncinsLabel            matlab.ui.control.Label
    SpinnerTime                 matlab.ui.control.Spinner
    UIAxesOverview
    UIAxesResidual
    UISpinnerLabels
    UISpinners
    
  end % properties (Access = public)
  
  properties (Access = private)
    
    eecmodel
    
    tU_s
    U_V
    tI_s
    I_A
    
    I_idx;
    U_idx;
    
    timeSync_s = 0
    Ifit_A
    
  end % properties (Access = private)
  
  % Methods
  methods (Access = private, Hidden = true)
    
    function plotTimeseries(app)
      
      % Calculate dc response
      response = dc(app.eecmodel,app.tU_s(app.U_idx),app.Ifit_A);
      
      % Timeseries plot
      cla(app.UIAxesOverview);
      co = app.UIAxesOverview.ColorOrder;
      
      % Measurement points
      plot(app.UIAxesOverview,app.tU_s(app.U_idx),app.U_V(app.U_idx),'.');
      
      % EEC model result
      plot(app.UIAxesOverview,response.t_s(:),response.Umodel_V(:),'-k');
      
      % EEC model element voltages
      app.UIAxesOverview.ColorOrder = co;
      app.UIAxesOverview.ColorOrderIndex = 1;
      ocvIdx = contains(lower(response.Ulabels),'ocv');%strcmpi(response.Ulabels, 'ocv');
      if any(ocvIdx)
        hline = plot(app.UIAxesOverview,response.t_s(:),...
          [response.Uelements_V(:,ocvIdx) response.Uelements_V(:,~ocvIdx)+response.Uelements_V(:,ocvIdx)],'-');
      else
        hline = plot(app.UIAxesOverview,response.t_s(:),...
        response.Uelements_V(:,:),'-');
      end
      
      % Highlight selected element
      hline(strcmpi(app.ElementsDropDown.Value, app.ElementsDropDown.Items)).LineWidth = 2;
      
      % Format timeseries plot
      xtickformat(app.UIAxesOverview, 'mm:ss.SSS');
      legend(app.UIAxesOverview,['Reference','EEC model',response.Ulabels],...
        'Orientation','horizontal','Location','best','NumColumns',5);
      
      % Residual plot
      cla(app.UIAxesResidual);
      
      % Measurement - (EEC model voltage)
      plot(app.UIAxesResidual,response.t_s(:),app.U_V(app.U_idx)-response.Umodel_V(:),'.-');
      xtickformat(app.UIAxesResidual, 'mm:ss.SSS');
      
    end % function plotTimeseries(app)
    
  end % methods (Access = private, Hidden = true)
  
  methods (Access = private, Hidden = false)
    
    % Create UIFigure and components
    function createComponents(app)
      
      % Create UIFigure and hide until all components are created
      app.UIFigure = uifigure('Visible', 'off');
      app.UIFigure.Position = [100 100 640 480];
      app.UIFigure.Name = 'Time domain';
      
      % Create TabGroup
      app.TabGroup = uitabgroup(app.UIFigure);
      app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);
      app.TabGroup.Position = [1 1 640 480];
      
      % Create OverviewTab
      app.OverviewTab = uitab(app.TabGroup);
      app.OverviewTab.Title = 'Overview';
      
      % Create UIAxesOverview
      app.UIAxesOverview = axes(app.OverviewTab);
      title(app.UIAxesOverview, '')
      xlabel(app.UIAxesOverview, '')
      ylabel(app.UIAxesOverview, 'voltage in V')
      app.UIAxesOverview.Box = 'on';
      app.UIAxesOverview.NextPlot = 'add';
      app.UIAxesOverview.XGrid = 'on';
      app.UIAxesOverview.YGrid = 'on';
      app.UIAxesOverview.Units = 'points';
      % app.UIAxesOverview.Position = [1 171 534 284];
      % app.UIAxesOverview.ActivePositionProperty = 'outerposition';
      app.UIAxesOverview.OuterPosition = [1 171 534 284];
      
      % Create UIAxesResidual
      app.UIAxesResidual = axes(app.OverviewTab);
      title(app.UIAxesResidual, '')
      xlabel(app.UIAxesResidual, 'duration')
      ylabel(app.UIAxesResidual, 'residual voltage in V')
      app.UIAxesResidual.Box = 'on';
      app.UIAxesResidual.NextPlot = 'add';
      app.UIAxesResidual.XGrid = 'on';
      app.UIAxesResidual.YGrid = 'on';
      app.UIAxesResidual.Units = 'points';
      % app.UIAxesResidual.Position = [1 1 534 171];
      % app.UIAxesResidual.ActivePositionProperty = 'outerposition';
      app.UIAxesResidual.OuterPosition = [1 1 534 171];
      
      % Sync x axes of both graphs
      linkaxes([app.UIAxesOverview app.UIAxesResidual],'x');
      
      % Create TimesyncinsLabel
      app.TimesyncinsLabel = uilabel(app.OverviewTab);
      app.TimesyncinsLabel.Position = [534 433 86 22];
      app.TimesyncinsLabel.Text = 'Time sync in s';
      
      % Create SpinnerTime
      app.SpinnerTime = uispinner(app.OverviewTab);
      app.SpinnerTime.Step = 1e-6;
      app.SpinnerTime.Value = app.timeSync_s;
      app.SpinnerTime.ValueChangedFcn = createCallbackFcn(app, @SyncTimeValueChanged, true);
      app.SpinnerTime.Position = [534 412 100 22];
            
      % Create ElementsDropDownLabel
      app.ElementsDropDownLabel = uilabel(app.OverviewTab);
      app.ElementsDropDownLabel.Position = [534 354 91 22];
      app.ElementsDropDownLabel.Text = 'Model elements';
      
      % Create ElementsDropDown
      app.ElementsDropDown = uidropdown(app.OverviewTab);
      app.ElementsDropDown.ValueChangedFcn = createCallbackFcn(app, @ElementsDropDownValueChanged, true);
      app.ElementsDropDown.Position = [534 333 100 22];
      app.ElementsDropDown.Items = getElementLabels(app.eecmodel);
      
      % Create SpinnerLabel and Spinners
      n = max(cellfun(@length, {app.eecmodel.Elements.Description}));
      for i = 1:n
        % Create SpinnerLabel
        label = uilabel(app.OverviewTab);
        label.Visible = 'off';
        label.Position = [534 354-i*42 100 22];
        label.Text = '';
        app.UISpinnerLabels{i} = label;
        
        % Create Spinner
        spinner = uispinner(app.OverviewTab);
        spinner.Visible = 'off';
        spinner.Limits = [-Inf Inf];
        spinner.ValueChangedFcn = createCallbackFcn(app, @ElementValueChanged, true);
        spinner.Position = [534 333-i*42 100 22];
        app.UISpinners{i} = spinner;
      end
      
      % Show the figure after all components are created
      app.UIFigure.Visible = 'on';
      
    end % function createComponents(app)
    
    % Code that executes after component creation
    function startupFcn(app, varargin)
            
      % Run DropDown value change event
      ElementsDropDownValueChanged(app);
      
    end % function startupFcn(app, varargin)
    
    % Value changed function: SpinnerTime
    function SyncTimeValueChanged(app, ~)
      
      % Set time sync value
      app.timeSync_s = app.SpinnerTime.Value;
      
      % Trim time vector to overlapping area only
%       app.I_idx = (app.tI_s(:)+seconds(app.timeSync_s) >= max(app.tI_s(1)+seconds(app.timeSync_s),app.tU_s(1)) & app.tI_s(:)+seconds(app.timeSync_s) <= min(app.tI_s(end)+seconds(app.timeSync_s),app.tU_s(end)));
      app.U_idx = find(app.tU_s(:) >= max(app.tI_s(1)+seconds(app.timeSync_s),app.tU_s(1)) & app.tU_s(:) <= min(app.tI_s(end)+seconds(app.timeSync_s),app.tU_s(end)));
      
      % Interp1 to match new time vector
      app.Ifit_A = interp1(app.tI_s+seconds(app.timeSync_s),app.I_A,app.tU_s(app.U_idx));
%       if length(app.U_idx) ~= length(app.Ifit_A)
%         app.U_idx = app.U_idx();
%       end
      
      % Update plot
      plotTimeseries(app);
      
    end % function SyncTimeValueChanged(app, event)
    
    % Value changed function: ElementsDropDown
    function ElementsDropDownValueChanged(app, ~)
      
      % get selected element
      order = strcmpi(app.ElementsDropDown.Value, app.ElementsDropDown.Items);
      element = app.eecmodel.Elements(order);
      
      for i = 1:length(element.Description)
        % adjust label and spinner value
        app.UISpinnerLabels{i}.Text = element.Description{i};
        app.UISpinners{i}.Tooltip = element.Description{i};
        app.UISpinners{i}.Value = element.value(i);
        % nice stepsize depending on value
        [y,e,~] = engunits(abs(element.value(i)));
        app.UISpinners{i}.Step = 1e-1/e*10^(max(0,floor(log(y)/log(10))));
        
        % show label and spinner
        app.UISpinnerLabels{i}.Visible = 'on';
        app.UISpinners{i}.Visible = 'on';
      end
      % hide all other labels and spinners
      for j = i+1:length(app.UISpinnerLabels)
        app.UISpinnerLabels{j}.Visible = 'off';
        app.UISpinners{j}.Visible = 'off';
      end
      
      % refresh plot
      plotTimeseries(app);
      
    end % function ElementsDropDownValueChanged(app, event)
    
    % Callback function
    function ElementValueChanged(app, ~)
      
      % get selected element
      order = strcmpi(app.ElementsDropDown.Value, app.ElementsDropDown.Items);
      element = app.eecmodel.Elements(order);
      
      for i = 1:length(element.value)
        element.value(i) = app.UISpinners{i}.Value;
      end
      
      % set new element values
      setElementValues(app.eecmodel,element.value,'Element',find(order))
      
      % refresh plot
      plotTimeseries(app);
      
    end % function ElementValueChanged(app, event)
    
    % Callback function
    function SaveModelButtonPushed(app, ~)
      
      % export adjusted eecmodel values to "base" workspace
      title = 'Save adjusted eecmodel';
      answer = [];
      while isempty(answer)
        answer = inputdlg('Enter the variable name',title,[1 50],{'eecmodel'});
        if isempty(answer)
          return;
        end
        % check for valid variable name
        if isvarname(answer{1})
          % save eecmodel to base workspace
          assignin("base",answer{1},app.eecmodel);
        else
          warning('Invalid filename %s.', answer{1});
          answer = [];
        end
      end
      
    end % function SaveModelButtonPushed(app, event)
    
    % Selection change function: TabGroup
    function TabGroupSelectionChanged(app, ~)
      
      selectedTab = app.TabGroup.SelectedTab;
      
      switch lower(selectedTab.Title)
        case {'overview'}
          % refresh plot
          plotTimeseries(app);
          
      end
      
    end % function TabGroupSelectionChanged(app, event)
    
  end % methods (Access = private, Hidden = false)
  
  % App creation and deletion
  methods (Access = public, Hidden = false)
    
    % Construct app
    function app = ShowTimeseries(eecmodel,tU_s,U_V,tI_s,I_A,varargin)
      
      % Parse input parameters
      p = inputParser;
      p.addRequired("app",@(x)validateattributes(x,"ShowTimeseries",{"scalar"}));
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      p.addRequired("tU_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("U_V",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(tU_s)}));
      p.addRequired("tI_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("I_A",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(tI_s)}));
      p.addParameter("TimeSync",0,@(x)validateattributes(x,["numeric" "duration"],{"scalar"}));
      p.parse(app,eecmodel,tU_s,U_V,tI_s,I_A,varargin{:});
      
      % Copy eecmodel to properties
      app.eecmodel = p.Results.eecmodel;
      
      % Copy data vectors to properties
      app.U_V = p.Results.U_V;
      app.I_A = p.Results.I_A;
      
      % Convert double to seconds
      if ~isduration(p.Results.tU_s)
        app.tU_s = seconds(p.Results.tU_s);
      else
        app.tU_s = p.Results.tU_s;
      end
      if ~isduration(p.Results.tI_s)
        app.tI_s = seconds(p.Results.tI_s);
      else
        app.tI_s = p.Results.tI_s;
      end
      
      % Default time sync
      if isduration(p.Results.TimeSync)
        app.timeSync_s = seconds(p.Results.TimeSync);
      else
        app.timeSync_s = p.Results.TimeSync;
      end
      
      % Interpolate current waveform initially
%       app.I_idx = (app.tI_s(:)+seconds(app.timeSync_s) >= max(app.tI_s(1)+seconds(app.timeSync_s),app.tU_s(1)) & app.tI_s(:)+seconds(app.timeSync_s) <= min(app.tI_s(end)+seconds(app.timeSync_s),app.tU_s(end)));
      app.U_idx = (app.tU_s(:) >= max(app.tI_s(1)+seconds(app.timeSync_s),app.tU_s(1)) & app.tU_s(:) <= min(app.tI_s(end)+seconds(app.timeSync_s),app.tU_s(end)));
      
      app.Ifit_A = interp1(app.tI_s+seconds(app.timeSync_s),app.I_A,app.tU_s(app.U_idx));
      
      % Create UIFigure and components
      createComponents(app)
      
      % Register the app with App Designer
      registerApp(app, app.UIFigure)
      
      % Execute the startup function
      runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
      
      if nargout == 0
        clear app
      end
      
    end % Construct app
    
    % Destructor app
    function delete(app)
      
      % Delete UIFigure when app is deleted
      delete(app.UIFigure)
      
    end % Destructor app
    
  end % methods (Access = public, Hidden = false)
  
end