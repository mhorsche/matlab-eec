classdef eec < handle & matlab.mixin.CustomDisplay & matlab.mixin.SetGet
  % EQUIVALENT-ELECTRIC-CIRCUIT (EEC) model class
  %
  %   eecmodel = eec(ModelString);
  %   returns EEC model object
  %
  %   eecmodel = eec(ModelString,'NAME','VALUE')
  %   returns EEC model object with additional name-value arguments. For
  %   instance, you can specify InitialValues for default values of all EEC
  %   model elements.
  %
  %   Input Arguments:
  %   ----------------
  %
  %   ModelString      - List of EEC elements, can be separated (, + - _)
  %                      string or cell array containing one or multiple
  %                      elements of 'OCV' 'R' 'C' 'L' 'RC' or 'RL'.
  %
  %   Name-Value arguments:
  %   ---------------------
  %
  %   InitialValues    - Initial values for each EEC element. 'OCV', 'R', 
  %                      'C' and 'L' consist of one value (voltage, resistance, 
  %                      capacity or inductivity) while 'RC' and 'RL' have 
  %                      three (resistance, time constant, initial voltage/charge)
  %                      values.
  %   SpicePath        - Path to LTspice XVII executable to run spice
  %                      simulation in batch mode.
  %
  %   Methods:
  %   --------
  %
  %   dc               - Run DC analyse and return time domain response
  %   dcfit            - Optimize EEC elements with given time domain
  %                      signals (t_s, U_V, I_A)
  %   ac               - Run AC analyse and return frequency domain response
  %   acfit            - Optimize EEC elements with given frequency domain
  %                      signals (f_Hz, ZRe_Ohm, ZIm_Ohm)
  %   spice            - Run DC response in spice batch mode to veriify
  %                      analytic results
  %
  %   timeseries       - Open frontend application to adjust a predefined
  %                      model manually in the time domain to a given
  %                      timeseries.
  %   nyquist          - Open frontend application to adjust a predefined
  %                      model manually in the frequency domain to a given
  %                      nyquist plot.
  %
  %   Examples:
  %   ---------
  %
  %     Run without additional arguments opens input dialog to adjust initial values
  %   eecmodel = eec("R+RC");
  %
  %     Predefined initial values without input prompt
  %   eecmodel = eec("OCV,R,L,RC,RC,RC",...
  %     ...%             OCV     R       L        RC1             RC2             RC3
  %     'InitialValues',[3.6500  0.0300  500e-9   5.0e-3 1.0e-3 0 10e-3 100e-3 0  20e-3 10 0 ]);
  
  properties(Constant, Access = private, Hidden = false)
    
    % Valid EEC elements
    ValidElements = {'OCV','R','C','L','RC','RL'};

  end % properties(Constant, Access = private, Hidden = false)
  
  properties(Access = private, Hidden = false)
    
    % EEC model string
    ModelString;
    
    % Model got fitted
    IsFitted;
    
    % Processing time of last query
    QueryDuration;
    
    % EEC spice simulation results
    SpiceSimulation;
    
  end % properties(Access = private, Hidden = false)
  
  properties(Access = public, Hidden = false)
    
    % EEC elements
    Elements = struct('Order',{},'Type',{},'Description',{},'Equation',{},'TransferFunction',{},'TimeResponse',{},...
      'value',{},'initial',{},'min',{},'max',{});
    
    % LTspice XVII path
    SpicePath = "";
    % SpicePath = "/Applications/LTspice.app/Contents/MacOS/LTspice";
    % SpicePath = "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe";
    
  end % properties(Access = private, Hidden = true)
  
  % Static helper function: 
  % FILTEROUTLIER, PLOTRESPONSE, RLE
  methods(Static, Access = public, Hidden = false)
    
    function yy = filterOutlier(y,varargin)
      %FILTEROUTLIER
      % Filter outlier using multiple of standard deviation as bondaries.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   y        - Values (waveform) which should be filtered
      %
      %   Optional arguments:
      %   ---------------------
      %
      %   alpha    - Multiplier for standard deviation to calculate lower
      %              and upper boundaries. Alpha defaults to 0.5
      %
      %   Example:
      %   --------
      %
      %   yy = eec.filterOutlier(y,0.7);
      
      p = inputParser;
      p.KeepUnmatched = true;
      
      p.addRequired("y",@(x)validateattributes(x,["numeric"],{"nonempty"}));
      p.addOptional("alpha",0.5,@(x)validateattributes(x,["numeric"],{"nonempty"}));
      
      p.parse(y,varargin{:});
      
      % Get input values
      y = p.Results.y;
      yy = y;
      
      % Calculate lower and upper bondaries
      dY = p.Results.alpha * std(diff(y));
      lb = (smooth(y, 3) - dY);
      ub = (smooth(y, 3) + dY);
      
      % yy(yy > ub || yy < lb) = nan;
      tail = 0;
      for i = 2:length(yy)
        if (tail == 0 && (yy(i) > ub(i) || yy(i) < lb(i)))
          tail = i-1;
        end
        if (tail ~= 0 && (yy(i) < ub(i) && yy(i) > lb(i)))
          for j = tail+1:i-1
            yy(j) = y(tail) + (y(i) - y(tail)) / (i-tail) * (j-tail);
          end
          tail = 0;
        end
      end
        
    end % function yy = filterOutlier(y, alpha)
    
    function h = plotResponse(response,varargin)
      %PLOTRESPONSE
      % Simple overview figure to show any (DC or AC) response data.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   response     - Response structure from DC or AC analysis.
      %
      %   Optional arguments:
      %   ---------------------
      %
      %   handle       - Figure handle which should be used, otherwise new
      %                  figure will be created.
      %
      %   Example:
      %   --------
      %
      %   res = eecmodel.dc(time,current);
      %   hFig = eec.plotResponse(res);
      
      p = inputParser;
      p.KeepUnmatched = true;
      
      p.addRequired("response",@(x)validateattributes(x,["struct"],{"scalar"}));
      p.addParameter("handle",false,@(x)validateattributes(x,["logical" "handle"],{"scalar"}));
      p.parse(response,varargin{:});
      
      % Show DC response
      if all(isfield(response,["t_s" "I_A" "Umodel_V"]))
        
        % new figure or use given handle
        if ~ishandle(p.Results.handle)
          h = figure('Name', 'EEC - DC Response');
        else
          h = figure(p.Results.handle);
        end
        
        % new axes
        ax = axes();
        hold on; grid on; box on;
        
        % show current
        yyaxis left;
        ax.YAxis(1).Color = [0 0 0];
        plot(response.t_s,response.I_A,'-k');
        ylabel('current in A')
        
        % show voltage response
        yyaxis right;
        ax.YAxis(2).Color = [0 0 0];
        plot(response.t_s,response.Umodel_V, 'o-k');
        xtickformat('mm:ss.SSS');
        ylabel('voltage in V')
        
        % show voltage for each element
        if all(isfield(response,["Uelements_V" "Ulabels"]))
          ax.ColorOrder = lines(size(response.Uelements_V,2));
          ocvIdx = contains(lower(response.Ulabels),'ocv');
          if any(ocvIdx)
            % spice?
            if isfield(response,'tspice_s')
              plot(response.tspice_s(:),[response.Uelements_V(:,ocvIdx) response.Uelements_V(:,~ocvIdx)+response.Uelements_V(:,ocvIdx)]);
            else
              plot(response.t_s(:),[response.Uelements_V(:,ocvIdx) response.Uelements_V(:,~ocvIdx)+response.Uelements_V(:,ocvIdx)]);
            end
          else
            % spice?
            if isfield(response,'tspice_s')
              plot(response.tspice_s(:),response.Uelements_V(:,:));
            else
              plot(response.t_s(:),response.Uelements_V(:,:));
            end
          end
          legend(['Current' 'Voltage response' response.Ulabels],...
            'Orientation','horizontal','Location','northoutside');
        else
          legend("Current","Voltage response",...
            'Orientation','horizontal','Location','northoutside');
        end
          
      end
      
      % show AC response
      if all(isfield(response,["f_Hz" "ZRe_Ohm" "ZIm_Ohm"]))
        % new figure or use given handle
        if ~ishandle(p.Results.handle)
          h = figure('Name', 'EEC - AC Response');
        else
          h = figure(p.Results.handle);
        end
        
        % new axes
        ax = axes();
        hold on; grid on; box on;
        
        % show Nyquist
        p1 = plot(response.ZRe_Ohm,response.ZIm_Ohm,'o-k');
        xlabel('Re\{Z\} in Ohm')
        ylabel('Im\{Z\} in Ohm')
        ax.YDir = 'reverse';
        ax.DataAspectRatioMode = 'manual';
        ax.DataAspectRatio = [1 1 1];
        
        % adjust DataTip and show frequency for each point
        p1.DataTipTemplate.DataTipRows(1).Label = 'Re\{Z\}';
        p1.DataTipTemplate.DataTipRows(2).Label = 'Im\{Z\}';
        p1.DataTipTemplate.DataTipRows(3) = dataTipTextRow('f in Hz',response.f_Hz) ;
      
        % show Nyquist for each element
        if all(isfield(response,["Zelements_Ohm" "Zlabels"]))
          ax.ColorOrder = lines(size(response.Zelements_Ohm,2));
          
          % Do not show OCV elements, no influence in frequency domain
          ocvIdx = ~contains(lower(response.Zlabels),'ocv');
          pAll = plot(response.Zelements_Ohm(:,ocvIdx),'.-');
          
          for p1 = pAll'
            % adjust DataTip and show frequency for each point
            p1.DataTipTemplate.DataTipRows(1).Label = 'Re\{Z\}';
            p1.DataTipTemplate.DataTipRows(2).Label = 'Im\{Z\}';
            p1.DataTipTemplate.DataTipRows(3) = dataTipTextRow('f in Hz',response.f_Hz) ;
          end
          
          legend(['AC response' response.Zlabels(ocvIdx)],...
            'Orientation','horizontal','Location','northoutside');
        else
          legend("AC response",...
            'Orientation','horizontal','Location','northoutside');
        end
        
      end
      
    end % function h = plotResponse(response,varargin)
    
    function data = rle(val)
      %RLE
      % Run-length encoding (RLE) to compress a vector by splitting into
      % two vectors, one containing the values and the second the run
      % length. Additionally a third vector is returned, containing the 
      % index. This Index can also be retrieved by the cumulated sum of the 
      % length.
      % 
      % 
      % Based on code by Stefan Eireiner (2020). RLE de/encoding 
      % (https://www.mathworks.com/matlabcentral/fileexchange/4955-rle-de-encoding)
      % MATLAB Central File Exchange. Retrieved November 14, 2020.
      % 
      %   Input Arguments:
      %   ----------------
      %
      %   val      - Numeric vector which should be compressed (encoded)
      %              or cell array which should be decompressed (decoded)
      %
      %   Example:
      %   --------
      %
      %   res = eecmodel.rle([1 1 1 1 0 0 1 1 2 2 3 3])
      
      if iscell(val)
        
        % decoding
        i = cumsum([ 1 val{2} ]);
        j = zeros(1,i(end)-1);
        j(i(1:end-1)) = 1;
        data = val{1}(cumsum(j));
        
      else
        
        % encoding
        val = val(:)';
        i = [ find(val(1:end-1) ~= val(2:end)) length(val) ];
        data{3} = [1 i(1:end-1)+1];
        data{2} = diff([ 0 i ]);
        data{1} = val(i);
        
      end
      
    end % function data = rle(x)
    
  end % methods(Static, Access = public, Hidden = false)
  
  % Static private helper functions:
  % EXTRACTEXCEPTIONMESSAGE, SHOWEXCEPTIONASWARNING, SIZEDISP
  methods(Static, Access = private, Hidden = true)
    
    function str = extractExceptionMessage(e)
      %EXTRACTEXCEPTIONMESSAGE
      % Return exception message of given exception structure.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   exception     - Any exception structure.
      
      if strcmpi(e.identifier,'MATLAB:Java:GenericException')
        exceptionObj = e.ExceptionObject;
        str = char(exceptionObj.getMessage);
      else
        str = e.message;
      end
      
    end % function errmessage = extractExceptionMessage(e)
    
    function showExceptionAsWarning(e)
      %SHOWEXCEPTIONASWARNING
      % Print given exception as warning with backtrace off
      %
      %   Input Arguments:
      %   ----------------
      %
      %   exception     - Any exception structure.
      
      warning off backtrace
      warning(eec.extractExceptionMessage(e))
      warning on backtrace
      
    end % function showExceptionAsWarning(e)
    
    function str = sizeDisp(value)
      %SIZEDISP
      % Return size (eg. 1x1) as nice formatted string.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   value     - Variable which size should be displayed.
      
      str = strjoin(arrayfun(@num2str,size(value),'UniformOutput',false),'x');
      
    end % function str = sizeDisp(value)
    
  end % methods(Static, Access = private, Hidden = true)
  
  % Static private analytic functions: 
  % RCNUMERICAL, RLNUMERICAL
  methods(Static, Access = private, Hidden = false)
    
    function U = rcNumerical(values,t_s,I_A)
      %RCNUMERICAL
      % Helper function to calculate voltage response of RC-elements for
      % given current waveform numerically:
      %
      %           rising: U = U_max * (1 - exp( -t/\tau )) = (R * I) * (1 - exp( -(t+dt)/RC ))
      % inflection point: dt = -log((1 - exp( -(dt_rising) / \tau ))) * \tau == -log((1 - exp( -(dt_rising) / RC ))) * RC
      %            pause: U = U_end * exp( -(t+dt)/\tau )) = U_end * exp( -(t+dt)/RC )
      %
      %   Input Arguments:
      %   ----------------
      %
      %   values    - Values of RC-element, must be a 3x1 vector with
      %               following parameters:
      %                 (1) Ohmic resistance
      %                 (2) Time constant
      %                 (3) Precharge
      %   t_s       - Time vector in Seconds
      %   I_A       - Belonging current vector in Ampere
      
      U = values(1) * I_A; % baselevel: R * I
      U(1) = values(3);    % starting value
      for j = 2:length(t_s)
        U(j) = U(j) - (U(j) - U(j-1)) * exp( -( (t_s(j)-t_s(j-1)) / values(2) ) );
      end
      
    end % RC calculation: U = rcNumerical(values,t_s,I_A)
     
    function U = rlNumerical(values,t_s,I_A)
      %RLNUMERICAL
      % Helper function to calculate voltage response of RL-elements for
      % given current waveform numerically:
      %
      %           rising: U = L * dI/dt * (1 - exp( -t/\tau )) = L * dI/dt * (1 - exp( -t/(L/R) ))
      % inflection point: dt = -log((1 - exp( -(dt_rising) / \tau ))) * \tau == -log((1 - exp( -(dt_rising) / (L/R) ))) * (L/R)
      %            pause: U = L * dI/dt * exp( -(t+dt)/\tau )) = L * dI/dt * exp( -(t+dt)/(L/R) )
      %
      %   Input Arguments:
      %   ----------------
      %
      %   values    - Values of RL-element, must be a 3x1 vector with
      %               following parameters:
      %                 (1) Inductivity
      %                 (2) Time constant
      %                 (3) Precharge
      %   t_s       - Time vector in Seconds
      %   I_A       - Belonging current vector in Ampere
      
      U = values(1) * [0; diff(I_A)]./[1; diff(t_s)]; % baselevel: L * dI/dt
      U(1) = values(3);                               % starting value
      for j = 2:length(t_s)
        U(j) = U(j) - (U(j) - U(j-1)) * exp( -( (t_s(j)-t_s(j-1)) / values(2) ) );
      end
      
    end % RL calculation: U = rlNumerical(values,t_s,I_A)
    
  end % methods(Static, Access = private, Hidden = false)
  
  % Private helper functions: 
  % INPUTPROMPTVALUES, FITTINGBOUNDARIES, SETELEMENTSVALUES, SORTELEMENTS, 
  % CALCACRESPONSE, CALCDCRESPONSE, RUNSPICENETLIST
  methods(Access = private, Hidden = true)
    
    function inputPromptValues(eecmodel,definput,type,description,showDialog)
      %INPUTPROMPTVALUES
      % Show input prompt or set given values directly as output values.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   definput      - Default input dialog values
      %   type          - Which values should be set: 'initial', 'value',
      %                   'min', 'max'
      %   description   - Input dialog description
      %   showDialog    - Show input dialog or skip and use default values
      %                   without user confirmation
      
      % show input dialog for all element values
      prompt = cellfun(@(x)sprintf('%s for %s',description,x),[eecmodel.Elements.Description],'UniformOutput',false);
      
      dims = [1 75];
      opts.Interpreter = 'tex';
      title = sprintf('%s values',description);
      answer = [];
      while isempty(answer)
        
        if showDialog
          % show input dialog to define lower bondary values
          answer = inputdlg(prompt,title,dims,definput,opts);
        else
          % use input parameters without dialog
          answer = definput;
        end
        
        % check for cancelation
        if isempty(answer)
          return;
        end
        
        try
          % try to parse answers
          tail = 1;
          for i = 1:length(eecmodel.Elements)
            idx = [eecmodel.Elements.Order] == i;
            head = tail +( length(eecmodel.Elements(idx).(type)) - 1);
            val = str2double(answer(tail:head));
            eecmodel.Elements(idx).(type) = val(:)';
            tail = head +1;
            % set NaN to initial value
            eecmodel.Elements(idx).(type)(isnan(eecmodel.Elements(idx).(type))) = eecmodel.Elements(idx).initial(isnan(eecmodel.Elements(idx).(type)));
          end
          
        catch e
          % otherwise show warning and reopen input dialog
          warning(eec.extractExceptionMessage(e))
          definput = answer;
          answer = [];
        end
        
      end
      
    end % function inputPromptValues(eecmodel,definput,type,description,skipDialog)
    
    function [algorithm, lb, ub] = fittingBoundaries(eecmodel,varargin)
      %FITTINGBOUNDARIES
      % Get algorithm and boundaries depending on relative/absolute 
      % boundaries given by input parameters. Use either relative or
      % absolute boundaries.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %
      %   Name-Value arguments:
      %   ---------------------
      %
      %   Algorithm               - Which optimizer algorithm should be
      %                             used:
      %                             'trust-region-reflective' (default) or 'levenberg-marquardt'
      %   LowerBoundary           - Absolute values for lower boundary
      %   LowerBoundaryRelative   - Relative values for lower boundary, the
      %                             reference is Element.values
      %   UpperBoundary           - Absolute values for upper boundary
      %   UpperBoundaryRelative   - Relative values for upper boundary, the
      %                             reference is Element.values
      
      p = inputParser;
      p.KeepUnmatched = true;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      
      p.addParameter("Algorithm",'trust-region-reflective',@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      
      p.addParameter("LowerBoundary",[eecmodel.Elements.min],@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.min])}));
      p.addParameter("LowerBoundaryRelative",ones(size([eecmodel.Elements.min])),@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.min])}));
      p.addParameter("UpperBoundary",[eecmodel.Elements.max],@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.max])}));
      p.addParameter("UpperBoundaryRelative",ones(size([eecmodel.Elements.max])),@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.max])}));
      
      p.parse(eecmodel,varargin{:});
      
      % Absolute OR relative boundaries
      if all(~contains({'LowerBoundary', 'LowerBoundaryRelative'}, p.UsingDefaults)) || ...
          all(~contains({'UpperBoundary', 'UpperBoundaryRelative'}, p.UsingDefaults))
        error('eec:eec:InvalidBoundaries', ...
          'Can only handle absolute or relative boundaries.');
      end
      
      % check for selected algorithm {'levenberg-marquardt' 'trust-region-reflective'}
      switch lower(p.Results.Algorithm)
        case {'levenberg' 'levenberg-marquardt' 'marquardt'}
          algorithm = 'levenberg-marquardt';
          lb = [];
          ub = [];
          
        otherwise
          algorithm = 'trust-region-reflective';
          
          % handle lower boundary
          if ~any(contains(p.UsingDefaults, {'LowerBoundaryRelative'}))
            % use relative boundary with initial values as reference
            definput = arrayfun(@(x)sprintf('%g',x),[eecmodel.Elements.value] .* p.Results.LowerBoundaryRelative,'UniformOutput',false);
          else
            % use given (or default) lower boundary values
            definput = arrayfun(@(x)sprintf('%g',x),p.Results.LowerBoundary,'UniformOutput',false);
          end
          showDialog = sum(contains(p.UsingDefaults, {'LowerBoundary', 'LowerBoundaryRelative'})) == 2;
          inputPromptValues(eecmodel,definput,'min','Lower bondary',showDialog);
          
          % handle upper boundary
          if ~any(contains(p.UsingDefaults, {'UpperBoundaryRelative'}))
            % use relative boundary with initial values as reference
            definput = arrayfun(@(x)sprintf('%g',x),[eecmodel.Elements.value] .* p.Results.UpperBoundaryRelative,'UniformOutput',false);
          else
            % use given (or default) lower boundary values
            definput = arrayfun(@(x)sprintf('%g',x),p.Results.UpperBoundary,'UniformOutput',false);
          end          
          showDialog = sum(contains(p.UsingDefaults, {'UpperBoundary', 'UpperBoundaryRelative'})) == 2;
          inputPromptValues(eecmodel,definput,'max','Upper bondary',showDialog);
          
          % get boundaries (min/max values)
          lb = [eecmodel.Elements.min];
          ub = [eecmodel.Elements.max];
          
      end
      
    end % function [algorithm, lb, ub] = fittingBoundaries(eecmodel,varargin)
    
    function setElementsValues(eecmodel,values,type)
      %SETELEMENTSVALUES
      % Set given values to Elements.(type) values.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   values        - New element values
      %   type          - Which values should be set: 'initial', 'value',
      %                   'min', 'max'
      
      % always save as column vector
      values = values(:)';
      
      % assign values to each element
      tail = 1;
      for i = 1:length(eecmodel.Elements)
        idx = [eecmodel.Elements.Order] == i;
        head = tail +( length(eecmodel.Elements(idx).(type)) - 1);
        eecmodel.Elements(idx).(type) = values(tail:head);
        tail = head +1;
      end
      
    end % function setElementsValues(eecmodel,values,type)
    
    function sortElements(eecmodel)
      %SORTELEMENTS
      % Sort elements by its name and time constants, RC1's time constant
      % will be smaller than RC2's.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      
      % cache elements
      el = eecmodel.Elements;
      
      % sort 'RC' elements with its time constants (2nd value) in ascending order
      idx = find(strcmpi({el.Type}, 'RC'));
      val = [el(idx).value];
      [~, srt] = sort(val(2:3:end));
      
      % sort all fields except 'Order'
%       fnames = fieldnames(el);
      for i = 1:length(idx)
        eecmodel.Elements(idx(i)).value = el(idx(srt(i))).value;
%         for field = fnames(~strcmpi(fnames, 'Order'))'
%           eecmodel.Elements(idx(i)).(field{1}) = el(idx(srt(i))).(field{1});
%         end
      end
      
    end % function sortElements(eecmodel)
    
    function [Z_Ohm, Zelements_Ohm] = calcAcResponse(eecmodel,omega,values,varargin)
      %CALCACRESPONSE
      % Calcualte frequency domain response using each elements transfer 
      % function. Impedance results are splitted in real and complex part.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   omega         - Angular frequency vector in radians per second
      %                   (2*pi*f)
      %   values        - Element values
      
      % To use the trust-region-reflective algorithm, such as when you want to
      % include bounds, you must split the real and complex parts of the
      % coefficients into separate variables.
      % http://de.mathworks.com/help/optim/ug/fit-model-to-complex-data.html
      
      % Run transfer function and split results
      Zelements_Ohm = cellfun(@(tf)tf(values,omega),{eecmodel.Elements.TransferFunction},'UniformOutput',false);
      Zelements_Ohm = [Zelements_Ohm{:}];
      
      % Sum of all elements for real and imaginary part
      Z_Ohm = [sum(real(Zelements_Ohm),2) sum(imag(Zelements_Ohm),2)];
      
    end % function [Z_Ohm, Zelements_Ohm] = calcAcResponse(eecmodel,omega,values,varargin)
    
    function [U_V, Uelements_V] = calcDcResponse(eecmodel,t_s,I_A,values,varargin)
      %CALCDCRESPONSE
      % Calcualte time domain response using each elements transfer 
      % function.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   t_s           - Time vector in seconds
      %   I_A           - Belonging input current vector in Ampere
      %   values        - Element values

      % Run transfer function and split results
      Uelements_V = cellfun(@(tf)tf(values,t_s,I_A),{eecmodel.Elements.TimeResponse},'UniformOutput',false);
      Uelements_V = [Uelements_V{:}];
      
      % Filter outlier
      for i = 1:size(Uelements_V,2)
        Uelements_V(:,i) = eec.filterOutlier(Uelements_V(:,i),0.5);
      end
      
      % Sum of all over potentials of each single element
      U_V = sum(Uelements_V, 2, 'omitnan');
      
    end % function [U_V, Uelements_V] = calcDcResponse(eecmodel,t_s,I_A,values,varargin)
    
    function [U_V, sim] = runSpiceNetlist(eecmodel,t_s,I_A,varargin)
      %RUNSPICENETLIST
      % Generate spice netlist (eec.cir) and run Spice .tran command. The
      % results are read back into Matlab using LTspice2Matlab 
      % (https://github.com/PeterFeicht/ltspice2matlab) and parsed to
      % match response structure. The implementation is very basic for now.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   t_s           - Time vector in seconds
      %   I_A           - Belonging input current vector in Ampere
      %
      %   Optional arguments:
      %   ---------------------
      %
      %   values        - Element values
      % 
      %   Name-Value arguments:
      %   ---------------------
      %
      %   MaxStepsize           - Max stepsize parameter [dTmax] of .tran
      %                           command
      %   UseInitialConditions  - Set .uic flag of .tran command
      
      p = inputParser;
      p.KeepUnmatched = true;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      p.addRequired("t_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("I_A",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(t_s)}));
      p.addOptional("values",[],@(x)validateattributes(x,["numeric"],{"size" size([eecmodel.Elements.value])}));
      
      p.addParameter("MaxStepsize",[],@(x)validateattributes(x,["numeric"],{"scalar"}));
      p.addParameter("UseInitialConditions",false,@(x)validateattributes(x,["logical"],{"scalar"}));
      
      p.parse(eecmodel,t_s,I_A,varargin{:});
      
      % set new values for this iteration if necessary
      if ~isempty(p.Results.values)
        setElementsValues(eecmodel,p.Results.values,'value');
      end
      
      % generate new spice netlist *.cir, LTSpice XVII needs UTF16-LE encoding
      % all commands are documented here: http://ltwiki.org/
      netlistName = "eec.cir";
      fid = fopen(netlistName, "w+");
      fwrite(fid, ['* Simple EEC battery model: ' char(strjoin(eecmodel.ModelString, ' + ')) newline], "uint16");
      
      % current sorce
      % Syntax: Ixxx n+ n- <current> [AC=<amplitude>] [load]
      %         Ixxx n+ n- PULSE(Ioff Ion Tdelay Trise Tfall Ton Tperiod Ncycles)
      %         Ixxx n+ n- SINE(Ioffset Iamp Freq Td Theta Phi Ncycles)
      %         Ixxx n+ n- tbl=(<voltage, current>, <voltage, current>, ...)
      %         Ixxx n+ n- PWL(t1 i1 t2 i2 t3 i3...)
      fwrite(fid, sprintf('Ibat V001 GND PWL( %s)%s',sprintf("%e %e ",[t_s -I_A]'),newline), "uint16");
      
      % handle each element separately
      for i = 1:length(eecmodel.Elements)
        
        % get netname +/- and connect last element to ground
        
        if i < length(eecmodel.Elements)
          net = sprintf('V%03d V%03d',i,i+1);
        else
          net = sprintf('V%03d GND',i);
        end
        
        element = eecmodel.Elements([eecmodel.Elements.Order] == i);
        switch upper(element.Type)
          case {'OCV'}
            % Open circuit voltage
            if length(element.value) == 1
              % constant voltage
              % Syntax: Vxxx n+ n- <voltage> [AC=<amplitude>]
              %         + [Rser=<value>] [Cpar=<value>]
              fwrite(fid, sprintf('Vocv%d %s %g %s',i,net,element.value(1),newline), "uint16");
            else
              % arbitrary piece-wise linear voltage
              % Syntax: Vxxx n+ n- PWL(t1 v1 t2 v2 t3 v3...)
              %         For times before t1, the voltage is v1. For times between
              %         t1 and t2, the voltage varies linearly between v1 and v2.
              %         There can be any number of time, voltage points given. For
              %         times after the last time, the voltage is the last voltage.
              fwrite(fid, sprintf('Vocv%d %s PWL( %s)%s',i,net,sprintf("%e %e ",element.value),newline), "uint16");
            end
            
          case {'R'}
            % Resistor
            % Syntax: Rxxx n+ n- <value> [tc=tc1, tc2, ...]
            %         + [temp=<value>]
            fwrite(fid, sprintf('R%d %s %e %s',i,net,max(1e-12,element.value(1)),newline), "uint16");
            
          case {'C'}
            % Capacitor
            % Syntax: Cnnn n+ n- <capacitance> [ic=<value>]
            %         + [Rser=<value>] [Lser=<value>] [Rpar=<value>]
            %         + [Cpar=<value>] [m=<value>]
            %         + [RLshunt=<value>] [temp=<value>]
            fwrite(fid, sprintf('C%d %s %e Rser=0 %s',i,net,element.value(1),newline), "uint16");
            
          case {'L'}
            % Inductivity
            % Syntax: Lxxx n+ n- <inductance> [ic=<value>]
            %         + [Rser=<value>] [Rpar=<value>]
            %         + [Cpar=<value>] [m=<value>] [temp=<value>]
            fwrite(fid, sprintf('L%d %s %e Rser=0 %s',i,net,element.value(1),newline), "uint16");
            
          case {'RC'}
            % RC element: R = \tau / R
            fwrite(fid, sprintf('R%d %s %e %s',i,net,max(1e-12,element.value(1)),newline), "uint16");
            fwrite(fid, sprintf('C%d %s %e ic=%e %s',i,net,min(1e12,element.value(2)/element.value(1)),element.value(3),newline), "uint16");
            
          case {'RL'}
            % RL element: L = \tau * R
            fwrite(fid, sprintf('R%d %s %e %s',i,net,max(1e-12,element.value(1)),newline), "uint16");
            fwrite(fid, sprintf('L%d %s %e ic=%e %s',i,net,element.value(2)*element.value(1),element.value(3),newline), "uint16");
            
        end
        
      end
      
      % Perform a Nonlinear Transient Analysis
      % Syntax: .tran <Tstep> <Tstop> [Tstart [dTmax]] [modifiers]
      %   or    .tran <Tstop> [modifiers]
      uic = '';
      if p.Results.UseInitialConditions
        uic = 'uic ';
      end
      fwrite(fid, sprintf('.tran %e %s%s',max(t_s),uic,newline), "uint16");
      
      % .OPTIONS -- Set Simulator Options
      %     Keyword       | Data Type | Default Value | Description
      % ------------------+-----------+---------------+---------------------------------------------------------------------------------------------------------------
      %     baudrate      |    Num    |     (none)    | Used for eye diagrams. Tells the waveform viewer how to wrap the abscissa time to overlay the bit transitions.
      %     abstol        |    Num    |      1pA      | Absolute current error tolerance
      %     chgtol        |    Num    |      10fC     | Absolute charge tolerance
      %     cshunt        |    Num    |       0.      | Optional capacitance added from every node to ground
      %     cshuntintern  |    Num    |     cshunt    | Optional capacitance added from every device internal node to ground.
      %     defad         |    Num    |       0.      | Default MOS drain diffusion area
      %     defas         |    Num    |       0.      | Default MOS source diffusion area
      %     defl          |    Num    |     100um     | Default MOS channel length
      %     defw          |    Num    |     100um     | Default MOS channel width
      %     delay         |    Num    |       0.      | Used for eye diagrams. Shifts the bit transitions in the diagram.
      %     fastaccess    |   flag    |     false     | Convert to fastaccess file format at end of simulation.
      %     flagloads     |   flag    |     false     | Flags external current sources as loads.
      %     Gmin          |    Num    |     1e-12     | Conductance added to every PN junction to aid convergence.
      %     gminsteps     |    Num    |       25      | Set to zero to prevent gminstepping for the initial DC solution.
      %     gshunt        |    Num    |       0.      | Optional conductance added from every node to ground.
      %     itl1          |    Num    |      100      | DC iteration count limit.
      %     itl2          |    Num    |       50      | DC transfer curve iteration count limit.
      %     itl4          |    Num    |       10      | Transient analysis time point iteration count limit
      %     itl6          |    Num    |       25      | Set to zero to prevent source stepping for the initial DC solution.
      %     srcsteps      |    Num    |       25      | Alternative name for itl6.
      %     maxclocks     |    Num    |     Infin.    | maximum number of clock cycles to save
      %     maxstep       |    Num    |     Infin.    | Maximum step size for transient analysis
      %     meascplxfmt   |  string   |      bode     | Complex number format of .meas statement results. One of "polar", "cartesian", or "bode".
      %     measdgt       |    Num    |       6       | Number of significant figures used for .measure statement output.
      %     method        |  string   |      trap     | Numerical integration method, either trapezoidal or Gear
      %     minclocks     |    Num    |       10      | minimum number of clock cycles to save
      %     MinDeltaGmin  |    Num    |      1e-4     | Sets a limit for termination of adaptive gmin stepping.
      %     nomarch       |   flag    |     false     | Do not plot marching waveforms
      %     noopiter      |   flag    |     false     | Go directly to gmin stepping.
      %     numdgt        |    Num    |       6       | Historically "numdgt" was used to set the number of significant figures used for output data. In LTspice, if "numdgt" is set to be > 6, double precision is used for dependent variable data.
      %     pivrel        |    Num    |     1e-3      | Relative ratio between the largest column entry and an acceptable pivot value.
      %     pivtol        |    Num    |     1e-13     | Absolute minimum value for a matrix entry to be accepted as a pivot.
      %     reltol        |    Num    |      001      | Relative error tolerance.
      %     srcstepmethod |    Num    |       0       | Which source stepping algorithm to start with.
      %     sstol         |    Num    |      001      | Relative error for steady-state detection.
      %     startclocks   |    Num    |       5       | Number of clock cycles to wait before looking for steady state.
      %     temp          |    Num    |     27uC      | Default temperature for circuit element instances that don't specify temperature.
      %     tnom          |    Num    |     27uC      | Default temperature at which device parameters were measured for models that don't specify this temperature.
      %     topologycheck |    Num    |       1       | Set to zero to skip check for floating nodes, loops of voltage sources, and non-physical transformer winding topology
      %     trtol         |    Num    |      1.0      | Set the transient error tolerance. This parameter is an estimate of the factor by which the actual truncation error is overestimated.
      %     trytocompact  |    Num    |       1       | When non-zero, the simulator tries to condense LTRA transmission lines' history of input voltages and currents.
      %     vntol         |    Num    |      1uV      | Sets the absolute voltage error tolerance.
      %     plotreltol    |    Num    |     0025      | Sets the relative error tolerance for waveform compression.
      %     plotvntol     |    Num    |     10uV      | Sets the absolute voltage error tolerance for waveform compression.
      %     plotabstol    |    Num    |      1nA      | Sets the absolute current error tolerance for waveform compression.
      %     plotwinsize   |    Num    |      300      | Number of data points to compress in one window. Set to zero to disable compression.
      %     ptrantau      |    Num    |       1       | Characteristic source start-up time for a damped pseudo transient analysis to find the operating point. Set to zero to disable pseudo transient.
      %     ptranmax      |    Num    |       0       | If set non-zero, that time of the damped pseudo transient analysis is used as the operating point whether the circuit has settled or not.
      if ~isempty(p.Results.MaxStepsize)
        fwrite(fid, sprintf('.options maxstep=%e %s',p.Results.MaxStepsize,newline), "uint16");
      end
      
      % Annotate the Subcircuit Pin Names to the Port Currents
      fwrite(fid, ['.backanno' newline], "uint16");
      fwrite(fid, ['.end' newline], "uint16");
      
      % close and save file
      fclose(fid);
      
      % run spice simulation
      tSim = tic;
      res = system(sprintf('"%s" -b %s',eecmodel.SpicePath,netlistName));
      if res > 0
        eecmodel.SpiceSimulation = [];
         ME = MException('eec:runSpiceNetlist:SpiceSimulationError', ...
          'Spice simulation failed (error %d).',res);
        throw(ME);
      end
      sim = LTspice2Matlab(strrep(netlistName,'.cir','.raw'));
      sim.ComputationTime = toc(tSim);
      
      % extract voltage with equal time vector for lsqnonlin and other
      % fitting algorithms
      U_V = interp1(sim.time_vect, sim.variable_mat(strcmpi(sim.variable_name_list, 'V(v001)'),:),t_s,'spline');
      
      % lsqcurvefit does not accept NaN values so keep them zero
      U_V(isnan(U_V)) = 0;
    end % function [U_V, sim] = runSpiceNetlist(eecmodel,t_s,I_A,varargin)
    
  end % methods(Static, Access = public, Hidden = false)
  
  % Public helper functions: 
  % GETELEMENTLABELS, SETELEMENTVALUES, SYNC
  methods(Access = public, Hidden = false)
    
    function labels = getElementLabels(eecmodel)
      %GETELEMENTLABELS
      % Get element labels as cell array in the same order as the result
      % vectors.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      
      labels = {eecmodel.Elements.Type};
      for i = 1:length(labels)
        labels{i} = sprintf('%s%d',labels{i},sum(contains(labels(1:i), labels{i})));
      end
      
    end % function labels = getElementLabels(eecmodel)
    
    function setElementValues(eecmodel,values,varargin)
      %SETELEMENTVALUES
      % Helper function to simply set element values.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   values        - New element value(s)
      % 
      %   Name-Value arguments:
      %   ---------------------
      %
      %   type          - Which values should be set, defaults to 'value': 
      %                   'initial', 'value', 'min', 'max'
      %   Element       - Element name
      %
      %   Example:
      %   --------
      %
      %   setElementValues(eecmodel,3.532,"Element",'OCV',"Type",'value');
      %     or shorthand:
      %   setElementValues(eecmodel,'OCV',3.532)
      
      p = inputParser;
      p.KeepUnmatched = true;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      p.addRequired("values",@(x)validateattributes(x,["numeric" "string" "char"],{"nonempty"}));
      p.addOptional("value",[],@(x)validateattributes(x,["numeric"],{"nonempty"}));
      p.addParameter("Type","value",@(x)validateattributes(x,["string" "char"],{"scalartext" "nonempty"}));
      p.addParameter("Element",[],@(x)validateattributes(x,["numeric" "string" "char"],{"nonempty"}));
      
      p.parse(eecmodel,values,varargin{:});
      
      % check for short input version
      if ~isnumeric(p.Results.values)
        element = char(p.Results.values);
        values = p.Results.value;
      else
        element = p.Results.Element;
        values = p.Results.values;
      end
      
      % Type can be 'value', 'initial', 'min' or 'max' or a unique substr
      valid = {'value' 'initial' 'min' 'max'};
      if sum(contains(valid,p.Results.Type)) ~= 1
        ME = MException('eec:setElementValues:IncorrectDataValue', ...
          'Valid values for type are "%s".', strjoin(valid, '", "'));
        throw(ME);
      else
        type = lower(valid{contains(valid,p.Results.Type)});
      end
      
      % check for single element
      if ~isempty(element)
        % get single element
        if isnumeric(element)
          % number must be order value
          idx = find([eecmodel.Elements.Order] == element);
        else
          % element type (+ consecutive number), eg. 'R', 'RC2'
          tmp = regexp(char(element),'(?<type>[^\d]*)(?<nth>\d*)','names');
          if isempty(tmp.nth)
            tmp.nth = 1;
          else
            tmp.nth = str2double(tmp.nth);
          end
          idx = find(strcmpi({eecmodel.Elements.Type},tmp.type));
          if isempty(idx)
            ME = MException('eec:setElementValues:InvalidElement', ...
              'Invalid element name %s. Please use order number or type name (+ consecutive number), eg. "R", "RC2".',...
              element);
            throw(ME);
          end
          idx = idx(tmp.nth);
        end
        % check for valid index
        if isempty(idx) || idx <= 0 || idx > length(eecmodel.Elements)
          ME = MException('eec:setElementValues:InvalidElement', ...
            'Invalid index %d for given element %s. Please use order number or type name (+ consecutive number), eg. "R", "RC2".',...
            idx,element);
          throw(ME);
        end
        % check values length
        if length(values) ~= length(eecmodel.Elements(idx).(type))
          ME = MException('eec:setElementValues:IncorrectDataSize', ...
            'To set %s for %s%d value length [%d] must match with elements parameter size [%d].',...
            type,eecmodel.Elements(idx).type,eecmodel.Elements(idx).Order,length(values),length(eecmodel.Elements(idx).(type)));
          throw(ME);
        end
        % always save as row vector
        eecmodel.Elements(idx).(type) = values(:)';
      else
        % set all values at once if value size matches
        if length(values) ~= length([eecmodel.Elements.(type)])
          ME = MException('eec:setElementValues:IncorrectDataSize', ...
            'To set %s for all elements value length [%d] must match with elements parameter size [%d].',type,length(values),length([eecmodel.Elements.(type)]));
          throw(ME);
        end
        setElementsValues(eecmodel,values,type);
      end
      
    end % function setElementValues(eecmodel,values,varargin)
  
    function res = sync(eecmodel,tI_s,I_A,tU_s,U_V,varargin)
      %SYNC 
      % Synchronize submitted current/voltage waveform. Returns current/
      % voltage waveform with combined time vector. This function is very
      % basic and should be used carefully.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   tI_s          - Current time vector in seconds
      %   I_A           - Current vector (must match with time vector) or
      %                   constant current value in Ampere
      %   tU_s          - Voltage time vector in seconds
      %   U_V           - Voltage vector as observed output vector in Volt
      %
      %   Name-Value arguments:
      %   ---------------------
      %
      %   Iterations    - Number of iterations to minimize error. Default
      %                   value is 100
      %   doPlot        - Plot DC response. This can be either a logical
      %                   value or a figure handle which should be used
      %                   instead of a new figure.
      %
      %   Example:
      %   --------
      %
      %   res = eecmodel.sync(tI_s,I_A,tU_s,U_V)
      
      p = inputParser;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      
      p.addRequired("tI_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("I_A",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(tI_s)}));
      p.addRequired("tU_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("U_V",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(tU_s)}));

      p.addParameter("Iterations",100,@(x)validateattributes(x,["numeric"],{"nonempty" "nonzero" "nonnegative"}));
      p.addParameter("doPlot",false,@(x)validateattributes(x,["logical" "handle"],{"scalar"}));
      
      p.parse(eecmodel,tI_s,I_A,tU_s,U_V,varargin{:});
      
      U_V = p.Results.U_V(:);
      I_A = p.Results.I_A(:);
      
      % convert duration to seconds
      if isduration(p.Results.tI_s)
        tI_s = seconds(p.Results.tI_s(:));
      else
        tI_s = p.Results.tI_s(:);
      end
      if isduration(p.Results.tU_s)
        tU_s = seconds(p.Results.tU_s(:));
      else
        tU_s = p.Results.tU_s(:);
      end
      
%       idxU = (tU_s >= max(tI_s(1),tU_s(1)) & tU_s <= min(tI_s(end),tU_s(end)));
      idxI = (tI_s >= max(tI_s(1),tU_s(1)) & tI_s <= min(tI_s(end),tU_s(end)));
      t_s = tI_s(idxI);
      
      % synchronize current/voltage signals using model parameters
      Umodel_V = calcDcResponse(eecmodel,t_s,I_A(idxI),[eecmodel.Elements.value]);
      Umeas_V = interp1(tU_s,U_V,t_s,'spline');
      
      if any(strcmp({eecmodel.Elements.Type},'L'))
        
        % If model has inductivity, use diff(U) as threshold
        dUmeas_V = smooth(diff(Umeas_V),ceil(1e-6/min(diff(t_s))));
        dUmodel_V = smooth(diff(Umodel_V),ceil(1e-6/min(diff(t_s))));
        
        % Get trigger values where signal drops below dU/2
        RLE = eec.rle(dUmeas_V < min(dUmeas_V)/2);
        index = RLE{3}(RLE{1} == 1);
        [~, idx] = max(RLE{2}(RLE{1} == 1));
        Umeas_shift = index(idx);
        
        RLE = eec.rle(dUmodel_V < min(dUmodel_V)/2);
        index = RLE{3}(RLE{1} == 1);
        [~, idx] = max(RLE{2}(RLE{1} == 1));
        Umodel_shift = index(idx);
        
        tsync_s = t_s(Umeas_shift)-t_s(Umodel_shift);
        
        % Plot results?
        if p.Results.doPlot ~= false
          % new figure or use given handle
          if ~ishandle(p.Results.doPlot)
            h = figure('Name', 'EEC - Synchronisation');
          else
            h = figure(p.Results.doPlot);
          end
                    
          % new axes
          ax = axes();
          hold on; grid on; box on;
          
          % show voltage response
          yyaxis left;
          ax.YAxis(1).Color = [0 0 0];
          plot(t_s,Umeas_V,'-','Color',[0 0.4470 0.7410]);
          plot(t_s,Umodel_V-mean(Umodel_V-Umeas_V),'-','Color',[0.8500 0.3250 0.0980]);
          ylabel('voltage in V')
          
          plot([t_s(Umeas_shift) t_s(Umeas_shift)],ylim,'-','Color',[0 0.4470 0.7410]);
          plot([t_s(Umodel_shift) t_s(Umodel_shift)],ylim,'-','Color',[0.8500 0.3250 0.0980]);
          
          % show differential voltage response
          yyaxis right;
          ax.YAxis(2).Color = [0 0 0];
          plot(t_s(1:end-1),dUmeas_V,'--','Color',[0 0.4470 0.7410]);
          hold on;
          plot(t_s(1:end-1),dUmodel_V,'--','Color',[0.8500 0.3250 0.0980]);
          ylabel('differential voltage in V/s')
                    
        end
        
      else
        
        % Synchronize current/voltage by minimizing error between voltage
        % response and given voltage waveform
        Umodel_V = Umodel_V - mean(Umodel_V);
        Umeas_V = Umeas_V - mean(Umeas_V);
        
        % Get rought trigger values where signal drops below (Umin + Umax)/2
        th = (max(Umeas_V)+min(Umeas_V))/2;
        Umeas_shift = find(Umeas_V < th,1,'first');
        Umodel_shift = find(Umodel_V < th,1,'first');
        tsync_s = t_s(Umeas_shift)-t_s(Umodel_shift);
        
        % Shift voltage response and voltage waveform to minimize error
        dir = -1;
        step = min(diff(t_s));
        dUrms = NaN(p.Results.Iterations,1);
        for i = 1:p.Results.Iterations
          if tsync_s > 0
            dU = (Umeas_V(1:sum(t_s > t_s(1)+tsync_s)) - Umodel_V(t_s > t_s(1)+tsync_s));
          else
            dU = (Umeas_V(t_s > t_s(1)-tsync_s)) - Umodel_V(1:sum(t_s > t_s(1)-tsync_s));
          end
          dUrms(i) = rms(dU);

          % Residium getting higher
          if i > 1 && dUrms(i) > dUrms(i-1)
            break;
          end
          tsync_s = tsync_s+step*dir;
        end
        
      end % if any(strcmp({eecmodel.Elements.Type},'L'))
      
      res.tsync_s = seconds(tsync_s);
      res.t_s = seconds(t_s);
      res.U_V = interp1(tU_s,U_V,t_s+tsync_s);
      res.I_A = I_A(idxI);
      
    end % Sync current/voltage: function res = sync(eecmodel,tI_s,I_A,tU_s,U_V,varargin)
    
  end % methods(Access = public, Hidden = true)
  
  % Public EEC functions: 
  % DC, DCFIT, AC, ACFIT, SPICE
  methods(Access = public, Hidden = false)
    
    function res = dc(eecmodel,t_s,I_A,varargin)
      %DC analysis
      % Calculate time domain response for given time and current.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   t_s           - Time vector in seconds (or duration)
      %   I_A           - Current vector (must match with time vector) or
      %                   constant current value in Ampere
      %
      %   Name-Value arguments:
      %   ---------------------
      %
      %   doPlot        - Plot DC response. This can be either a logical
      %                   value or a figure handle which should be used
      %                   instead of a new figure.
      % 
      %   Example:
      %   --------
      %
      %     Constant current value:
      %   res = eecmodel.dc(linspace(0,10,101),-2,"doPlot",true);
      %     or current waveform:
      %   res = eecmodel.dc(linspace(0,10,101),sin(2*pi*0:100));
      
      p = inputParser;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      p.addRequired("t_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("I_A",@(x)validateattributes(x,["numeric"],{"nonempty"}));
      
      p.addParameter("doPlot",false,@(x)validateattributes(x,["logical" "handle"],{"scalar"}));
      
      p.parse(eecmodel,t_s,I_A,varargin{:});
      
      % convert duration to seconds
      if isduration(p.Results.t_s)
        t_s = seconds(p.Results.t_s(:));
      else
        t_s = p.Results.t_s(:);
      end
      
      % check for continuous or contant value
      if isscalar(p.Results.I_A)
        I_A = ones(size(t_s)) .* p.Results.I_A;
      else
        if ~isequal(size(p.Results.I_A(:)), size(t_s))
          error('eec:eec:IncorrectDataSize', ...
            'Invalid vector length, "t_s" [%s] and "I_A" [%s] must match in size.', ...
            eec.sizeDisp(t_s), eec.sizeDisp(p.Results.I_A(:)));
        else
          I_A = p.Results.I_A(:);
        end
      end
      idx = isnan(I_A);
      I_A(idx) = 0;
      
      % calculate DC response
      [Umodel_V, Uelements_V] = calcDcResponse(eecmodel,...
        t_s,I_A,[eecmodel.Elements.value],varargin{:});
      
      % copy duration, voltage and current to result
      res.t_s = seconds(t_s(:));
      res.I_A = I_A(:);
      res.Umodel_V = Umodel_V(:);
      res.Uelements_V = Uelements_V;
      res.Ulabels = getElementLabels(eecmodel);
      
      % Plot results?
      if p.Results.doPlot ~= false
        eecmodel.plotResponse(res,'handle',p.Results.doPlot);
      end
      
    end % DC response: res = dc(eecmodel,t_s,I_A,varargin)
    
    function res = dcfit(eecmodel,tI_s,I_A,tU_s,Umodel_V,varargin)
      %DCFIT 
      % Run lsqcurvefit for given time/current and minimize error to
      % submitted voltage waveform.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   tI_s          - Current time vector in seconds
      %   I_A           - Current vector (must match with time vector) or
      %                   constant current value in Ampere
      %   tU_s          - Voltage time vector in seconds
      %   U_V           - Voltage vector as observed output vector in Volt
      %
      %   Name-Value arguments:
      %   ---------------------
      %
      %   Method                  - Use Matlab or Spice to calculate
      %                             waveform for each iteration: 
      %                             'analytical' (default) or 'spice'
      %   Repetitions             - Number of repetitions of lsqcurvefit 
      %                             command. The optimizer iterations are
      %                             set separately. Default value is 4
      %
      %   Algorithm               - Which optimizer algorithm should be
      %                             used:
      %                             'trust-region-reflective' (default) or 'levenberg-marquardt'
      %   LowerBoundary           - Absolute values for lower boundary
      %   LowerBoundaryRelative   - Relative values for lower boundary, the
      %                             reference is Element.values
      %   UpperBoundary           - Absolute values for upper boundary
      %   UpperBoundaryRelative   - Relative values for upper boundary, the
      %                             reference is Element.values
      %
      %   SetInitials             - Use current values (starting point) and
      %                             overwrite initial values
      %   UseInitials             - Use initial values as starting point
      %   Display                 - Display options of lsqcurvefit command,
      %                             default is 'off', use 'iter' or 'final'
      %                             for additional information
      %   
      %   Example:
      %   --------
      %
      %   res = eecmodel.fit(linspace(0,10,101),[0 repmat(1,1,100)],...
      %      linspace(0,10,101),exp(-linspace(0,10,101)/1));
      %
      %   res = eecmodel.fit(tI_s,I_A,tU_s,U_V,'Algorithm','levenberg-marquardt')
      
      p = inputParser;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      p.addRequired("tI_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("I_A",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(tI_s)}));
      p.addRequired("tU_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("U_V",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(tU_s)}));
      
      p.addParameter("Method",'analytical',@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      p.addParameter("Repetitions",4,@(x)validateattributes(x,["numeric"],{"nonempty" "nonzero" "nonnegative"}));
      
%       p.addParameter("MaxStepsize",[],@(x)validateattributes(x,["numeric"],{"scalar"}));
%       p.addParameter("UseInitialConditions",false,@(x)validateattributes(x,["logical"],{"scalar"}));
      
      p.addParameter("Algorithm",'trust-region-reflective',@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      p.addParameter("LowerBoundary",[eecmodel.Elements.min],@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.min])}));
      p.addParameter("LowerBoundaryRelative",ones(size([eecmodel.Elements.min])),@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.min])}));
      p.addParameter("UpperBoundary",[eecmodel.Elements.max],@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.max])}));
      p.addParameter("UpperBoundaryRelative",ones(size([eecmodel.Elements.max])),@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.max])}));
      
      p.addParameter("SetInitials",false,@(x)validateattributes(x,["logical"],{"scalar"}));
      p.addParameter("UseInitials",false,@(x)validateattributes(x,["logical"],{"scalar"}));
      
      p.addParameter("Display",'off',@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      
      p.parse(eecmodel,tI_s,I_A,tU_s,Umodel_V,varargin{:});
      
      Umodel_V = p.Results.U_V(:);
      I_A = p.Results.I_A(:);
      if isscalar(I_A)
        I_A = ones(size(tI_s)) .* I_A;
      else
        if ~isequal(size(I_A), size(tI_s))
          error('eec:eec:IncorrectDataSize', ...
            'Invalid vector length, "tI_s" [%s] and "I_A" [%s] must match in size.', ...
            eec.sizeDisp(tI_s), eec.sizeDisp(I_A));
        end
      end
      I_A(isnan(I_A)) = 0;
      
      % Set initials to current values
      if logical(p.Results.SetInitials)
        setElementsValues(eecmodel,[eecmodel.Elements.value],'initial');
      end
      
      % Use initials for first iteration
      if logical(p.Results.UseInitials)
        setElementsValues(eecmodel,[eecmodel.Elements.initial],'value');
      end
      
      % Get algorithm and boundaries if needed
      [algorithm, lb, ub] = fittingBoundaries(eecmodel,varargin{:});
      
      % Set lsqcurvefit options
      options = optimoptions(@lsqcurvefit,...
        'Algorithm', algorithm,...
        'TolX', 1e-20, ...
        'TolFun', 1e-20,...
        'Display', p.Results.Display, ...
        'MaxIter', 100, ...
        'MaxFunEvals', 10e3, ...
        'FunValCheck', 'on',...
        'UseParallel', false);
      
      tic;
      for i = 1:p.Results.Repetitions
        % synchronize current/voltage signals
        syn = sync(eecmodel,tI_s,I_A,tU_s,Umodel_V);

        % convert duration to seconds
        if isduration(syn.t_s)
          syn.t_s = seconds(syn.t_s);
        end
                
        % Omit NaN values, otherwise lsqnonline will fail
        idx = ~isnan(syn.I_A) & ~isnan(syn.U_V);
        
        % check for selected method {'spice' 'analytical'}
        switch lower(p.Results.Method)
          case {'spice'}
            % run spice simulation for each iteration
            fun = @(values,t_s)runSpiceNetlist(eecmodel,t_s,-syn.I_A(idx),values); %,varargin{:}

          otherwise
            % run step response calculation for each iteration
            fun = @(values,t_s)calcDcResponse(eecmodel,t_s,syn.I_A(idx),values); %,varargin{:}

        end
      
        % lsqcurvefit: Fit impedance in the time domain
        [fittedParameters, ~, ~, exitflag, output] = lsqcurvefit(fun,...
          [eecmodel.Elements.value],syn.t_s(idx),syn.U_V(idx),...
          lb,ub,options);
        output.Exitflag = exitflag;
      end
      output.ComputationTime = toc;
        
      % set fitted values
      setElementsValues(eecmodel,fittedParameters,'value');
      
      % sort RC elements time constants in ascending order
      sortElements(eecmodel);
      
      % calculate DC response for fitted model
      [Umodel_V, Uelements_V] = fun([eecmodel.Elements.value],syn.t_s(idx));
      
      % copy duration, voltage and current to result
      res.tsync_s = syn.tsync_s;
      res.t_s = seconds(syn.t_s(:));
      res.I_A = syn.I_A(:);
      res.U_V = syn.U_V(:);
      res.Umodel_V = Umodel_V(:);
      res.Uelements_V = Uelements_V;
      res.Ulabels = getElementLabels(eecmodel);
      
      % Show fitting results
      displayFittingResults(eecmodel,output);
      
    end % Fit time domain: function res = dcfit(eecmodel,t_s,I_A,U_V,varargin)
    
    function res = ac(eecmodel,f_Hz,varargin)
      %AC analysis
      % Calculate frequency domain response for given frequency vector.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   f_Hz          - Frequency vector in Hertz
      %
      %   Name-Value arguments:
      %   ---------------------
      %
      %   doPlot        - Plot AC response. This can be either a logical
      %                   value or a figure handle which should be used
      %                   instead of a new figure.
      % 
      %   Example:
      %   --------
      %
      %   res = eecmodel.ac(logspace(-2,4,1001),"doPlot",true);

      p = inputParser;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      p.addRequired("f_Hz",@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative"}));
      
      p.addParameter("doPlot",false,@(x)validateattributes(x,["logical" "handle"],{"scalar"}));
      
      p.parse(eecmodel,f_Hz,varargin{:});
      
      % copy valid input frequencies (f > 0)
      f_Hz = p.Results.f_Hz(p.Results.f_Hz > 0);
      omega = 1i*2*pi .* f_Hz(:);
      
      % calculate real/imaginary part given EEC model values
      [Z_Ohm, Zelements_Ohm] = calcAcResponse(eecmodel,omega,[eecmodel.Elements.value],varargin{:});
      
      % copy frequency, real and imaginary part to result
      res.f_Hz = f_Hz(:);
      res.ZRe_Ohm = Z_Ohm(:,1);
      res.ZIm_Ohm = Z_Ohm(:,2);
      res.Zelements_Ohm = Zelements_Ohm;
      res.Zlabels = getElementLabels(eecmodel);
      
      % Plot results?
      if p.Results.doPlot ~= false
        eecmodel.plotResponse(res,'handle',p.Results.doPlot);
      end
      
    end % AC response: res = ac(eecmodel,f_Hz,varargin)
    
    function res = acfit(eecmodel,f_Hz,ZRe_Ohm,ZIm_Ohm,varargin)
      %ACFIT
      % Run lsqcurvefit for given frequencies and minimize error to
      % submitted real and imaginary part.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   f_Hz          - Frequency vector in Hertz
      %   ZRe_Ohm       - Real part of referenz values (must match with 
      %                   frequency vector)
      %   ZIm_Ohm       - Imaginary part of referenz values (must match 
      %                   with frequency vector)
      %
      %   Name-Value arguments:
      %   ---------------------
      %
      %   Method                  - Use Matlab or Spice to calculate
      %                             waveform for each iteration: 
      %                             'analytical' (default) or 'spice'
      %   Repetitions             - Number of repetitions of lsqcurvefit 
      %                             command. The optimizer iterations are
      %                             set separately. Default value is 4
      %
      %   Algorithm               - Which optimizer algorithm should be
      %                             used:
      %                             'trust-region-reflective' (default) or 'levenberg-marquardt'
      %   LowerBoundary           - Absolute values for lower boundary
      %   LowerBoundaryRelative   - Relative values for lower boundary, the
      %                             reference is Element.values
      %   UpperBoundary           - Absolute values for upper boundary
      %   UpperBoundaryRelative   - Relative values for upper boundary, the
      %                             reference is Element.values
      %
      %   SetInitials             - Use current values (starting point) and
      %                             overwrite initial values
      %   UseInitials             - Use initial values as starting point
      %   Display                 - Display options of lsqcurvefit command,
      %                             default is 'off', use 'iter' or 'final'
      %                             for additional information
      %   
      %   Example:
      %   --------
      %
      %   res = eecmodel.acfit(Freq_Hz,ZRe_Ohm,ZIm_Ohm)
      
      p = inputParser;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      p.addRequired("f_Hz",@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative"}));
      p.addRequired("ZRe_Ohm",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(f_Hz)}));
      p.addRequired("ZIm_Ohm",@(x)validateattributes(x,["numeric"],{"nonempty" "size" size(f_Hz)}));
      
      p.addParameter("Method",'analytical',@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      p.addParameter("Repetitions",4,@(x)validateattributes(x,["numeric"],{"nonempty" "nonzero" "nonnegative"}));
      
      p.addParameter("Algorithm",'trust-region-reflective',@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      p.addParameter("LowerBoundary",[eecmodel.Elements.min],@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.min])}));
      p.addParameter("LowerBoundaryRelative",ones(size([eecmodel.Elements.min])),@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.min])}));
      p.addParameter("UpperBoundary",[eecmodel.Elements.max],@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.max])}));
      p.addParameter("UpperBoundaryRelative",ones(size([eecmodel.Elements.max])),@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative" "size" size([eecmodel.Elements.max])}));
      
      p.addParameter("SetInitials",false,@(x)validateattributes(x,["logical"],{"scalar"}));
      p.addParameter("UseInitials",false,@(x)validateattributes(x,["logical"],{"scalar"}));
      
      p.addParameter("Display",'off',@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      
      p.parse(eecmodel,f_Hz,ZRe_Ohm,ZIm_Ohm,varargin{:});
    
      % copy input values and skip invalid frequencies
      f_Hz = p.Results.f_Hz(:);
      ZRe_Ohm = p.Results.ZRe_Ohm(f_Hz > 0);
      ZIm_Ohm = p.Results.ZIm_Ohm(f_Hz > 0);
      f_Hz = f_Hz(f_Hz > 0);

      % Reset to initial values
      if logical(p.Results.SetInitials)
        setElementsValues(eecmodel,[eecmodel.Elements.value],'initial');
      end
      
      % Use initials for first iteration
      if logical(p.Results.UseInitials)
        setElementsValues(eecmodel,[eecmodel.Elements.initial],'value');
      end
      
      % Get algorithm and boundaries if needed
      [algorithm, lb, ub] = fittingBoundaries(eecmodel,varargin{:});
      
      % Set lsqcurvefit options
      options = optimoptions(@lsqcurvefit,...
        'Algorithm', algorithm,...
        'TolX', 1e-20, ...
        'TolFun', 1e-20,...
        'Display', p.Results.Display, ...
        'MaxIter', 1e3, ...
        'MaxFunEvals', 20e3, ...
        'FunValCheck', 'on',...
        'UseParallel', false);
      
      tic;
      for i = 1:p.Results.Repetitions
        % split input data in real and imaginary parts
        omega = 1i*2*pi .* f_Hz;
        Z = [ZRe_Ohm(:) ZIm_Ohm(:)];
        
        % check for selected method {'spice' 'analytical'}
        switch lower(p.Results.Method)
          case {'spice'}
            % calculate impedance using ac analysis, NOT implemented yet
            ME = MException('eec:acfit:UseSpice', ...
              'Fitting methode "spice" is not implemented yet, please use "analytical" (default) instead.');
            throw(ME);
            
          otherwise
            % calculate impedance using transfer function
            fun = @(values,omega)calcAcResponse(eecmodel,omega,values);
            
        end
        
        % lsqcurvefit: Fit impedance in the frequency domain
        [fittedParameters, ~, ~, exitflag, output] = lsqcurvefit(fun,...
          [eecmodel.Elements.value],omega,Z,...
          lb,ub,options);
        output.Exitflag = exitflag;
      end
      output.ComputationTime = toc;
      
      % set fitted values
      setElementsValues(eecmodel,fittedParameters,'value');
      
      % sort RC elements time constants in ascending order
      sortElements(eecmodel);
      
      % calculate real/imaginary part for fitted parameters
      Z = fun([eecmodel.Elements.value],omega);
      
      % copy frequency, real and imaginary part to result
      res.f_Hz = f_Hz;
      res.ZRe_Ohm = Z(:,1);
      res.ZIm_Ohm = Z(:,2);
      
      % Show fitting results
      displayFittingResults(eecmodel,output);
      
    end % Fit frequency domain: function res = acfit(eecmodel,f_Hz,ZRe_Ohm,ZIm_Ohm,varargin) 
    
    function res = spice(eecmodel,t_s,I_A,varargin)
      %SPICE analysis
      % Run spice simulation for given current waveform. The range of 
      % functions is very limited at the moment. Only simple .tran (DC
      % response) command is implemented.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   t_s           - Time vector in seconds (or duration)
      %   I_A           - Current vector (must match with time vector) or
      %                   constant current value in Ampere
      %
      %   Name-Value arguments:
      %   ---------------------
      %
      %   SpicePath             - Path to LTspice XVII executable to run spice
      %                           simulation in batch mode.
      %   MaxStepsize           - Max stepsize parameter [dTmax] of .tran
      %                           command
      %   UseInitialConditions  - Set .uic flag of .tran command
      %   doPlot                - Plot DC response. This can be either a logical
      %                           value or a figure handle which should be used
      %                           instead of a new figure.
      % 
      %   Example:
      %   --------
      %
      %     Constant current value:
      %   res = eecmodel.spice(linspace(0,10,101),-2,"doPlot",true);
      %     or current waveform:
      %   res = eecmodel.spice(linspace(0,10,101),sin(2*pi*0:100));
      
      p = inputParser;
      
      p.addRequired("eecmodel",@(x)validateattributes(x,"eec",{"scalar"}));
      p.addRequired("t_s",@(x)validateattributes(x,["numeric" "duration"],{"nonempty"}));
      p.addRequired("I_A",@(x)validateattributes(x,["numeric"],{"nonempty"}));
      
      p.addParameter("SpicePath","",@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      p.addParameter("MaxStepsize",1e-3,@(x)validateattributes(x,["numeric"],{"scalar"}));
      p.addParameter("UseInitialConditions",false,@(x)validateattributes(x,["logical"],{"scalar"}));
      
      p.addParameter("doPlot",false,@(x)validateattributes(x,["logical" "handle"],{"scalar"}));
      
      p.parse(eecmodel,t_s,I_A,varargin{:});
      
      % LTspice XVII path
      if (p.Results.SpicePath ~= "")
        eecmodel.SpicePath = char(p.Results.SpicePath);
      end
      
      % Check for valid LTspice executable
      if ~exist(eecmodel.SpicePath, 'file')
        if ismac
            % default path on MAC
            defname = "/Applications/LTspice.app/Contents/MacOS/LTspice";
        elseif isunix
            % default path on Linux
        elseif ispc
            % default path on Windows
            defname = "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe";
        end
        [file,path] = uigetfile('*.*','Select LTspice executable',defname);
        eecmodel.SpicePath = fullfile(path,file);
      end
      
      % convert duration to seconds
      if isduration(p.Results.t_s)
        t_s = seconds(p.Results.t_s(:));
      else
        t_s = p.Results.t_s(:);
      end
      
      % check for continuous or contant value
      if isscalar(p.Results.I_A)
        I_A = ones(size(t_s)) .* p.Results.I_A;
      else
        if ~isequal(size(p.Results.I_A(:)), size(t_s))
          error('eec:eec:IncorrectDataSize', ...
            'Invalid vector length, "t_s" [%s] and "I_A" [%s] must match in size.', ...
            eec.sizeDisp(t_s), eec.sizeDisp(p.Results.I_A(:)));
        else
          I_A = p.Results.I_A(:);
        end
      end
      I_A(isnan(I_A)) = 0;
      
      % run spice simulation and handle results
      [U_V, sim] = runSpiceNetlist(eecmodel,t_s,I_A,varargin{:});
      
      % copy duration, voltage and current to result
      res.t_s = seconds(t_s);
      res.Umodel_V = U_V;
      res.I_A = I_A;
      
      % copy simulation results to class
      eecmodel.SpiceSimulation = sim;
      
      % get detailed spice results and time, voltage and current for each element
      Uelements_V = nan(length(sim.time_vect),sum(contains(sim.variable_type_list, 'voltage'))+1);
      UelementsName = cell(1,size(Uelements_V,2));
      Uelements_V(:,1) = sim.variable_mat(strcmpi(sim.variable_name_list, 'V(v001)'),:)';
      UelementsName{1} = 'U(bat)';
      cnt = 2;
      for i = 1:length(eecmodel.Elements)
        UelementsName{cnt} = sprintf('U(%s%d)',eecmodel.Elements([eecmodel.Elements.Order] == i).Type,i);
        np = sim.variable_mat(strcmpi(sim.variable_name_list, sprintf('V(v%03d)',i)),:)';
        nn = sim.variable_mat(strcmpi(sim.variable_name_list, sprintf('V(v%03d)',i+1)),:)';
        % last element is connected to negative netname GND which is zero at all
        % time
        if isempty(nn)
          Uelements_V(:,i+1) = np;
        else
          Uelements_V(:,i+1) = np-nn;
        end
        cnt = cnt +1;
      end
      
      Ielements_A = nan(length(sim.time_vect),sum(contains(sim.variable_type_list, 'current')));
      IelementsName = cell(1,size(Ielements_A,2));
      Ielements_A(:,1) = sim.variable_mat(strcmpi(sim.variable_name_list, 'I(Ibat)'),:)';
      IelementsName{1} = 'I(bat)';
      cnt = 2;
      for i = 1:length(eecmodel.Elements)
        % some elements might have more than one current
        idx = find(cellfun(@(x)~isempty(x), regexp(sim.variable_name_list, ['^I(.+' num2str(i) ')'], 'start')));
        for j = 1:length(idx)
          IelementsName{cnt} = sim.variable_name_list{idx(j)};
          Ielements_A(:,cnt) = sim.variable_mat(idx(j),:)';
          cnt = cnt +1;
        end
      end
      
      % detailed spice results
      res.tspice_s = seconds(sim.time_vect(:));
      res.Uelements_V = Uelements_V;
      res.Ulabels = UelementsName;
      res.Ielements_A = Ielements_A;
      res.Ilabels = IelementsName;
      
      % Plot results?
      if p.Results.doPlot ~= false
        eecmodel.plotResponse(res,'handle',p.Results.doPlot);
      end
      
    end % Spice: res = spice(eecmodel,t_s,I_A,varargin)
    
  end % methods(Access = public, Hidden = false)
  
  % Interactive UI wrapper functions:
  % TIMESERIES, NYQUIST
  methods(Access = public, Hidden = false)
    
    function app = timeseries(eecmodel,tI_s,I_A,tU_s,U_V,varargin)
      %TIMESERIES
      % Open interactive UI to investigate EEC time response.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   tI_s          - Current time vector in seconds
      %   I_A           - Current vector (must match with time vector) or
      %                   constant current value in Ampere
      %   tU_s          - Voltage time vector in seconds
      %   U_V           - Voltage vector as observed output vector in Volt
      % 
      %   Example:
      %   --------
      %
      %   app = eecmodel.timeseries(tI_s,I_A,tU_s,U_V);
      
      I_idx = (tI_s >= max(tI_s(1),tU_s(1)) & tI_s <= min(tI_s(end),tU_s(end)));
      U_idx = (tU_s >= max(tI_s(1),tU_s(1)) & tU_s <= min(tI_s(end),tU_s(end)));
      
      % Run MATLAB app as separate object
      app = ShowTimeseries(eecmodel,tU_s(U_idx),U_V(U_idx),tI_s(I_idx),I_A(I_idx),varargin{:});
      
      if nargout == 0
        clear app
      end
      
    end % function app = timeseries(eecmodel,tI_s,I_A,tU_s,U_V,varargin)
    
    function app = nyquist(eecmodel,f_Hz,ZRe_Ohm,ZIm_Ohm,varargin)
      %TIMESERIES
      % Open interactive UI to investigate EEC frequency response.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   f_Hz          - Frequency vector in Hertz
      %   ZRe_Ohm       - Real part of referenz values (must match with 
      %                   frequency vector)
      %   ZIm_Ohm       - Imaginary part of referenz values (must match 
      %                   with frequency vector)
      % 
      %   Example:
      %   --------
      %
      %   app = eecmodel.nyquist(f_Hz,ZRe_Ohm,ZIm_Ohm);
      
      % Run MATLAB app as separate object
      app = ShowNyquist(eecmodel,f_Hz,ZRe_Ohm,ZIm_Ohm,varargin{:});
      
      if nargout == 0
        clear app
      end
      
    end % function app = nyquist(eecmodel,f_Hz,ZRe_Ohm,ZIm_Ohm,varargin)
    
  end
  
  % Constructor/Destructor methods:
  % EEC, DELETE
  methods(Access = public, Hidden = true)
    
    function eecmodel = eec(ModelString,varargin)
      %EEC
      % Constructor class returns EEC model object.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   ModelString    - List of EEC elements, can be separated (, + - _)
      %                    string or cell array containing one or multiple
      %                    elements of 'OCV' 'R' 'C' 'L' 'RC' or 'RL'.
      %
      %   Name-Value arguments:
      %   ---------------------
      %
      %   InitialValues  - Initial values for each EEC element. 'OCV', 'R',
      %                    'C' and 'L' consist of one value (voltage, resistance,
      %                    capacity or inductivity) while 'RC' and 'RL' have
      %                    three (resistance, time constant, initial voltage/charge)
      %                    values.
      %   SpicePath      - Path to LTspice XVII executable to run spice
      %                    simulation in batch mode.
      %
      % eecmodel = eec("R+RC");
      % eecmodel = eec("OCV,R,L,RC,RC,RC",...
      %   ...%             OCV     R       L        RC1             RC2             RC3
      %   'InitialValues',[3.6500  0.0300  500e-9   5.0e-3 1.0e-3 0 10e-3 100e-3 0  20e-3 10 0 ]);
      
      if nargin == 0
        return;
      end
      
      p = inputParser;
      
      p.addRequired("ModelString",@(x)validateattributes(x,["string" "char" "cell"],{"nonempty" "scalartext"}));
      p.addParameter("InitialValues",[],@(x)validateattributes(x,["numeric"],{"nonempty" "nonnegative"}));
      p.addParameter("SpicePath","",@(x)validateattributes(x,["string" "char"],{"scalartext"}));
      
      p.parse(ModelString,varargin{:});
      
      ModelString = p.Results.ModelString;
      
      if isempty(ModelString)
        ME = MException('eec:eec:IncorrectDataType', ...
          'ModelString is empty');
        throw(ME);
      end
      
      switch class(ModelString)
        
        case {'string','char'}
          ModelString = strsplit(ModelString,{',' '+' '-' '_'});
          if any(cellfun(@isempty, ModelString) == 1)
            error(message('eec:eec:IncorrectDataType','ModelString'));
          end
          
        case {'cell'}
          if ~all(cellfun(@ischar, ModelString) == 1) || any(cellfun(@isempty, ModelString) == 1)
            error(message('eec:eec:IncorrectDataType','ModelString'));
          end
          ModelString = cellstr(ModelString);
          
      end
      
      if any(cellfun(@(val)sum(strcmpi(val, eec.ValidElements)), ModelString) ~= 1)
        error('eec:eec:IncorrectDataValue', ...
          'Invalid model string, valid EEC elements are "%s".', strjoin(eec.ValidElements, '", "'));
      end
      
      % LTspice XVII path
      SpicePath = char(p.Results.SpicePath);
      if ~isempty(SpicePath)
        eecmodel.SpicePath = SpicePath;
      end
      
      % get initial parameters for given EEC elements
      cnt = 1;
      for i = 1:length(ModelString)
        
        element.Order = i;
        element.Type = upper(ModelString{i});
        
        switch element.Type
          case 'OCV'
            % not visible in frequency domain
            element.TransferFunction = @(~,omega)zeros(size(omega));
            
            % time domain: U = OCV(SOC)
            element.Equation{1} = 'U = OCV(SOC)';
            element.TimeResponse = eval(sprintf('@(values,t_s,I_A)values(%d) .* ones(size(t_s))',cnt));
            
            element.Description{1} = sprintf('Voltage of OCV%d [V]',i);
            element.initial(1) = 3.65;
            
          case 'R'
            % frequency domain: tf = R
            element.TransferFunction = eval(sprintf('@(values,omega)ones(size(omega)).*values(%d)',cnt));
            
            % time domain: U = R * I
            element.Equation{1} = 'U = R * I';
            element.TimeResponse = eval(sprintf('@(values,t_s,I_A)values(%d) .* I_A',cnt));
            
            element.Description{1} = sprintf('Ohmic resistance R%d [Ohm]',i);
            element.initial(1) = 35e-3;
            
          case 'L'
            % frequency domain: tf = s*L
            element.TransferFunction = eval(sprintf('@(values,omega) omega .* values(%d)',cnt));
            
            % time domain: U = L * dI/dt
            element.Equation{1} = 'U = L * dI/dt';
            element.TimeResponse = eval(sprintf('@(values,t_s,I_A)values(%d).*smooth(smooth([0;diff(I_A)]./[1;diff(t_s)],2e-6/min(diff(t_s))),2e-6/min(diff(t_s)))',cnt));
            
            element.Description{1} = sprintf('External Inductivity L%d [H]',i);
            element.initial(1) = 15e-9;
            
          case 'CL'
            % frequency domain: tf = s*L
            element.TransferFunction = eval(sprintf('@(values,omega) omega .* values(%d)',cnt));
            
            % time domain: U = L * dI/dt
            element.Equation{1} = 'U = L * dI/dt';
            element.TimeResponse = eval(sprintf('@(values,t_s,I_A)values(%d).*smooth(smooth([0;diff(I_A)]./[1;diff(t_s)],2e-6/min(diff(t_s))),2e-6/min(diff(t_s)))',cnt));
            
            element.Description{1} = sprintf('External Inductivity L%d [H]',i);
            element.initial(1) = 15e-9;
            
          case 'C'
            % frequency domain: tf = 1/(s*C)
            element.TransferFunction = eval(sprintf('@(values,omega) 1 ./ (omega .* values(%d))',cnt));
            
            % time domain: U = Q/C = (\int I dt)/C
            element.Equation{1} = 'U = Q / C';
            element.TimeResponse = eval(sprintf('@(values,t_s,I_A)cumtrapz(t_s, I_A) ./ values(%d)',cnt));
            
            element.Description{1} = sprintf('Capacitor C%d [F]',i);
            element.initial(1) = 100;
            
          case 'RC'
            % frequency domain: tf = R / (1 + s*\tau)
            element.TransferFunction = eval(sprintf('@(values,omega) values(%d) ./ (1 + omega .* values(%d) )',cnt,cnt+1));
            
            %           rising: U = U_max * (1 - exp( -t/\tau )) = (R * I) * (1 - exp( -(t+dt)/RC ))
            element.Equation{1} = 'U = (R * I) * (1 - exp( -(t+dt) / \tau )) where \tau = R*C';
            
            % inflection point: dt = -log((1 - exp( -(dt_rising) / \tau ))) * \tau == -log((1 - exp( -(dt_rising) / RC ))) * RC
            %            pause: U = U_end * exp( -(t+dt)/\tau )) = U_end * exp( -(t+dt)/RC )
            element.Equation{2} = 'U = (R * I) * exp( -(t+dt) / \tau )) where dt = -log( 1 - exp( -(t_end)/\tau) ) * \tau';
            element.TimeResponse = eval(sprintf('@(values,t_s,I_A)eec.rcNumerical(values(%d:%d),t_s,I_A)',cnt,cnt+2));
            
            element.Description{1} = sprintf('Resistor of RC%d element [Ohm]',i);
            element.initial(1) = 5e-3;
            
            element.Description{2} = sprintf('Time constant of RC%d element [s]',i);
            element.initial(1,2) = 1;
            
            element.Description{3} = sprintf('Initial voltage of RC%d element [V]',i);
            element.initial(1,3) = 0;
            
          case 'RL'
            % frequency domain: tf = sR / (s - 1/\tau)
            element.TransferFunction = eval(sprintf('@(values,omega) values(%d) .* omega ./ (omega - 1 / values(%d) )',cnt,cnt+1));
            
            %           rising: U = L * dI/dt * (1 - exp( -t/\tau )) = L * dI/dt * (1 - exp( -t/(L/R) ))
            element.Equation{1} = 'U = (L * dI/dt) * (1 - exp( -t/\tau )) where \tau = L/R';
            
            % inflection point: dt = -log((1 - exp( -(dt_rising) / \tau ))) * \tau == -log((1 - exp( -(dt_rising) / (L/R) ))) * (L/R)
            %            pause: U = L * dI/dt * exp( -(t+dt)/\tau )) = L * dI/dt * exp( -(t+dt)/(L/R) )
            element.Equation{2} = 'U = (L * dI/dt) * exp( -(t+dt)/\tau )) where dt = -log(1-exp(-(t_end)/\tau))*\tau';
            element.TimeResponse = eval(sprintf('@(values,t_s,I_A)eec.rlNumerical(values(%d:%d),t_s,I_A)',cnt,cnt+2));
            
            element.Description{1} = sprintf('Resistor of RL%d element [Ohm]',i);
            element.initial(1) = 5e-3;
            
            element.Description{2} = sprintf('Time constant of RL%d element [s]',i);
            element.initial(1,2) = 10e-6;
            
            element.Description{3} = sprintf('Initial current of RL%d element [A]',i);
            element.initial(1,3) = 0;
            
        end
        
        % append element
        element.value = zeros(size(element.initial));
        element.min = zeros(size(element.initial));
        element.max = inf(size(element.initial));
        eecmodel.Elements(end+1) = element;
        
        % adjust cnt variable to match [eecmodel.Elements.value] order 
        cnt = cnt +length(element.value);
        
      end % for i = 1:length(ModelString)
      
      % check initial parameters input
      if ~isempty(p.Results.InitialValues)
        % initial parameters must match in size, otherwise keep auto values
        if ~isequal(length(p.Results.InitialValues(:)), length([eecmodel.Elements.initial]))
          warning('Initial parameters [%d] must match the EEC element parameters length [%d], using default values instead.',...
            length(p.Results.InitialValues(:)),length([eecmodel.Elements.initial]));
        else
          % convert input parameters to cell array
          definput = arrayfun(@(x)sprintf('%g',x), p.Results.InitialValues, 'UniformOutput', false);
        end
      else
        definput = arrayfun(@(x)sprintf('%g',x), [eecmodel.Elements.initial], 'UniformOutput', false);
      end
      
      % set initial values using input dialog
      prompt = [eecmodel.Elements.Description];
      answer = [];
      title = 'Initial values';
      dims = [1 50];
      opts.Interpreter = 'tex';
      while isempty(answer)
        
        if any(contains(p.UsingDefaults, 'InitialValues'))
          % show input dialog to define initial values
          answer = inputdlg(prompt,title,dims,definput,opts);
        else
          % use input parameters without dialog
          answer = definput;
        end
        
        % check for cancelation
        if isempty(answer)
          % remove all elements
          eecmodel.Elements(:) = [];
          return;
        end
        
        try
          % try to parse answers
          tail = 1;
          for i = 1:length(ModelString)
            idx = [eecmodel.Elements.Order] == i;
            head = tail +( length(eecmodel.Elements(idx).initial) - 1);
            val = str2double(answer(tail:head));
            eecmodel.Elements(idx).initial = val(:)';
            eecmodel.Elements(idx).value = val(:)';
            tail = head +1;
          end
          
        catch e
          % otherwise show warning and reopen input dialog
          warning(eec.extractExceptionMessage(e))
          definput = answer;
          answer = [];
        end
        
      end
      
      % sort RC elements time constants in ascending order
      sortElements(eecmodel);
      
      % set model string for display
      eecmodel.ModelString = ModelString;
      
    end % function eecmodel = eec(ModelString,varargin)
    
    function delete(eecmodel)
      %DELETE
      % Destructor class deletes EEC model instance.
      
      if isvalid(eecmodel)
        close(eecmodel);
      end
      
    end % function delete(obj)
    
  end % methods(Access = public, Hidden = true)
  
  % Display class methods:
  % DISPLAYSCALAROBJECT, DISPLAYFITTINGRESULTS
  methods(Access = protected)
    
    function displayScalarObject(eecmodel)
      %DISPLAYSCALAROBJECT
      % Default command window output whenever an EEC instance is
      % displayed.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      
      header = matlab.mixin.CustomDisplay.getSimpleHeader(eecmodel);
      disp(header);
      
      % Default Properties
      disp(['    ', '          SpicePath: ''', char(eecmodel.SpicePath), '''']);
      disp(['    ', '        ModelString: {''', char(strjoin(string(eecmodel.ModelString), "', '")), '''}']);
      fprintf('\n');
      
      % EEC elements
      for i = 1:length(eecmodel.Elements)
        
        element = eecmodel.Elements([eecmodel.Elements.Order] == i);
        
        % show elements as table
        fprintf("+-------------------------------------------------------------------------------+\n");
        fprintf("| <strong>%2d. element %s</strong> %s |\n",i,element.Type,repmat(' ',1,64-numel(element.Type)));
        switch upper(element.Type)
          case {'R', 'C', 'L'}
            fprintf("| continuous: %s %s |\n",element.Equation{1},repmat(' ',1,64-numel(element.Equation{1})));
          case {'RC', 'RL'}
            fprintf("|     rising: %s %s |\n",element.Equation{1}(1:42),repmat(' ',1,64-numel(element.Equation{1}(1:42))));
            fprintf("|             %s %s |\n",element.Equation{1}(43:end),repmat(' ',1,64-numel(element.Equation{1}(43:end))));
            fprintf("|    falling: %s %s |\n",element.Equation{2}(1:37),repmat(' ',1,64-numel(element.Equation{2}(1:37))));
            fprintf("|             %s %s |\n",element.Equation{2}(38:end),repmat(' ',1,64-numel(element.Equation{2}(38:end))));
        end
        fprintf("|%s|\n", repmat(' ',1,79));
        for j = 1:length(element.Description)
          str = sprintf(['%' num2str(max(cellfun(@numel,element.Description))) 's: %g'],element.Description{j},element.value(j));
          fprintf("| %s %s |\n",str,repmat(' ',1,76-numel(str)));
        end
      end
      fprintf("+-------------------------------------------------------------------------------+\n");      
      fprintf('\n');
      
    end % function displayScalarObject(eecmodel)
    
    function displayFittingResults(eecmodel,output)
      %DISPLAYFITTINGRESULTS
      % Show fitting results for each EEC element formatted as table. Also
      % display dense information about lsqcurvefit optimization process.
      %
      %   Input Arguments:
      %   ----------------
      %
      %   eecmodel      - EEC model object
      %   output        - lsqcurvefit optimization process output structure
      
      % Get detailed exitflag information
      switch output.Exitflag
        case 1
          output.Status = 'Function converged to a solution x.';
        case 2
          output.Status = 'Change in x was less than the specified tolerance.';
        case 3
          output.Status = 'Change in the residual was less than the specified tolerance.';
        case 4
          output.Status = 'Magnitude of search direction was smaller than the specified tolerance.';
        case 0
          output.Status = 'Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations.';
        case -1
          output.Status = 'Output function terminated the algorithm.';
        case -2
          output.Status = 'Problem is infeasible: the bounds lb and ub are inconsistent.';
      end
      
      % show results
      if output.Exitflag > 0
        fprintf("\nFitting ('%s') <strong>successful</strong> after %d iterations, computation time %fs\n%s\n",output.algorithm,output.iterations,output.ComputationTime,output.Status);
      else
        fprintf("\nFitting ('%s') <strong>failed</strong> after %d iterations, computation time %fs\n%s\n",output.algorithm,output.iterations,output.ComputationTime,output.Status);
      end
      % show elements as table
      fprintf("\n+-------------------------------------------------------------------------------+\n");
      fprintf("|                                     Initial value    Final value       Diff.  |\n");
      fprintf("+-------------------------------------------------------------------------------+\n");
      n = max(cellfun(@numel,[eecmodel.Elements.Description]));
      for i = 1:length(eecmodel.Elements)
        element = eecmodel.Elements([eecmodel.Elements.Order] == i);
        fprintf("| <strong>%2d. element %- 65s</strong> |\n",i,element.Type);
        for j = 1:length(element.initial)
          boundary = '';
          if element.value(j) <= element.min(j)*1.0001
            boundary = 'lb';
          elseif element.value(j) >= element.max(j)*0.9999
            boundary = 'ub';
          end
          str = sprintf(['%' num2str(n) 's: %12g  %12g %2s  %+7.1f %%'],element.Description{j},element.initial(j),element.value(j),boundary,(element.value(j)-element.initial(j))/element.initial(j)*100);
          fprintf("| %-77s |\n",str);
          %           fprintf("| %s                  %e              %e         %+7.2f %%   |\n",element.initial(j),element.value(j),(element.value(j)-element.initial(j))/element.initial(j)*100);
        end
      end
      fprintf("+-------------------------------------------------------------------------------+\n");
      fprintf('\n');
      
    end % function displayFittingResults(eecmodel,output)
    
  end % methods(Access = protected)
  
end