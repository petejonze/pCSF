classdef Feedback < graphic
    % dfdfd
    %
    %   dfdfdf
    %
    % See Also:
    %   none
    %
    % Example:
    %   none
    %
    % Author:
    %   Pete R Jones <petejonze@gmail.com>
    %
    % Verinfo:
    %   1.0 PJ 19/04/2018 : first_build\n
    %
    % Todo:
    %   none
    %
    % Copyright 2018 : P R Jones
    % *********************************************************************
    %    

  	%% ====================================================================
    %  -----PROPERTIES-----
    %$ ====================================================================
                      
    properties (Constant)
        IMG = '/coins/gold_coin.png';
    end

    properties (GetAccess = public, SetAccess = private)
        graphic
    end
    
    %% ====================================================================
    %  -----CONSTRUCTOR/DESTRUCTOR METHODS-----
    %$ ====================================================================    

    methods (Access = public)
        
        function obj = Lives(winhandle, imgDir)
            % initialise graphic(s)
            obj.graphic = PtrVisual(fullfile(imgDir, obj.IMG), winhandle, 1);
            
            % register with graphics manager
            dfdfdf
        end
        
    end
    
    
    %% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
          
    methods (Access = public)
        
    end
    

    %% ====================================================================
    %  -----STATIC METHODS-----
    %$ ====================================================================    
    
    methods (Static, Access = public)
        
    end
    
end