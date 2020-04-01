classdef (Sealed) DataManager < Singleton
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
    %   1.0 PJ 20/04/2018 : first_build\n
    %   1.1 PJ 26/11/2019 : tweak to prevent score < 0\n
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
    
    properties (GetAccess = public, SetAccess = private)
        score = 0;
    end
    
    %% ====================================================================
    %  -----SINGLETON CONSTRUCTOR-----
    %$ ====================================================================    

    methods (Access = ?Singleton)

        function obj = DataManager()
            % DataManager Constructor.
            %
            % @date     20/04/18
            % @author   PRJ
            %       
        end

    end
    
    
  	%% ====================================================================
    %  -----STATIC METHODS-----
    %$ ====================================================================

    methods (Static, Access = public)   
        
        function [] = addPoints(nPoints)
            if nargin < 1 || isempty(nPoints)
                nPoints = 1;
            end
            
            obj = DataManager.getInstance();
            obj.score = max(0, obj.score + nPoints); % cannot go below 0
        end
        
    end
    
    
  	%% ====================================================================
    %  -----SINGLETON BLURB-----
    %$ ====================================================================

    methods (Static, Access = ?Singleton)
        function obj = getSetSingleton(obj)
            persistent singleObj
            if nargin > 0, singleObj = obj; end
            obj = singleObj;
        end
    end
    methods (Static, Access = public)
        function obj = init(varargin)
            obj = Singleton.initialise(mfilename('class'), varargin{:});
        end
        function obj = getInstance()
            obj = Singleton.getInstanceSingleton(mfilename('class'));
        end
        function [] = finishUp()
            Singleton.finishUpSingleton(mfilename('class'));
        end
    end
    
end