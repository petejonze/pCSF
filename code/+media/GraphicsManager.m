classdef (Sealed) GraphicsManager < Singleton
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
    
    properties (GetAccess = public, SetAccess = private)
        winhandle
    end
        
    properties (GetAccess = private, SetAccess = private)
        graphicObjs
    end
    
    %% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
    
    methods (Access = public)
        
        %% == CONSTRUCTOR =================================================
        
        function obj = GraphicsManager(winhandle)
            % GraphicsManager Constructor.
            %
            % @date     19/04/18
            % @author   PRJ
            %       
            
            % init
            obj.graphicObjs = {};
            
            % store
            obj.winhandle = winhandle;
        end
        
        function [] = delete(obj)
            % GraphicsManager Destructor.
            %
            % @date     19/04/18
            % @author   PRJ
            %       
            
            % release buffer:
            for i = 1:length(obj.graphicObjs)
                delete(obj.graphicObjs{i});
                obj.graphicObjs{i} = [];
            end
        end

        %% == METHODS =====================================================
        
        function [] = register(obj, gfx)
            obj.graphicObjs{end+1} = gfx;
        end
        
        function [] = remove(obj, gfx)
            % search for item
            for i = length(obj.graphicObjs):-1:1 % start from end (usually FIFO)
                if obj.graphicObjs{i} == gfx
                    obj.graphicObjs(i) = [];
                    return; % stop immediately if object has already been found
                end
            end
            
            % defensive
            error('Requested object could not be found?? Nothing removed');
        end
        
        function [] = drawAll(obj)
            gfx = obj.graphicObjs; % in case some get deleted as part of update
            for i = 1:length(gfx)
                gfx{i}.update(obj.winhandle);
            end
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
        function obj = getInstance()
            obj = Singleton.getInstanceSingleton(mfilename('class'));
        end
        function [] = finishUp()
            Singleton.finishUpSingleton(mfilename('class'));
        end
    end
    
end