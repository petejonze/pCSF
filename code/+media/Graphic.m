classdef (Abstract) Graphic < handle
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
       
    properties (GetAccess = public, SetAccess = protected)
    end

    
    %% ====================================================================
    %  -----ABSTRACT PUBLIC METHODS-----
    %$ ====================================================================
    
    methods(Abstract, Access = public)

        % Draw.
        %
        % @date     19/04/18
        % @author   PRJ
        %
        update(obj, winhandle)
    end
    
    
    %% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
      
    methods(Access = public)

        %% == METHODS =====================================================

        
    end
    
	%% ====================================================================
    %  -----PROTECTED METHODS-----
    %$ ====================================================================
      
    methods(Access = protected)

        function obj = Graphic()
            % register with graphics manager
            media.GraphicsManager.getInstance().register(obj);
        end
        
    end
end