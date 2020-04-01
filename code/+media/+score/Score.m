classdef Score < media.Graphic
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
        X = 20;
        Y = 10% + 144;
        COLOR = [255 255 255];
        FONTSIZE = 144;
    end

    properties (GetAccess = private, SetAccess = private)
    end
    
    
    %% ====================================================================
    %  -----CONSTRUCTOR/DESTRUCTOR METHODS-----
    %$ ====================================================================    

    methods (Access = public)
        
        function obj = Score()
        end
        
    end
    
    
    %% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
          
    methods (Access = public)
        
        function [] = update(obj, winhandle)
            % formatting
            Screen('TextFont', winhandle, 'Courier New');
            Screen('TextSize', winhandle, obj.FONTSIZE);
            Screen('TextStyle', winhandle, 1+2);
            
            % compute text
            txt = sprintf('%i', DataManager.getInstance().score);
            
            % draw:
            % NB: formatted text only works for windows displays with
            % scaling >100% IF gstreamer 1.4 installed
            % DrawFormattedText(winhandle, txt, obj.X, obj.Y, obj.COLOR, [],[],[], 3);
            % safer:
            Screen('DrawText', winhandle, txt, obj.X, obj.Y, obj.COLOR);
        end
        
    end
    

    %% ====================================================================
    %  -----STATIC METHODS-----
    %$ ====================================================================    
    
    methods (Static, Access = public)
        
    end
    
end