classdef Feedback < media.Graphic
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
        DURATION_SECS = 1;
        internalResources = media.feedback.FeedbackMedia();
    end

    properties (GetAccess = private, SetAccess = private)
        img
        tStart_secs
    end
    
    
    %% ====================================================================
    %  -----CONSTRUCTOR/DESTRUCTOR METHODS-----
    %$ ====================================================================    

    methods (Access = public)
        
        function obj = Feedback(xy_px, wasHit)
        	fprintf('Feedback.m :: Feedback starting\n');
                        
            % initialise timer
            obj.tStart_secs = GetSecs();
            
            % place graphic
            if wasHit
                obj.img = copy(obj.internalResources.img_hit);
            else
                obj.img = copy(obj.internalResources.img_miss);
            end
          	obj.img.setXY(xy_px(1), xy_px(2), true);
            
                
            % play audio
            if wasHit
                obj.internalResources.snd_hit.play();
            else
                obj.internalResources.snd_miss.play();
            end
        end
        
    end
    
    
    %% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
          
    methods (Access = public)
        
        function [] = update(obj, winhandle)
            obj.img.draw();
            
            if (GetSecs() - obj.tStart_secs) > obj.DURATION_SECS
                media.GraphicsManager.getInstance().remove(obj);
                delete(obj);
            end
            
        end
        
    end
    

    %% ====================================================================
    %  -----STATIC METHODS-----
    %$ ====================================================================    
    
    methods (Static, Access = public)
        
    end
    
end