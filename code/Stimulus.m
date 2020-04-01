classdef Stimulus < handle
    % dfdfd
    %
    %   dfdfdf
    %
    % Matlab:           v2012(?)
    %
    % See also:         none
    %
    % Example:          pop_CSF_v0_0_1()
    %
    % Author(s):    	Pete R Jones <petejonze@gmail.com>
    %
    % Version History:  0.0.1	PJ  18/04/2018    Initial build
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
          
    properties (Constant) %@TODO: make variable and set at start
        ONRAMP_SECS = 1;
        STIM_DURATION_SECS = 6;
        OFFRAMP_SECS = 1;
        POST_STIM_GRACE_SECS = 1;
    end
    
    properties (GetAccess = public, SetAccess = private)
        % gabor params
        phase
        freq_cpd
        freq_cpp
        sc_px
        contrast
        aspectratio
        estHPoseAtOnset
        %
        log_x_px
        log_y_px
        log_t_secs
        % time
        tStart_secs
    end
    
    
    %% ====================================================================
    %  -----CONSTRUCTOR/DESTRUCTOR METHODS-----
    %$ ====================================================================    

    methods (Access = public)
        
        function obj = Stimulus(phase, freq_cpd, freq_cpp, sc_px, contrast, aspectratio, estHPoseAtOnset)
            % init
            obj.tStart_secs = GetSecs();
            
            % store
            obj.phase       = phase;
            obj.freq_cpd    = freq_cpd;
            obj.freq_cpp    = freq_cpp;
            obj.sc_px       = sc_px;
            obj.contrast    = contrast;
            obj.aspectratio = aspectratio;
            obj.estHPoseAtOnset = estHPoseAtOnset;
        end
    end
    
    
    %% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
          
    methods (Access = public)
        
        function p = getPhase(obj)
            p = obj.phase;
        end
        
        function f_cpd = getFreq_cpd(obj)
            f_cpd = obj.freq_cpd;
        end
        
        function f_cpp = getFreq_cpp(obj)
            f_cpp = obj.freq_cpp;
        end
        
        function sc = getSpatialConstant(obj)
            sc = obj.sc_px;
        end
                    
        function c = getContrast(obj)
        	c = obj.getWindow * obj.contrast;
        end
        
        function ar = getAspectRatio(obj)
            ar = obj.aspectratio;
        end
        
        function [] = logLocation(obj, x_px, y_px, t_secs)
            obj.log_x_px(end+1) = x_px;
            obj.log_y_px(end+1) = y_px;
            obj.log_t_secs(end+1) = t_secs;
        end
        function xyt = getLocationLog(obj)
            xyt = [obj.log_x_px; obj.log_y_px; obj.log_t_secs];
        end
        
        function istimedout = checkForTimeout(obj)
            istimedout = (GetSecs() - obj.tStart_secs) > (obj.STIM_DURATION_SECS + obj.POST_STIM_GRACE_SECS);
        end
        
    end

    
    %% ====================================================================
    %  -----PRIVATE METHODS-----
    %$ ====================================================================
          
    methods (Access = private)

        function w = getWindow(obj)
            tElapsed_secs = GetSecs() - obj.tStart_secs;
            if tElapsed_secs < obj.ONRAMP_SECS
            	w = .5*(1 - cos(2*pi*tElapsed_secs/(obj.ONRAMP_SECS*2)));
            elseif tElapsed_secs < (obj.STIM_DURATION_SECS - obj.OFFRAMP_SECS)
                w = 1;
            elseif tElapsed_secs < obj.STIM_DURATION_SECS
                tRemaining_secs = obj.STIM_DURATION_SECS - tElapsed_secs;
                w = 1 - .5*(1 - cos(2*pi*tRemaining_secs/(obj.ONRAMP_SECS*2)));
            else
                w = 0;
            end
        end
    end


    %% ====================================================================
    %  -----STATIC METHODS-----
    %$ ====================================================================    
    
    methods (Static, Access = public)
    end
    
end