classdef FeedbackMedia < handle
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
        IMG_HIT_FN = 'gold_coin.png';
        IMG_MISS_FN = '200px-Sad_face.svg.png';
        SND_HIT_FN = 'pop.wav';
        SND_MISS_FN = 'bad.wav';
    end

    properties (GetAccess = public, SetAccess = private)
        img_hit
        img_miss
        snd_hit
        snd_miss
    end
    
    
    %% ====================================================================
    %  -----CONSTRUCTOR/DESTRUCTOR METHODS-----
    %$ ====================================================================    

    methods (Access = public)
        
        function obj = FeedbackMedia()
            localDir = fileparts(mfilename('fullpath')); % get directory in which this class file resides
            
            % initialise graphic(s)
            gm = media.GraphicsManager.getInstance();
            obj.img_hit = PtrVisual(fullfile(localDir, obj.IMG_HIT_FN), gm.winhandle, 1, [0 0]);
            obj.img_miss = PtrVisual(fullfile(localDir, obj.IMG_MISS_FN), gm.winhandle, 1, [0 0]);

            % initialise sound
            obj.snd_hit = PtrSound(fullfile(localDir, obj.SND_HIT_FN));
            obj.snd_miss = PtrSound(fullfile(localDir, obj.SND_MISS_FN));
        end
        
    end        

end