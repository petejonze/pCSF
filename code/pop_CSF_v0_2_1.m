function [] = pop_CSF_v0_2_1(pid, DOB, eyeStr)
% Estimating CSFs by popping gabors on a touchscreen
%
% Requires:         Mario Kleiner's Psychtoolbox (v3+)
%                   my OpenFace toolbox (for head tracking)
%                   my QUEST+ toolbox
%
% Matlab:           v2012+
%
% See also:         none
%
% Example:          pop_CSF_v0_2_1()
%
% Author(s):    	Pete R Jones <petejonze@gmail.com>
%
% Version History:  0.0.1	PJ  18/04/2018	Initial build
%                   0.0.2	PJ  19/04/2018	added Gabor class
%                   0.0.3	PJ  19/04/2018	misc
%                   0.0.4	PJ  20/04/2018	added feedback media
%                   0.0.5	PJ  20/04/2018	added QUEST+, OpenFace, and DataManager, Score (display)
%                   0.0.6	PJ  17/05/2019	for piloting
%                   0.1.0	PJ  17/09/2019	Updated data entry, for Doaa's supression study
%                   0.2.0	PJ  06/11/2019	now saving trial by trial data (with explict timestamp linkage between OpenFace and trials). Also now you lose points for incorrect answers. This update also involved tweaks to OpenFace.m and Stimulus.m
%                   0.2.1   PJ  26/11/2019  stopped score ever going below 0. Changed bug in saved file name
%
%
% Copyright 2020 : P R Jones <petejonze@gmail.com>
% *********************************************************************
% 

% @TODO randomly add easy trials (for motivation)
% @TODO log every screen press (including timestamps)
% @TODO start logging IPD?

    %%%%%%%%
    %% -1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% pre-init  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % pre-init
        clc
        close all
        clearJavaMem();
        % clear all;
        
        % suppress unnecessary warnings in Editor
        %#ok<*UNRCH>


    %%%%%%%
    %% 0 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Params  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        DO_HEADTRACKING = true;
        USE_HEADTRACKING_TO_PERFORM_REALTIME_SCALING = false;
        DO_PLOT = true;
        IS_DEBUG_MODE = false;
        
        graphicParams.screenNum = min(1, max(Screen('Screens'))); % 0 or 1
        graphicParams.monitorWidth_cm = 26.0;
        graphicParams.testScreenWidth = 2736;
        graphicParams.viewDist_cm = 50;
        background_CL = [0.5 0.5 0.5];
        %background_CL = [0.1 0.1 0.1]; <-- not sure this actualy works
        nFramesPerSecond = 60;
        ifi_hz = 1/nFramesPerSecond;
        gabor_supportSize_px = 350;
        
        newStimPerSecond = 1;
        maxNStim = 5;
        
        % internal
        pNewStim_baseRate = newStimPerSecond/nFramesPerSecond;
  
        [screenWidth_px, screenHeight_px]=Screen('WindowSize', graphicParams.screenNum);
        [xPix_px, yPix_px] = meshgrid(gabor_supportSize_px:10:(screenWidth_px-gabor_supportSize_px), gabor_supportSize_px:10:(screenHeight_px-gabor_supportSize_px)); % sparse sample

        % validate
        if USE_HEADTRACKING_TO_PERFORM_REALTIME_SCALING && ~DO_HEADTRACKING
            error('inconsistent inputs');
        end
            
    %%%%%%%
    %% 1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Init  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % make sure in this file's directory
        fn = mfilename('fullpath');
        pathstr = fileparts(fn);
        cd(pathstr);

        % init head tracking
        if DO_HEADTRACKING
            deviceID = 0; % 1 for my SPro, 0 for Doaa's SPro
            OpenFace.init(deviceID, false); % true to enable automatic plotting
            headPose_updateTimeLog_secs = [];
        end

        % assert PTB installed
        AssertOpenGL();
        
        % parse inputs; get participant ID (integer)
        if nargin < 1 || isempty(pid)
            while 1
                pid = input('Participant ID (integer): ', 's');
                if ~isempty(pid)
                    pid = str2double(pid);
                    if ~isempty(pid) && ceil(pid) == floor(pid)
                        break
                    else
                        warning('input must be numeric, and must be an integer');
                    end
                end
            end
        end
        
        % parse inputs; get participant DOB (string)
        if nargin < 2 || isempty(DOB)
            while 1
                DOB = input('DOB: ', 's');
                if ~isempty(DOB)
                    break
                end
            end
        end
        
        % parse inputs; get eye
        if nargin < 3 || isempty(eyeStr)
            while 1
                eyeStr = upper(input('[L]eft or [R]ight Eye? ', 's'));
                if ismember(eyeStr, {'L','R'})
                    break
                end
            end
        end

        % init timers
        start_dateTime = datestr(now(),30);
        
        % init data manager
        DataManager.init();
                
        % set windows screen brightness to max (though NB: may still be below
        % max if not fully charged and/or plugged in, depending on various
        % arcane windows settings
        system('powershell -inputformat none (Get-WmiObject -Namespace root/WMI -Class WmiMonitorBrightnessMethods).WmiSetBrightness(1,100)');
        
   
    %%%%%%%
    %% 2 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialise QuestPlus (compute stimulus range) %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % screen
        monitorWidth_cm = graphicParams.monitorWidth_cm;
        screenWidth_px = graphicParams.testScreenWidth;
        viewDist_cm = graphicParams.viewDist_cm;
        pixel_per_cm = screenWidth_px/monitorWidth_cm
        screenWidth_dg = 2*rad2deg(atan(monitorWidth_cm/(2*viewDist_cm)))
        pixel_per_dg = screenWidth_px/screenWidth_dg
        maxFreq_cpd = pixel_per_dg/2 % cycles per degree

        % gabor
        gabor_rotation_deg = 90;
        gabor_phase = 0;
        gabor_contrast = 1; % max 1
        gabor_supportWidth_px = gabor_supportSize_px;
        gabor_supportHeight_px = gabor_supportSize_px;
        gabor_sc_px = round(gabor_supportSize_px/6.667); % 8 seems about right
        gabor_aspectratio = 1;

        % stimuli
        ppc = round(logspace(log10(8),log10(gabor_sc_px), 10)) % minimum of 4 pixels per cycle (i.e., 2 per band), maximum of N pixels per cycle (i.e., minimum of ~4 bands per gabor)
        ppc = unique(ppc + mod(ppc, 2)) % round up to nearest even number
        cpp = 1./ppc
        cpd = cpp .* pixel_per_dg
        cpd = fliplr(cpd)
        
        % QuestPlus
        QP = QuestPlusFactory.createQuestPlus(cpd);
        QP_startGuess_mean = QP.getParamEsts('mean')
       
        
    %%% TRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	try
       
        
    %%%%%%%
    %% 4 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Launch OpenGL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Disable screen checks
        Screen('Preference', 'SkipSyncTests', 1);

        % set font size
        Screen('Preference', 'DefaultFontSize', 16);
        
        % Define requirements for onscreen window - Setup imaging pipeline:
        PsychImaging('PrepareConfiguration');

        % Want a full 32 bit floating point precision framebuffer:
        % This will provide an effective resolution of 23 bits in the displayable
        % luminance range -- Plenty of headroom for the up-to 16 bits output devices.
        % Hardware older than NVidia Geforce 8000 or ATI Radeon HD 2000 won't be
        % able to do alpha-blending at this mode though. Not an issue here, as we
        % don't need alpha-blending...
        PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');

        % Ensure that we have sufficient contrast resolution to exceed the
        % human threshold (~9 bits)
        % PsychImaging('AddTask', 'General', 'EnableNative10BitFramebuffer');
        PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput');

        % Open the onscreen window, get its handle and bounding rectangle. We open
        % with a black background ( == 0 ) on display screen 'screenid':
        winhandle = PsychImaging('OpenWindow', graphicParams.screenNum, background_CL);
        %winhandle = PsychImaging('OpenWindow', graphicParams.screenNum, background_CL, [0 0 2735 1823]); % for screen recording
        %winhandle = PsychImaging('OpenWindow', graphicParams.screenNum, background_CL, [0 0 1000 1000]); % for debugging


        % Allow transparency
        Screen('BlendFunction', winhandle,  GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        % initialise graphics manager
        media.GraphicsManager(winhandle);

      	% hide mouse cursor
        HideCursor();
        
        % calib
        originalGammaTable = Screen('ReadNormalizedGammaTable', winhandle);
        calib = load('../calib/rawMeasures-bitsteal_CVL-SP4-CRS_28May2019_20181128T061607.mat');
        % tmp hack: for now we shall just assume uniformity, and average
        % across locations
        n = length(calib.in_CL);   
        in_CL = calib.in_CL;
        out_cdm2 = mean(reshape(calib.out_cdm2,[],n));
        %figure(); plot(in_CL, out_cdm2, '-o');
        % fit
        fittedmodel = fit(out_cdm2(:)./max(out_cdm2(:)), in_CL(:), 'splineinterp'); % 'linear');
        estVals = fittedmodel(linspace(0,1,256));
        estVals(estVals<0) = 0;
        Screen('LoadNormalizedGammaTable', winhandle, estVals*[1 1 1],1);

        
    %%%%%%%
    %% 5 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Init other hardware  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % keyboard
        KbName('UnifyKeyNames');
        spaceKey = KbName('SPACE');
        %enterKey = KbName('return');
        escapeKey = KbName('escape');
        %finKey = KbName('f');
        
        % touchscreen
        touchScreen_deviceNumber = 1;
        try
            PsychHID('KbQueueCreate', touchScreen_deviceNumber)
            PsychHID('KbQueueStart', touchScreen_deviceNumber)
            PsychHID('KbQueueFlush', touchScreen_deviceNumber)
        catch ME
            LoadPsychHID()
            rethrow(ME)
        end
            
        
        % double-tap suppression
        lastPressTime_secs = -inf;
        lastPressXY_px = [0 0];
        
        % pre-cache audioplayer
        x = getPureTone(1000, 44100, 0.1, 0.01) * .01;
        pahandle = audioplayer(x, 44100);
    	pahandle.play();
        WaitSecs(0.1);
        
        media.feedback.Feedback([0 0], true);
        

    %%%%%%%
    %% 6 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Init media  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        % Build a procedural gabor texture for a gabor with a support of tw x
        % th pixels and the 'nonsymetric' flag set to 0 to force only symmetric
        % aspect ratios. Since the texture is procedurally generated, the
        % precise parameters can be set dynamically at runtime (i..e during
        % the call to DrawTexture)
        nonSymmetric                = 0;
        backgroundColorOffset       = [background_CL 1]; % e.g., [.5 .5 .5 1]; % 1 for full opacity
        gabor_modulateColor         = [1 1 1 0];
        disableNorm                 = 1;
        contrastPreMultiplicator    = mean(background_CL); % e.g., 0.5; % CHANGE
        validModulationRange        = [-1 1];
        [gaborid, gaborrect] = CreateProceduralGabor(winhandle, gabor_supportWidth_px, gabor_supportHeight_px, nonSymmetric, backgroundColorOffset, disableNorm, contrastPreMultiplicator, validModulationRange);

        % init LJP: movement handler
        x = []; % randsample(1:screenWidth_px, 1);
        y = []; % randsample(1:screenHeight_px, 1);
        v = []; % 500 % (rand()-0.5)*4;
        u = []; % 500 % (rand()-0.5)*4;
        m = []; % 1;
        borders = [0 screenWidth_px 0 screenHeight_px];
      	d = gabor_supportWidth_px*2;                    % particle diameter
        dt = ifi_hz;                         	% time step in seconds
        epsilon = [];                           %material parameter (?)
        maxVelocity = 500;
        myLJP = LennardJonesPotential(x,y,v,u,m,borders,d,dt,epsilon,maxVelocity);
        Stimuli = {};
                
        % data to store
        stimLog                 = [];
        stimLog.phase           = [];
        stimLog.freq_cpd       	= [];
        stimLog.freq_cpp       	= [];
        stimLog.spatialconstant = [];
        stimLog.contrast        = [];
        stimLog.aspectRatio   	= [];
        stimLog.tOnset_secs 	= [];
        stimLog.tOffset_secs 	= [];
        stimLog.response        = [];
        stimLog.estHPoseAtOnset	= [];
        stimLog.estHPoseAtOffset= [];
        stimLog.xyt_px          = {};
        % counter
        nIncorrectPresses = 0;
        nCorrectPresses = 0;
        
        % HUD (??)
        media.score.Score();
        
        
    %%%%%%%
    %% 7 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Run Experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % wait for input to start
        Screen('TextSize', winhandle, 96);
        Screen('DrawText', winhandle, 'Press Screen To Start', 100, 800, [0 0 0]);
        Screen('Flip', winhandle);
        while 1
            % Check for key press (quit)
            [keyIsDown, ~, keyCode ] = KbCheck;
            keyCode = find(keyCode, 1);            
            if keyIsDown
                switch keyCode
                    case escapeKey
                        error('Aborted By User');
                    otherwise
                        % do nothing
                end
            end
            
            % get press
            GetMouse();
            if PsychHID('KbQueueCheck', touchScreen_deviceNumber)
                break
            end
            % short pause
            WaitSecs(ifi_hz);
        end

        % init
        tStart = tic();
        if DO_HEADTRACKING
            OpenFace.reset(); % clear OpenFace buffer
            headPose_startTime_secs = GetSecs();
        end

        % run
        while 1
            % update OpenFace buffer
            if DO_HEADTRACKING
                newDat = OpenFace.update();
                if ~isempty(newDat)
                    n = size(newDat,1);
                    headPose_updateTimeLog_secs = [headPose_updateTimeLog_secs; ones(n,1)*GetSecs()]; %#ok
                end
            end
            
            % Check for key press (quit)
            [keyIsDown, ~, keyCode ] = KbCheck;
            keyCode = find(keyCode, 1);
            if keyIsDown
                switch keyCode
                    case escapeKey
                        break
                    case spaceKey
                        imageArray = Screen('GetImage', winhandle);
                        imwrite(imageArray, fullfile('screenshots', 'printscreen.jpg'));
                    otherwise
                        % do nothing
                end
            end
            
            % get any screen press and process it accordingly
         	[x,y] = GetMouse(); % query mouse (& get tap location if any)
            if PsychHID('KbQueueCheck', touchScreen_deviceNumber)

                % check if a double tap
                t_secs = GetSecs()-lastPressTime_secs;
                d_px = sqrt((lastPressXY_px(1) - x).^2 + (lastPressXY_px(2) - y).^2);
                if t_secs<0.75 && d_px<200
                    % probably double tap, ignore
                else
                    % ok, proceed
                    lastPressTime_secs = GetSecs();
                    lastPressXY_px = [x y];
                    
                    % check for Hit
                    d = sqrt((myLJP.x - x).^2 + (myLJP.y - y).^2);
                    idx = find(d < (gabor_sc_px*3*1.4142)); % a catch radius of 3*sqrt(2) seems about right
                    if any(idx)
                        idx1 = find(idx);
                        for i = 1:length(idx1)
                            nCorrectPresses = nCorrectPresses + 1;
                            % feedback
                            media.feedback.Feedback([x(idx1(i)) y(idx1(i))], true);
                            % update score
                            DataManager.addPoints(1);
                            % update psychophysical algorithm
                            % < !!! UPDATE QUEST+ WITH A HIT HERE !!! >
                            stim = [Stimuli{idx1(i)}.freq_cpp*pixel_per_dg; Stimuli{idx1(i)}.contrast];
                            anscorrect = true;
                            QP.update(stim, anscorrect);
                            
                            % log Hit
                            stimLog.phase(end+1)            = Stimuli{idx1(i)}.getPhase();
                            stimLog.freq_cpd(end+1)       	= Stimuli{idx1(i)}.getFreq_cpd();
                            stimLog.freq_cpp(end+1)       	= Stimuli{idx1(i)}.getFreq_cpp();
                            stimLog.spatialconstant(end+1)  = Stimuli{idx1(i)}.getSpatialConstant();
                            stimLog.contrast(end+1)         = Stimuli{idx1(i)}.getContrast();
                            stimLog.aspectRatio(end+1)   	= Stimuli{idx1(i)}.getAspectRatio();
                            stimLog.tOnset_secs(end+1)     	= Stimuli{idx1(i)}.tStart_secs;
                            stimLog.tOffset_secs(end+1)     = GetSecs();
                            stimLog.response(end+1)         = 1;
                            if DO_HEADTRACKING
                                estHPoseAtOnset = Stimuli{idx1(i)}.estHPoseAtOnset;
                                if isempty(estHPoseAtOnset)
                                    estHPoseAtOnset = nan(1,6);
                                end
                                stimLog.estHPoseAtOnset(end+1,:) 	= estHPoseAtOnset;
                                estHPoseAtOffset = OpenFace.getBuffer(true);
                                if isempty(estHPoseAtOffset)
                                    estHPoseAtOffset = nan(1,6);
                                end
                                stimLog.estHPoseAtOffset(end+1,:)	= estHPoseAtOffset;
                            else
                                stimLog.estHPoseAtOnset(end+1,:) 	= nan(1,6);
                                stimLog.estHPoseAtOffset(end+1,:) 	= nan(1,6);
                            end
                            stimLog.xyt_px{end+1}            = Stimuli{idx1(i)}.getLocationLog();
                        end
                        
                        % remove (all)
                        myLJP.removeMolecule(idx);
                        Stimuli(idx) = []; %#ok
                    else
                        media.feedback.Feedback([x y], false);
                      	DataManager.addPoints(-1);
                        nIncorrectPresses = nIncorrectPresses + 1;
                    end
                end
            end

            % draw background
            Screen('FillRect', winhandle , background_CL); % gray background
                 
            % draw gabor(s) [assuming there are any]
            if (myLJP.N > 0)
                % position params
                dstRect = CenterRectOnPointd(gaborrect, myLJP.x, myLJP.y)';
                angles_deg = ones(1,myLJP.N)*gabor_rotation_deg;
                
                % gabor params: [gabor_phase; gabor_freq_cpp; gabor_sc_px; gabor_contrast; gabor_aspectratio; 0; 0; 0]
                pars = zeros(8, myLJP.N); 
                for j = 1:myLJP.N
                    pars(1, j) = Stimuli{j}.getPhase();             % get phase
                    pars(2, j) = Stimuli{j}.getFreq_cpp();        	% get frequency
                    pars(3, j) = Stimuli{j}.getSpatialConstant(); 	% get spatial constant
                    pars(4, j) = Stimuli{j}.getContrast();          % get contrast
                    pars(5, j) = Stimuli{j}.getAspectRatio();       % get aspect ratio
                    
                    % log stimulus locations
                    Stimuli{j}.logLocation(myLJP.x(j), myLJP.y(j), GetSecs());
                end
                
                % draw
                Screen('DrawTextures', winhandle, gaborid, [], dstRect, angles_deg, [], [], gabor_modulateColor, [], kPsychDontDoRotation, pars);
            end
            
            % draw graphics (feedback, HUD, etc.)
            media.GraphicsManager.getInstance().drawAll();
            
            % Mark drawing ops as finished, so the GPU can do its drawing job while
            % we can compute updated parameters for next animation frame. This
            % command is not strictly needed, but it may give a slight additional
            % speedup, because the CPU can compute new stimulus parameters in
            % Matlab, while the GPU is drawing the stimuli for this frame.
            % Sometimes it helps more, sometimes less, or not at all, depending on
            % your system and code, but it only seldomly hurts.
            % performance...
            Screen('DrawingFinished', winhandle);

            % update graphic params
            [x,y] = myLJP.calcNextTimestep();
            
            % check for expired stims
            myLJP.N
            if (myLJP.N > 0)
                idx = zeros(myLJP.N, 1);
                for i = 1:myLJP.N
                    idx(i) = Stimuli{i}.checkForTimeout();
                end
                idx = idx == 1;
                
                % remove any time outs
                if any(idx)
                    % < !!! UPDATE QUEST+ WITH A MISS HERE !!! >
                    idx1 = find(idx);
                    for i = 1:length(idx1)
                        stim = [Stimuli{idx1(i)}.freq_cpp*pixel_per_dg; Stimuli{idx1(i)}.contrast];
                        anscorrect = false;
                        QP.update(stim, anscorrect);
                        
                        % log Miss
                        stimLog.phase(end+1)            = Stimuli{idx1(i)}.getPhase();
                        stimLog.freq_cpd(end+1)        	= Stimuli{idx1(i)}.getFreq_cpd();
                        stimLog.freq_cpp(end+1)        	= Stimuli{idx1(i)}.getFreq_cpp();
                        stimLog.spatialconstant(end+1)  = Stimuli{idx1(i)}.getSpatialConstant();
                        stimLog.contrast(end+1)         = Stimuli{idx1(i)}.getContrast();
                        stimLog.aspectRatio(end+1)   	= Stimuli{idx1(i)}.getAspectRatio();
                        stimLog.tOnset_secs(end+1)    	= Stimuli{idx1(i)}.tStart_secs;
                        stimLog.tOffset_secs(end+1) 	= GetSecs();
                        stimLog.response(end+1)         = 0;
                        if DO_HEADTRACKING
                            estHPoseAtOnset = Stimuli{idx1(i)}.estHPoseAtOnset;
                            if isempty(estHPoseAtOnset)
                                estHPoseAtOnset = nan(1,6);
                            end
                            stimLog.estHPoseAtOnset(end+1,:) 	= estHPoseAtOnset;
                            estHPoseAtOffset = OpenFace.getBuffer(true);
                            if isempty(estHPoseAtOffset)
                                estHPoseAtOffset = nan(1,6);
                            end
                            stimLog.estHPoseAtOffset(end+1,:)	= estHPoseAtOffset;
                        else
                            stimLog.estHPoseAtOnset(end+1,:) 	= nan(1,6);
                            stimLog.estHPoseAtOffset(end+1,:) 	= nan(1,6);
                        end
                        stimLog.xyt_px{end+1}            = Stimuli{idx1(i)}.getLocationLog();
                    end
                    
                    % remove (all)
                    myLJP.removeMolecule(idx);
                    Stimuli(idx) = []; %#ok
                    
                    % if in debug mode play a beep to alert user that a
                    % stim was missed
                    if IS_DEBUG_MODE
                        beep();
                    end
                end
            end
            
            % check for new stim
            if (myLJP.N < maxNStim) && (rand() < (pNewStim_baseRate/myLJP.N))
                fprintf('Placing new stimulus\n');

                % find possible locations
                if (myLJP.N > 0)
                    d_all = sqrt(bsxfun(@minus, xPix_px, permute(x,[3 2 1])).^2 + bsxfun(@minus, yPix_px, permute(y,[3 2 1])).^2); % ALT: reshape(x,[1 1 length(x)
                    d = min(d_all, [], 3);
                    idx_possible = d > gabor_supportSize_px; % non-overlapping locations not allowed
                else % all possible
                    idx_possible = xPix_px > 0;
                end
                
                % place stimulus
                if sum(idx_possible) == 0
                    warning('no valid locations! Stimuli not placed');
                else
                   	% pick random location
                    idx = randsample(find(idx_possible), 1);
                    x = xPix_px(idx);
                    y = yPix_px(idx);
                    % random starting velocity
                    v = 500 * (1 - rand()*2); % -x to x
                    u = 500 * (1 - rand()*2); % -x to x
                    % add (1 of 2)
                    myLJP.addMolecule(x,y, v,u, 1);

                    % add (2 of 2)
                    % < !!! QUERY QUEST+ FOR FREQ AND CONSTRAST HERE !!! >
                    %gabor_freq_cpp = randsample(cpp, 1); % TMP HACK
                    %gabor_contrast = rand(); % TMP HACK
                    stim = QP.getTargetStim();
                    gabor_freq_cpd = stim(1);
                    gabor_contrast = stim(2);
                    
                    % convert cpd to cpp
                    if DO_HEADTRACKING
                        estHPoseAtOnset = OpenFace.getBuffer(true);
                        % validate (1/2)
                        if USE_HEADTRACKING_TO_PERFORM_REALTIME_SCALING && ~isempty(estHPoseAtOnset)
                            estDist_cm = estHPoseAtOnset(3)/100;
                            % validate (2/2)
                            tAge_secs = GetSecs()-headPose_updateTimeLog_secs(end);
                            if ~isnan(estDist_cm) && estDist_cm>40 && estDist_cm<80 && tAge_secs<5
                                % update viewing distance info
                                screenWidth_dg = 2*rad2deg(atan(monitorWidth_cm/(2*viewDist_cm)));
                                pixel_per_dg = screenWidth_px/screenWidth_dg;
                            end
                        end
                    else
                        estHPoseAtOnset = [];
                    end
                    % <<convert>>
                    gabor_freq_cpp = gabor_freq_cpd / pixel_per_dg;
   
                    % create stimulus
                    Stimuli{end+1} = Stimulus(gabor_phase, gabor_freq_cpd, gabor_freq_cpp, gabor_sc_px, gabor_contrast, gabor_aspectratio, estHPoseAtOnset); %#ok
                end
            end
        
            % flip
            Screen('Flip', winhandle);

            % short pause
            WaitSecs(ifi_hz);
        end

        if DO_HEADTRACKING
            headPoseDat = OpenFace.getBuffer();
            headPose_stopTime_secs = GetSecs();
        end
        
        % Close up
        Screen('LoadNormalizedGammaTable', winhandle, originalGammaTable);
        sca();
        PsychHID('KbQueueStop', touchScreen_deviceNumber)
        PsychHID('KbQueueRelease', touchScreen_deviceNumber)
        SingletonManager.clearAll();
        ShowCursor();
        if DO_HEADTRACKING
            % kill
            OpenFace.finishUp();
        end
        
    catch ME
        try
            Screen('LoadNormalizedGammaTable', winhandle, originalGammaTable);
        catch
        end;
        sca();
        try
            PsychHID('KbQueueStop', touchScreen_deviceNumber)
            PsychHID('KbQueueRelease', touchScreen_deviceNumber)
        catch
        end;
        try
            SingletonManager.clearAll();
        catch
        end;
        ShowCursor();
        if DO_HEADTRACKING
            % kill
            OpenFace.finishUp();
        end
                
        rethrow(ME);
    end
    
    
    %%%%%%%
    %% 8 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Compute results  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        results = [];
        
        if DO_HEADTRACKING
            results.headPoseDat                     = headPoseDat;
            results.headPose_startTime_secs         = headPose_startTime_secs;
            results.headPose_stopTime_secs          = headPose_stopTime_secs;
            results.headPose_updateTimeLog_secs     = headPose_updateTimeLog_secs;
        end
        
        % compute test duration
        totalTestDuration_secs = toc(tStart);
                        
        % get final parameter estimates
        QP_endGuess_mean = QP.getParamEsts('mean');

        % display
        QP.disp();
        fprintf('     Gmax estimate: %1.2f	[start: %1.2f]\n', QP_endGuess_mean(1), QP_startGuess_mean(1));
        fprintf('     Fmax estimate: %1.2f	[start: %1.2f]\n', QP_endGuess_mean(2), QP_startGuess_mean(2));
        fprintf('        B estimate: %1.2f	[start: %1.2f]\n', QP_endGuess_mean(3), QP_startGuess_mean(3));
        fprintf('  pf_beta estimate: %1.2f	[start: %1.2f]\n', QP_endGuess_mean(4), QP_startGuess_mean(4));
        fprintf(' pf_gamma estimate: %1.2f	[start: %1.2f]\n', QP_endGuess_mean(5), QP_startGuess_mean(5));
        fprintf('pf_lambda estimate: %1.2f	[start: %1.2f]\n', QP_endGuess_mean(6), QP_startGuess_mean(6));


        % extract data ------------------------------------------------
        % compute S_start
        Gmax    = QP_startGuess_mean(1);	% peak gain (sensitivity): 2 -- 2000
        Fmax    = QP_startGuess_mean(2); 	% peak spatial frequency: 0.2 to 20 cpd
        B       = QP_startGuess_mean(3);	% bandwidth (full width at half maximum): 1 to 9 octaves
        f = logspace(log10(cpd(1)), log10(cpd(end)), 1000);
        S_start = log10(Gmax) - log10(2) * ( (log10(f) - log10(Fmax)) / (log10(2*B)/2) ).^2;
        S_start(f<Fmax) = log10(Gmax);
        % store
        results.start.Gmax = Gmax;
        results.start.Fmax = Fmax;
        results.start.B = B;
        results.start.f = f;
        results.start.S = S_start;
        
        % history
        results.QP_history_stim	= QP.history_stim;
        results.QP_history_resp	= QP.history_resp;
        
        % further history (manual logs)
        results.stimLog             = stimLog;
        results.nCorrectPresses     = nCorrectPresses;
        results.nIncorrectPresses   = nIncorrectPresses;
        
        % compute S_end
        Gmax    = QP_endGuess_mean(1);     % peak gain (sensitivity): 2 -- 2000
        Fmax    = QP_endGuess_mean(2); 	% peak spatial frequency: 0.2 to 20 cpd
        B       = QP_endGuess_mean(3);     % bandwidth (full width at half maximum): 1 to 9 octaves
        f = logspace(log10(cpd(1)), log10(cpd(end)), 1000);
        S_end = log10(Gmax) - log10(2) * ( (log10(f) - log10(Fmax)) / (log10(2*B)/2) ).^2;
      	S_end(f<Fmax) = log10(Gmax);
        % store
        results.end.Gmax = Gmax;
        results.end.Fmax = Fmax;
        results.end.B = B;
        results.end.f = f;
        results.end.S = S_end;
        
        % unlog units
        S_start = exp10(S_start);
        S_end	= exp10(S_end);

    %%%%%%%
    %% 5 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % construct file name and add path
        fn = sprintf('pCSF_Sup%i_%s_%s.mat', pid, eyeStr, start_dateTime);
        fullFn = fullfile('..','data', fn);

        % get any final comments
        notes = input('Experimenter comments (optional): ', 's');
        
        % organise data for saving
        tmp.test = results;
        tmp.totalTestDuration_secs = totalTestDuration_secs;
        tmp.pid = pid;
        tmp.DOB = DOB;
        tmp.testDateTime = start_dateTime;
        tmp.notes = notes;
        tmp.eyeStr = eyeStr; %#ok

        % save
        fprintf('saving data..\n');
        save(fullFn, '-struct', 'tmp');
        fprintf('..done\n\n');
        fprintf('   x = load(''%s'')\n', fullFn);
    


    %%%%%%%
    %% 6 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot results  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if DO_PLOT
            % plot
            figure()
            hold on
            plot(f, S_start, 'k');
            plot(f, S_end, 'b:');

            % annotate and format
            legend('Prior','Empirical', 'Location','South');
            set(gca, 'XScale','log', 'YScale','log');
            set(gca, 'XTick',[.5 1 2 5 10 20], 'YTick',[2 10 50 300 2000]);
            xlabel('Spatial Frequency (cpd)'); ylabel('Contrast Sensitivity (1/C)')
            xlim([.25 40]); ylim([1 2000]);
        end
        
    %%%%%%%
    %% 7 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Finishup  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('All done. Thanks for playing!\n');
    
end