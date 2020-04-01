classdef (Sealed) QuestPlusFactory < handle
    % Sstatic class for initialising QuestPlus
    %
    %   dfdfdf
    %
    % See Also:
    %   mess_3param_v2.m
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
        USE_PRIOR = true; % Uniform if false
    end
    
    %% ====================================================================
    %  -----STATIC METHODS (public)-----
    %$ ====================================================================
    
    methods (Static, Access = public)

        function pC = qCSF_getPC3(sf,c, Gmax,Fmax,B, pf_beta,pf_gamma,pf_lambda)
            % CSF
            S = log10(Gmax) - log10(2) * ( (log10(sf) - log10(Fmax)) / (log10(2*B)/2) ).^2;
            if (sf<Fmax)
                S = log10(Gmax);
            end

            % convert 'sensitivity' to 'contrast'
            S = exp10(S); % linearize (convert dB to linear units)
            alpha = 1/S; % convert sensitivity to %contrast threshold (e.g., 2=>50; 10=>10; 200=>0.5)
            
            % PF
            pC = pf_gamma+(1 - pf_gamma - pf_lambda).*(1-exp(-10.^(pf_beta.*(log10(c)-log10(alpha)))));
        end
        
        function QP = createQuestPlus(f_cpd)
            
            % set model
            F = @QuestPlusFactory.qCSF_getPC3;
            
            % define parameter domain
            Gmax        = logspace(log10(3), log10(300), 15);	% peak gain (sensitivity): 2 -- 2000
            Fmax        = logspace(log10(1), log10(10),  10); 	% peak spatial frequency: 0.2 to 20 cpd
            B           = linspace(1, 9, 9);                    % bandwidth (full width at half maximum): 1 to 9 octaves
            pf_beta     = 3;                                    % [fixed] psychometric function: slope
            pf_gamma    = 0.05;                                 % [fixed] psychometric function: guess rate
            pf_lambda   = [0.05 0.1 0.2];                    	% [rigid] psychometric function: lapse rate
            paramDomain = {Gmax, Fmax, B, pf_beta, pf_gamma, pf_lambda};

            % define testing domain
%             stimDomain =    {logspace(log10(2), log10(30), 15)      ... % spatial frequency
%                             ,logspace(log10(0.01), log10(100), 15)	... % contrast
%                             };
            stimDomain =    {f_cpd                                  ... % spatial frequency
                            ,logspace(log10(0.01), log10(100), 15)	... % contrast
                            };
                        
            % define response domain
            respDomain = [0 1];
            
            % define priors
          	priors = cell(size(paramDomain));
            if QuestPlusFactory.USE_PRIOR
                priors{1} = normpdf(paramDomain{1}, 3, 100); % low-ball
                priors{2} = normpdf(paramDomain{2}, 7, 2*2);
                priors{3} = normpdf(paramDomain{3}, 3, 1*2);
                % normalise so sum to 1
                priors{1} = priors{1}./sum(priors{1});
                priors{2} = priors{2}./sum(priors{2});
                priors{3} = priors{3}./sum(priors{3});
                priors{4} = 1;
                priors{5} = 1;
                priors{6} = ones(1,length(paramDomain{6}))/length(paramDomain{6}); % uniform
            else
                priors{1} = ones(1,length(paramDomain{1}))/length(paramDomain{1}); % uniform
                priors{2} = ones(1,length(paramDomain{2}))/length(paramDomain{2}); % uniform
                priors{3} = ones(1,length(paramDomain{3}))/length(paramDomain{3}); % uniform
                priors{4} = ones(1,length(paramDomain{4}))/length(paramDomain{4}); % uniform
                priors{5} = ones(1,length(paramDomain{5}))/length(paramDomain{5}); % uniform
                priors{6} = ones(1,length(paramDomain{6}))/length(paramDomain{6}); % uniform
            end
            
            % define other parameters (blank for default)
            stopRule        = [];
            stopCriterion   = [];
            minNTrials      = 50;
            maxNTrials      = 100;
            
            % create QUEST+ object
            QP = QuestPlus(F, stimDomain, paramDomain, respDomain, stopRule, stopCriterion, minNTrials, maxNTrials);
            
            % initialise priors/likelihoods
            fn = './myLikelihoods.mat';
            if exist(fn, 'file')
                file = dir(fn);
                fprintf('Using precomputed likelihoods [from %s]\n', file.date)
                QP.initialise(priors, fn)
            else
                QP.initialise(priors);
                QP.saveLikelihoods(fn);
            end
        end
        
    end
    
end