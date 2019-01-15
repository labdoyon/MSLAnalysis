function [perf_duration,...
            all_keys,...
            correct_keys,...
            incorrect_keys,...
            sequence_keys,...
            error_keys...
            ]...
            = get_keys_stat(keys, onsets, sequence, n_start_trial)

% INPUT
% keys                      a vector of keys
% onsets                    a vector of key onsets, i.e., the time that the key was presssed
% sequence      [integer]   a vector of numbers representing the sequence
% n_start_trial [integer]   the number of the first keys of the sequence to search for a trial; the default value is 2 keys
% 
% OUTPUT
% perf_duration  	[integer]    time between the first and the last presses key
% 
% all_keys    
%   ....n       [integer]
%   ....mean  	[double]    mean transition duration
%   ....sd   	[double]
%
% correct_keys              keys of correctly performed and completed or incomleted sequences are correct; keys within "head" and "tail" are also correct
%   ....n       [integer]
%   ....mean  	[double]    mean transition duration between correct keys; transitions between correct and incorrect keys are not considered
%   ....sd   	[double]
%
% incorrect_keys          	icorrect keys
%   ....n       [integer]
%   ....mean  	[double]    mean transition duration between incorrect keys; transitions between correct and incorrect keys are not considered
%   ....sd   	[double]
%
% sequence_keys          	keys within correctly performed and completed sequences
%   ....n       [integer]
%   ....mean  	[double]    mean transition duration between sequence keys and transitions between sequences or attempt to perform a sequence;
%                           such attempt is identified by the <n_start_trial> first keys of the sequence
%   ....sd   	[double]
%
% error_keys                keys within errors
%   ....n       [integer]
%   ....mean  	[double]    mean transition duration between error keys or transitons between errors; keys between sequence and error are not considered
%   ....sd   	[double]

% Ella Gabitov, 14 January, 2019

if numel(keys) ~= numel(onsets)
    error('The number of keys and the number of onsets do not match.');
end

if nargin < 4, n_start_trial = 2; end;
if isempty(n_start_trial) || isnan(n_start_trial) || n_start_trial == 0, n_start_trial = 2; end;

perf_duration = onsets(end) - onsets(1);

%% ALL KEYS

transition_durations = [];
for i_trans = 1 : numel(keys)-1
    transition_durations = [transition_durations onsets(i_trans+1) - onsets(i_trans)];
end
all_keys.n = numel(keys);                           % # of all keys
all_keys.mean = nanmean(transition_durations); 	% mean transition duration
all_keys.sd = nanstd(transition_durations);        % sd for transition durations

%% CORRECT & INCORRECT KEYS

iscorrect_keys = get_keys_info(keys, sequence, n_start_trial);
correct_transition_durations = [];
incoorect_transition_durations = [];
for i_trans = 1 : numel(iscorrect_keys)-1
    if iscorrect_keys(i_trans) && iscorrect_keys(i_trans+1)
        correct_transition_durations = [correct_transition_durations onsets(i_trans+1) - onsets(i_trans)];
    elseif ~iscorrect_keys(i_trans) && ~iscorrect_keys(i_trans+1)
        incoorect_transition_durations = [incoorect_transition_durations onsets(i_trans+1) - onsets(i_trans)];
    end
end
correct_keys.n = numel(find(iscorrect_keys));                   % # of correct keys
correct_keys.mean = nanmean(correct_transition_durations);      % mean transition duration
correct_keys.sd = nanstd(correct_transition_durations);         % sd for transition durations

incorrect_keys.n = numel(find(~iscorrect_keys));              	% # of incorrect keys
incorrect_keys.mean = nanmean(incoorect_transition_durations);	% mean transition duration
incorrect_keys.sd = nanstd(incoorect_transition_durations);  	% sd for transition durations


%% SEQUENCE & ERROR KEYS

trials = get_trials_info(keys, sequence, n_start_trial);
% trials{i}.type
% trials{i}.i_start
% trials{i}.i_end

trials{end+1}.type = 'end';  % add one extra trial at the end;

sequence_keys.n = 0;
sequence_transition_durations = [];

error_keys.n = 0;
error_transition_durations = [];

for i_trial = 1 : numel(trials)-1
    
    check_btwn_seq = 0;
    trial_tmp = trials{i_trial};
    
    switch trial_tmp.type
        
%         case 'head'
%             check_btwn_seq = 1;
            
        case 'sequence'
            check_btwn_seq = 1;
            onsets_tmp = onsets(trial_tmp.i_start:trial_tmp.i_end);
            sequence_keys.n = sequence_keys.n + numel(onsets_tmp);
            for i_trans = 1 : numel(onsets_tmp)-1
                sequence_transition_durations = [sequence_transition_durations onsets_tmp(i_trans+1) - onsets_tmp(i_trans)];
            end
            
        case 'error'
            onsets_tmp = onsets(trial_tmp.i_start:trial_tmp.i_end);
            error_keys.n = error_keys.n + numel(onsets_tmp);
            for i_trans = 1 : numel(onsets_tmp)-1
                error_transition_durations = [error_transition_durations onsets_tmp(i_trans+1) - onsets_tmp(i_trans)];
            end
                        
    end % SWITCH
    
    % transition between sequences  or attempt to perform a sequence
    if check_btwn_seq
        
        trial_next = trials{i_trial+1};
        if ~strcmp(trial_next.type, 'end')
            next_keys_tmp = keys(trial_next.i_start : trial_next.i_end);
            
            if numel(next_keys_tmp) >= n_start_trial && ...
                    ~isempty(strfind(next_keys_tmp(1:n_start_trial), sequence(1:n_start_trial)))
                sequence_keys.n = sequence_keys.n + 1;
                sequence_transition_durations = [sequence_transition_durations onsets(trial_next.i_start) - onsets(trial_tmp.i_end)];
            end
            
        end

    end % IF check for transition betweeen sequences
    
end % FOR each trial

sequence_keys.mean = nanmean(sequence_transition_durations); 	% mean transition duration
sequence_keys.sd = nanstd(sequence_transition_durations);     	% sd for transition durations

error_keys.mean = nanmean(error_transition_durations);      % mean transition duration
error_keys.sd = nanstd(error_transition_durations);     	% sd for transition durations

end

