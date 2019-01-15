function [perf_duration,...
            n_keys,...
            sequences,...
            btwn_seq,...
            errors...
            ]...
            = get_trials_stat(keys, onsets, sequence, n_start_trial)

% keys that are part a sequence at the very beginning or the very end of the block (head and tail, respectively) are not considered
% 
% INPUT
% keys                      a vector of keys
% onsets                    a vector of key onsets, i.e., the time that the key was presssed
% sequence      [integer]   a vector of numbers representing the sequence
% n_start_trial [integer]   the number of the first keys of the sequence to search for a trial; the default value is 2 keys
% 
% OUTPUT
% perf_duration  	[integer]    time between the first and the last presses key 
% n_keys            [integer]
%
% sequences     
%   ....n       [integer]
%   ....mean  	[double]    transitions before or between sequences are not considered
%   ....sd   	[double]
%
% btwn_seq     
%   ....n       [integer]
%   ....mean  	[double]    only transitions between sequences or attempt to perform a sequence;
%                           such attempt is identified by the <n_start_trial> first keys of the sequence
%   ....sd   	[double]
%
% errors
%   ....n       [integer]
%   ....mean  	[double]    transitions before or between errors are not considered
%   ....sd   	[double]
%

% Ella Gabitov, 14 January, 2019

if numel(keys) ~= numel(onsets)
    error('The number of keys and the number of onsets do not match.');
end

if nargin < 4, n_start_trial = 2; end;
if isempty(n_start_trial) || isnan(n_start_trial) || n_start_trial == 0, n_start_trial = 2; end;

perf_duration = onsets(end) - onsets(1);
n_keys = numel(keys);

trials = get_trials_info(keys, sequence, n_start_trial);
% trials{i}.type
% trials{i}.i_start
% trials{i}.i_end

trials{end+1}.type = 'end';  % add one extra trial at the end;

sequence_durations = [];
btwn_seq_durations = [];
error_durations = [];

for i_trial = 1 : numel(trials)-1
    
    check_btwn_seq = 0;
    trial_tmp = trials{i_trial};
    
    switch trial_tmp.type
        
%         case 'head'
%             check_btwn_seq = 1;
            
        case 'sequence'
            check_btwn_seq = 1;
            dur_tmp = onsets(trial_tmp.i_end) - onsets(trial_tmp.i_start);
            sequence_durations = [sequence_durations dur_tmp];
            
        case 'error'
            dur_tmp = onsets(trial_tmp.i_end) - onsets(trial_tmp.i_start);
            error_durations = [error_durations dur_tmp];
                        
    end % SWITCH
    
    % transition between sequences  or attempt to perform a sequence
    if check_btwn_seq
        
        trial_next = trials{i_trial+1};
        if ~strcmp(trial_next.type, 'end')
            next_keys_tmp = keys(trial_next.i_start : trial_next.i_end);
            
            if numel(next_keys_tmp) >= n_start_trial && ...
                    ~isempty(strfind(next_keys_tmp(1:n_start_trial), sequence(1:n_start_trial)))
                dur_tmp = onsets(trial_next.i_start) - onsets(trial_tmp.i_end);
                btwn_seq_durations = [btwn_seq_durations dur_tmp];
            end
            
        end

    end % IF check for transition betweeen sequences
    
end % FOR each trial

sequences.n = numel(sequence_durations);        % # of sequences
sequences.mean = nanmean(sequence_durations); 	% mean sequence duration
sequences.sd = nanstd(sequence_durations);      % sd for sequences

errors.n = numel(error_durations);          % # of errors
errors.mean = nanmean(error_durations); 	% mean error duration
errors.sd = nanstd(error_durations);        % sd for errors

btwn_seq.n = numel(btwn_seq_durations);     	% # of btwn_seq transitions
btwn_seq.mean = nanmean(btwn_seq_durations); 	% mean btwn_seq duration
btwn_seq.sd = nanstd(btwn_seq_durations);       % sd for btwn_seq transitions

end
































