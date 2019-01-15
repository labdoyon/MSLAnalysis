function [sequences, transitions] = get_sequences_info(keys, onsets, sequence, n_start_trial, n_sd)

% keys that are part a sequence at the very beginning or the very end of the block (head and tail, respectively) are not considered
% all vectors are oriented horizontally
% 
% INPUT
% keys                      a vector of keys
% onsets                    a vector of key onsets, i.e., the time that the key was presssed
% sequence      [integer]   a vector of numbers representing the sequence
% n_start_trial [integer]   the number of the first keys of the sequence to search for a trial; the default value is 2 keys
% n_sd          [integer]   the number of standard deviation; is used to exclude outliers
% 
% OUTPUT
% sequences     [cell array]
%   ....duration  	[double]    sequence duration; transition between sequences is not considered
%   ....transitions [double]    a vector of duraitons for the 1st, 2nd, 3rd, ... transitions within a sequence and between the sequences
%
% transitions
%   ...n     	[integer]   a vector of the number of the 1st, 2nd, 3rd, ... transitions within a sequence and between the sequences
%   ...n_out  	[integer]	a vector of the number of the outliers for the 1st, 2nd, 3rd, ... transitions within a sequence and between the sequences
%   ...mean  	[double]  	a vector of mean duration of the 1st, 2nd, 3rd, ... transitions within a sequence and between the sequences
%   ...sd   	[double]  	a vector of sd for the 1st, 2nd, 3rd, ... transitions within a sequence and between the sequences

% Ella Gabitov, 14 January, 2019

if numel(keys) ~= numel(onsets)
    error('The number of keys and the number of onsets do not match.');
end

if nargin < 4, n_start_trial = 2; end;
if isempty(n_start_trial) || isnan(n_start_trial) || n_start_trial == 0, n_start_trial = 2; end;

if nargin < 5, n_sd = 0; end;
if isempty(n_start_trial), n_sd = 0; end;

trials = get_trials_info(keys, sequence, n_start_trial);
% trials{i}.type
% trials{i}.i_start
% trials{i}.i_end

trials{end+1}.type = 'end';  % add one extra trial at the end;

sequences = {};
single_transitions = [];
i_seq = 0;

for i_trial = 1 : numel(trials)-1
    
    transitions_tmp = NaN(1, numel(sequence));
    trial_tmp = trials{i_trial};
    
    if strcmp(trial_tmp.type, 'sequence')
        
        i_seq = i_seq + 1;
        
        onsets_tmp = onsets(trial_tmp.i_start:trial_tmp.i_end);
        sequences{i_seq}.duration = onsets_tmp(end) - onsets_tmp(1);
        
        % transitions within the sequence
        % --------------------------------
        for i_trans = 1 : numel(onsets_tmp)-1
            transitions_tmp(i_trans) = onsets_tmp(i_trans+1) - onsets_tmp(i_trans);
        end
        
        % check for the transition between sequences
        % -------------------------------------------
        trial_next = trials{i_trial+1};
        if ~strcmp(trial_next.type, 'end')
            next_keys_tmp = keys(trial_next.i_start : trial_next.i_end);
            
            if numel(next_keys_tmp) >= n_start_trial && ...
                    ~isempty(strfind(next_keys_tmp(1:n_start_trial), sequence(1:n_start_trial)))
                transitions_tmp(end) = onsets(trial_next.i_start) - onsets(trial_tmp.i_end);
            end
            
        end

        sequences{i_seq}.transitions = transitions_tmp;
        single_transitions(i_seq, :) = transitions_tmp;
        
    end % IF is a sequence
        
end % FOR each trial

transitions = [];
transitions.n = sum(~isnan(single_transitions))';

if ~isempty(n_sd) && n_sd > 0   % remove outliers
    [~,~, data_no_outliers, n_outliers] = remove_outliers(single_transitions, n_sd); 
else
    data_no_outliers = single_transitions;
    n_outliers = zeros(1, numel(sequence));
end 

transitions.n_out = n_outliers;
transitions.mean = nanmean(data_no_outliers);
transitions.sd = nanstd(data_no_outliers);

end


































