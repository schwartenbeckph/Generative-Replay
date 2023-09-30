% obtain best hyperparameters from decoding and plot accuracy

function [L1s, L1_mean, TSs, TS_mean] = ObtaindecodingAcc(decode,plot_decoding)

    if nargin==0
        decode = 'localiser'; % can also be 'relation' to decode relations
%         decode = 'relation';
        plot_decoding = true;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETUP THE MATLAB PATHS and variables
    % this will be your base directory
    based = '/Users/Philipp/Dropbox/StimConstrMEG/Results_Scanning/code_MEG_Final';

    scan_result_path    = fullfile(based,'data');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if plot_decoding    
        
        if strcmp(decode,'localiser')
            
            [~,element_accuracy,~,L1_element,ts_best_elements_subj,L1_element_mean]   = AccuracyDecoding_allTS(scan_result_path,decode,true);
            
            [~,ts_best_elements] = max(element_accuracy); % this is a time step where decoding of elements peaks
            
            fprintf('Peak element decoding is at %dms.\n',(ts_best_elements-1)*10) % -1 because we start at 0!

        elseif strcmp(decode,'relation')
            
            [~,relation_accuracy,~,L1_relation,ts_best_relations_subj,L1_relation_mean]   = AccuracyDecoding_allTS(scan_result_path,decode,true);
            
            [~,ts_best_relations] = max(relation_accuracy); % this is a time step where decoding of relation peaks
            
            fprintf('Peak relation decoding is at %dms.\n',(ts_best_relations-1)*10)

        end
        
    else % just obtain best decoding stats
        
        if strcmp(decode,'localiser')
            load('DecodingParams_element.mat')
        elseif strcmp(decode,'relation')
            load('DecodingParams_relation.mat')
        end
        
    end
    
    if strcmp(decode,'localiser')
        
        L1s     = L1_element;
        L1_mean = L1_element_mean;
        TSs     = ts_best_elements_subj;
        TS_mean = ts_best_elements;
        
    elseif strcmp(decode,'relation')
        
        L1s     = L1_relation;
        L1_mean = L1_relation_mean;
        TSs     = ts_best_relations_subj;
        TS_mean = ts_best_relations;
        
    end    


end