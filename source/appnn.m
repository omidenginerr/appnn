function [results,hotspots,peraminoacid] = appnn(sequences)
%NNAP - Predicts the amyloidogenicity propensity of peptides and proteins
%       from the polypeptide sequence alone based on a feed-forward neural
%		network.
%   
% Syntax: 
%
%   [results,hotspots] = nnap(sequences)
%
% Inputs:
%
%   sequences - n*1 cell array of sequences
%
% Outputs:
%
%   results   - n*1 array containing the prediction results as negative (0)
%               or positive (1) for all provided sequences
%   hotspots  - n*1 cell array containing the predicted hotspots for all
%               provided sequences, each cell contains a matrix where each
%               row has two values the first is corresponds to the begining
%               of the hotspot and the second to the end of the hotspot.
%

    % Classify each sequence
    results = zeros(length(sequences),1);
    hotspots = cell(length(sequences),1);
    peraminoacid = cell(length(sequences),1);
    for i = 1:length(sequences)
        rawdata = classify(sequences{i});
        results(i) = max(rawdata);%ge(sum(ge(sqnclssfd,0.5)),1);
        hotspots{i} = stretches(rawdata);
        peraminoacid{i} = paaprediction(rawdata);
    end
end

function [hotspots,sqncrp] = classify(sequence)
    % Encode the sequence
    sqncrp = encode(sequence);

    % Load neural network
    load('data.mat','nn')
    
    % Classify each six amino acids window
    hotspots= nn(sqncrp)';
end

function [encodedseq] = encode(sequence)
    load('data.mat','aaindex','mapping')
    
    % Set lengths
    ftslen = length(mapping.fts);
    wndslen = length(sequence)-5;
    
    % Set output array
    encodedseq = zeros(ftslen,wndslen);
    % For each six amino acids window (i is the index of the first aa)
    for i = 1:wndslen
        windowtemp = zeros(ftslen,1);
        % For each feature to be computed (j is the feature index)
        for j = 1:ftslen
            temp = zeros(6,1);
            % For each amino acid in the window
            for h = 1:6
                % Get amino acid feature value
                index = find(strcmp(sequence(h+i-1),aaindex.aaindex));
                temp(h) = aaindex.values(index,mapping.fts(j));
            end
            % Compute window feature through the appropriate funtion
            windowtemp(j) = mapping.functions{j}(temp);
        end
        % Map each window accordingly to the original std mapping
        encodedseq(:,i) = mapstd('apply',windowtemp,mapping.std);
    end
end

function [hotspots] = stretches(rawdata)
    % Define results variable
    hotspots = [];
    
    % Find hotspots
    j = 1;
    while j <= length(rawdata)
        if rawdata(j) > 0.5
           a = j;
           while and(j<length(rawdata),rawdata(j) > 0.5)
               j = j+1;
           end
           hotspots = vertcat(hotspots,[a,j+5]);
        end
        j = j+1;
    end
end

function [paapv] = paaprediction(rawdata)
    % Define results variable
    paapv = [];
    
    % Extract per amino acid prediction values
    for i=1:(length(rawdata)+5)
        j = i-5;
        if i < 6
            j = 1;    
        end
        h = i;
        if i > length(rawdata)
            h = length(rawdata);
        end
        paapv = vertcat(paapv,max(rawdata(j:h)));
    end
end