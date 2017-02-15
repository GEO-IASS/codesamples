% Generate a List of Edges from Data %
function [L] = gen_graph(X)
    % Data X, Rows are points / Values %
    % List of Edges L %
    
    n = size(X, 1);
    L = zeros((n * n - n) / 2, 3); % n choose 2 edges %
    index = 1;
    
    % Every node is connected to every other node %
    for (k = 1:1:n)
        for (l = (k + 1):1:n)
            % For Gaussian similarity: Weights are (1 / Similarity): Dissimilarity %
            % Thus, L(index, :) = [k, l, 1.0 ./ similarity(X(k, :), X(l, :))]; %
            % Instead, here we are using the Euclidean distance %
            L(index, :) = [k, l, dissimilarity2(X(k, :), X(l, :))];
            index = index + 1;
        end
    end
end

% Similarity Function %
function [s] = similarity(a, b)
    % Two points a, b %
    % Returns Gaussian Similarity s %
    
    sig_inv2 = 1.0 / (0.50 * 0.50); % sigma = 0.5; %
    
    % exp{-(1/2) * ||a - b||_2^2 / \sigma^2} %
    s = norm(a - b, 2);
    s = exp(-0.5 * sig_inv2 * s * s); 
end

% Dissimilarity Function %
function [s] = dissimilarity2(a, b)
    % Two points a, b %
    % Returns Euclidean Dissimilarity s %
    
    s = norm(a - b, 2); 
end
