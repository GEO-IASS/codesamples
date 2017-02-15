% Greedy Algorithm to Build Minimal Weight Spanning Tree %
function [Lt] = kruskal(L, n)
    % Depends on detect_cycle %
    % List of edges L (First, Second, Positive Weight) %
    % Number of Nodes n %
    % List of edges Lt forming MWS Tree %
    
    % Check for appropriate number of nodes %
    nodes = unique([L(:, 1); L(:, 2)]); % Nodes in L %
    if (n < max(nodes))
        disp(['Number of Nodes is Insufficient: Correcting']);
        n = max(nodes);
    end
    
    % Find Isolated Nodes %
    lone_nodes = setdiff((1:1:n)', nodes); % Isolated Nodes %
    if (~isempty(lone_nodes))
        disp(['There are ', num2str(length(lone_nodes)), ' unaccounted-for isolated nodes: ', num2str(lone_nodes'), '; Correcting']);
        
        % Add Isolated Nodes to Graph (loops with no weight) %
        for (k = 1:1:length(lone_nodes))
            L = [L; lone_nodes(k), lone_nodes(k), 0];
        end
    end
    
    L = sortrows(L, 3); % Sort Edges from Smallest Weight to Largest %
    
    Lt = []; index = 1; len = size(L, 1); 
    while (1)  
        
        % Add edges to tree as long as no cycle is found %
        Ltemp = [Lt; L(index, :)];
        % detect_cycle finds whether a cycle would be formed
        % if the arc was added to the tree% 
        found = detect_cycle(Ltemp, n);
        if (0 == found)
            Lt = Ltemp;
        end
        
        % Move on to next edge %
        index = index + 1;
        if (index > len)
            break; 
        end
    end
    
    Lt = sortrows(Lt, 3); % Just in case %
end
