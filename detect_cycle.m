% Detect Cycle in Undirected Graph %
% Handles Multiple Connected Components %
function [found] = detect_cycle(L, n)
    % List of edges L %
    % Number of nodes n %
    % found = 1 if cycle, 0 otherwise %
    
    % Three 'colors' for each node: 
    % 0--Unvisited, 1--Currently Traversing Under, 2--Visited %
    colors = zeros(n, 1); % All nodes are unvisited %
    
    nodes = unique([L(:, 1); L(:, 2)]); % Nodes in L %
    
    v = L(1, 1);
    while (1)
        % Do DFS on every connected component of Graph %
        [colors, found] = visit(L, v, colors);
        
        % If a cycle is detected, break--no need to keep going %
        if (1 == found)
            break;
        end
        
        % If every node in L is hit, end %
        if (0 < min(colors(nodes)))
            break;
        end
        
        % Find disconnected node (new connected component) and repeat %
        I = find(0 == colors(nodes));
        v = nodes(I(1));
    end
end

% Recursive DFS %
function [colors, found] = visit(L, v, colors)
    % List of Edges L %
    % Current Node v %
    % Vector of Node States colors %
    % found: Whether cycle was detected %
    
    % If a node has been visited, there is a cycle %
    if (2 == colors(v))
        found = 1;
        return;
    end
    
    % The current node is being visited %
    colors(v) = 1;
    
    % Visit neighbors that haven't been visited and aren't current %
    N = neighbors(L, v);
    
    for (k = 1:1:length(N)) 
        if (1 ~= colors(N(k)))
            [colors, found] = visit(L, N(k), colors);
            
            % If there is a cycle, we are done %
            if (1 == found)
                return;
            end
        end
    end
    
    colors(v) = 2; found = 0;
end

% Find Neighbor Nodes %
function [I] = neighbors(L, v)
    % List of edges L %
    % Current node v %
    % Neighbor Nodes I %
    
    I1 = find(L(:, 1) == v); I1 = I1(:);
    I1 = L(I1, 2);
    I2 = find(L(:, 2) == v); I2 = I2(:);
    I2 = L(I2, 1);
    I = unique([I1; I2]);
end
