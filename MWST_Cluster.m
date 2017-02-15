% Cluster Based on MWS (Minimum Weight Spanning) Tree %
function [Ls, RSEW, n_list, edges] = MWST_Cluster(L, n, desired_comp)
    % Depends on kruskal %
    % List of Edges L %
    % Number of Nodes n %
    % Number of Desired Non-Singleton Components %
    % Final Tree Ls, RSEW, List of Nodes, Edges %
    
    Lt = kruskal(L, n); % MWS Tree with edges in increasing order %
    
    % By convention, isolated points have edges to themselves of weight 0 %
    singletons = find(0 == Lt(:, 3));
    if (~isempty(singletons))
        disp(['There are ', num2str(length(singletons)), ' isolated points. They will not count as clusters.']);
    end
    
    % Find connected components and make cuts %
    Ls = Lt;
    count = 1;
    while (1)
        % Count number of non-singleton components %
        [comp, root_size, n_list] = connected_components(Ls, n);
        comp = comp - length(find(1 >= root_size(:, 2)));
        [RSEW, edges] = find_edges(Ls, root_size, n_list);
        
        % No Infinite Loops or Too Many Clusters %
        if (comp >= desired_comp || count >= desired_comp - length(singletons))
            break;
        end
    
        % Make one cut % 
        Ls = split_comp(Ls, RSEW, edges);
        count = count + 1;
        
        if (0 == size(Ls, 1))
            disp(['No Edges Left, Only Dust']);
            break;
        end
    end
end

% Delete One Edge in a Principled Manner %
function [Lt] = split_comp(L, RSEW, edges)
    % List of edges L %
    % RSEW from find_edges %
    % edges from find_edges %
    
    m = size(RSEW, 1);
    
    % Look at every component for one edge to delete %
    index = 1; possible_edges = [];
    for (k = 1:1:m)
        % If there's no weight, only one node (isolated), or only one edge %
        % Don't alter the component %
        if (0 >= RSEW(k, 4) || 1 >= RSEW(k, 2) || 1 >= RSEW(k, 3))
            index = index + RSEW(k, 3);
            continue;
        end
        
        [~, idx] = max(L(edges(index:1:(RSEW(k, 3) - 1)), 3));
        possible_edges = [possible_edges; edges(index + idx - 1)];
        
        index = index + RSEW(k, 3);
    end
    
    [~, idx] = max(L(possible_edges, 3)); 
    idx = possible_edges(idx);
    I1 = (1:1:(idx - 1))'; I2 = ((idx + 1):1:size(L, 1))';
    Lt = L([I1; I2], :);   
end

% Find Edges by Component %
function [RSEW, edges] = find_edges(L, root_size, n_list)
    % List of Edges L %
    % [Root Node, Size of Component] root_size %
    % List of nodes sorted by cluster n_list %
    % [Root Node, Size of Component, Number of Edges in L, Weight] %
    % Edges in L (Indices) edges %
    
    m = size(root_size, 1);
    index = 1; edges = []; len = zeros(m, 1); weight = zeros(m, 1);

    % Loop through components %
    for (k = 1:1:m)
        sz = root_size(k, 2); 
        
        % Loop through nodes in component %
        I = [];
        for (l = 1:1:sz)
            I1 = find(n_list(index + (l - 1)) == L(:, 1)); I1 = I1(:);
            I2 = find(n_list(index + (l - 1)) == L(:, 2)); I2 = I2(:);
            
            I = [I; I1; I2];
        end
        I = unique(I);
        len(k) = length(I);
        edges = [edges; I];
        weight(k) = sum(L(I, 3));
        
        index = index + sz;
    end
    RSEW = [root_size, len, weight];
end

% Find Connected Components in Undirected Graph %
function [comp, root_size, n_list] = connected_components(L, n)
    % List of edges L %
    % Number of nodes n %
    % Number of Connected Components comp %
    % Roots of Connected Components and Size of Connected Components root_size %
    % List of nodes grouped by component n_list %
    
    % Three 'colors' for each node: 
    % 0--Unvisited, 1--Currently Traversing Under, 2--Visited %
    colors = zeros(n, 1);
    
    nodes = unique([L(:, 1); L(:, 2)]); % Nodes in L %
    
    v = L(1, 1); comp = 1; 
    root = []; size_cl = []; n_list = [];
    while (1)
        % Store 'root' of each component %
        root = [root; v];
        
        num_before = length(find(0 == colors(nodes)));
        c_bef = colors;
        
        % Do DFS on every connected component of Graph %
        [colors] = visit(L, v, colors);
        
        % Find size of each component %
        num_after = length(find(0 == colors(nodes)));
        size_cl = [size_cl; num_before - num_after];
        
        % Find nodes in each component %
        I_cl = find(0 ~= colors - c_bef);
        n_list = [n_list; I_cl];
        
        % If every node in L is hit, end %
        if (0 == num_after)
            break;
        end
        
        % Find disconnected node (new connected component) and repeat %
        I = find(0 == colors(nodes));
        v = nodes(I(1));
        comp = comp + 1;
    end
    
    root_size = [root, size_cl];
end

% Recursive DFS %
function [colors] = visit(L, v, colors)
    % List of Edges L %
    % Current Node v %
    % Vector of Node States colors %
    
    % If a node has been visited, there is a cycle %
    if (2 == colors(v))
        return;
    end
    
    % The current node is being visited %
    colors(v) = 1;
    
    % Visit neighbors that haven't been visited and aren't current %
    N = neighbors(L, v);
    
    for (k = 1:1:length(N)) 
        if (1 ~= colors(N(k)))
            [colors] = visit(L, N(k), colors);  
        end
    end
    
    colors(v) = 2; 
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
