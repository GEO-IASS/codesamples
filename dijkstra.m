function [dist, prev] = dijkstra(A, start)
    % Node arc A, start node %
    
    dist = (size(A, 1))^3 * ones(size(A, 1), 1); 
    prev = (size(A, 1))^3 * ones(size(A, 1), 1);
    
    dist(start) = 0; prev(start) = 0;
    
    Q = (1:1:size(A, 1))';
    QQ = [Q, dist];
    
    while (1)
        [~, idx] = min(QQ(:, 2));
        idx_n = QQ(idx, 1); 
        
        % FOR UNDIRECTED GRAPH %
        % arc_idx = find(abs(A(idx_n, :)) > 0.5); % Incoming neighbors: arc indices %
        
        % FOR DIRECTED GRAPH %
        arc_idx = find(A(idx_n, :) > 0.5); % Incoming neighbors: arc indices %
        
        for (l = 1:1:length(arc_idx))
            idx2 = arc_idx(l);
            
            % Undirected %
            % I = find(abs(A(:, idx2)) > 0.5);
            % idx2 = setdiff(I, [idx_n]);
            
            % Directed %
            idx2 = find(A(:, idx2) < -0.5);
            
            if (dist(idx_n) + 1 < dist(idx2))
                dist(idx2) = dist(idx_n) + 1;
                prev(idx2) = idx_n;
            end
        end
        
        if (size(QQ, 1) < 2)
            break; 
        end
        QQ = [QQ(1:1:(idx - 1), :); QQ((idx + 1):1:end, :)];
    end
end
