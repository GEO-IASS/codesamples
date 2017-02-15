K = 2; % Clusters %

% % Graph on Wikipedia Page for Kruskal %
% L = [1, 2, 7; 2, 3, 8; 1, 4, 5; 4, 2, 9; 2, 5, 7; 3, 5, 5; 4, 5, 15; 4, 6, 6; 6, 7, 11; 6, 5, 8; 5, 7, 9];
% n = 7;
% [Ls, RSEW, n_list, edges] = MWST_Cluster(L, n, K);
% 
% % Testing Isolated Points %
% L1 = [L; 10, 11, 100; 12 12 0];
% n1 = 12;
% [Ls, RSEW, n_list, edges] = MWST_Cluster(L1, n1, K);

% Generating a Random Graph %
% nr = 50; 
% X = randn(nr, 2);
% Lr = gen_graph(X);
% [Ls, RSEW, n_list, edges] = MWST_Cluster(Lr, nr, K);

% Circles %
nr = 100;
r1 = 10 + rand(nr / 2, 1); 
r0 = 5 + 2 * rand(nr / 2, 1); 
r = [r1; r0]; t = 2.0 * pi * rand(nr, 1); 
x = r .* cos(t); y = r .* sin(t); 
X = [x, y];
Lr = gen_graph(X);
[Ls, RSEW, n_list, edges] = MWST_Cluster(Lr, nr, K);
 
figure(1);
c = 'kbrgycm'; 
for (k = 1:1:size(RSEW, 1))
    index = 1;
    nodes = n_list(index:1:(index + RSEW(k, 2) - 1));
    x = X(nodes, 1); y = X(nodes, 2); 
    plot(x, y, [c(k), '*']); hold on; 
    ylim([-15, 15]);
    xlim([-15, 15]);
end
xlabel('x'); ylabel('y'); 

% %From Slides % 
% Testing Kruskal's Algorithm %
% kruskal(graph0, 6);

% Using MWST_Cluster %
% [Ls, RSEW, n_list, edges] = MWST_Cluster(graph1, 14, 2)

