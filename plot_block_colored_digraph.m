% Plot function with edge colors for dynamic compensator
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================

function plot_block_colored_digraph(A, n1)

n = size(A,1);
n2 = n - n1;
G = digraph(A');
nodeColors = [repmat([1 0 0], n1, 1);
              repmat([0 0 1], n2, 1)];
[src, tgt] = find(A);
edgeColors = zeros(length(src), 3);

for k = 1:length(src)
    i = src(k);
    j = tgt(k);

    if i <= n1 && j <= n1
        edgeColors(k,:) = [1 0 0]; 
    elseif i <= n1 && j > n1
        edgeColors(k,:) = [0 1 0];
    elseif i > n1 && j <= n1
        edgeColors(k,:) = [0 0 1]; 
    else
        edgeColors(k,:) = [1 0 1];
    end
end

figure;
h = plot(G, 'Layout', 'force', 'EdgeColor', edgeColors, 'NodeColor', nodeColors, ...
         'ArrowSize', 10, 'LineWidth', 1.5, 'EdgeLabel', G.Edges.Weight);

title('Closed-loop system, reduced order compensator');
end