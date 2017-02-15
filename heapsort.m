% a = [10; 3; 8; 1; 12; 2; 5; 4; 5; 7; 15; 14]; %

% Do Heapsort %
function [A] = heapsort(A, d)
    A = A(:);
    [A, heapsize] = build_heap(A, d);
    
    for (k = length(A):(-1):2)
        A = swap(A, 1, k);
        heapsize = heapsize - 1;
        
        A = siftdown(A, 1, d, heapsize);
    end
end

% Build Heap at Start %
function [A, heapsize] = build_heap(A, d)
    heapsize = length(A);
    
    start = ceil((length(A) - 1) / d);
    for (k = start:(-1):1)
        A = siftdown(A, k, d, heapsize);
    end
end

% Recursive function that keeps a heap a heap %
function [A] = siftdown(A, idx, d, heapsize)
    largest = idx;
    
    % Sorting %
    for (k = 1:1:d)
        idx2 = child(idx, k, d);
        if (idx2 <= heapsize)
            if(A(idx2) >= A(largest))
                largest = idx2;
            end
        end
        
    end
    
    % Swap %
    if (idx ~= largest)
        A = swap(A, idx, largest);
        A = siftdown(A, largest, d, heapsize);
    end
end

% Swap two elements of A %
function [B] = swap(A, idx1, idx2)
    B = A;
    temp = B(idx1);
    B(idx1) = B(idx2);
    B(idx2) = temp;
end

% Child Index %
function[k] = child(l, q, d)
    k = d * (l - 1) + (q + 1);
end
