function cnt = nine_digcnt(k)
%nine_digcnt Logarithmic scale which counts nines
% ex: 0.9  |-> 1
%     0.99 |-> 2
cnt = log(1-9*(k/9))/log(1/10);
end