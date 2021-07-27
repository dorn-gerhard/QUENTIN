function cycle = f_cycle(A, start, path, pos)

n = size(A,1);

index = 1:n;
x = zeros(n,1);
x(pos) = 1;
A_temp = A;
A_temp(path,:) = 0;
y = A_temp * x;

path_new = path;
if pos ~= start
    path_new = [path_new, pos];
end
return_array = cell(0);
index = index(y==1);

for k = 1 : sum(y)
    if index(k) == start
        if length(path) > 0
            cycle = [start, path_new];
            return_array = [return_array, cycle];
        end
    else
        pos_new = index(k);
        return_array = [return_array, f_cycle(A, start, path_new, pos_new)];
    end
end
cycle = return_array;


            