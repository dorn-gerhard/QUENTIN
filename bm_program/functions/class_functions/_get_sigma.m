function [Sig, status] = get_sigma(K, Tab_cut, status, mes_flag, Energy_cut, Energy_cut_2)
if isempty(mes_flag), mes_flag = false; end

diag_index_shortened = Tab_cut.Diagonal_Index();    

%procedure to extract steady state from K

%general idea: if some criteria is not fullified (not converged) increase Energy_cut

%first: extract eigenvector around zero

if  size(K,1) <= 500
    [ev,ew] = eig(full(K));
    [evl,ewl] = eig(full(K.'));
else
    %alternative is to use contour integration method or time evolution method
    [ev,ew, flag] = eigs(K,10,10^-10);
    [evl,ewl] = eigs(K.', 10, 10^-10);
    if flag ~= 0
        warning(['Flag when using eigs: ', num2str(flag)]) 
    end
    
    
    
    if any(isnan(ev(:,:)))
        error('EIGS not converged')
    end
    
    
    
    
end

%second: count, how many there are
%cases: 
%   number not equal
%   numbers equal, but more than 1
%   one left and right eigenvalue meeting the requirement

eww = diag(ew); 
ewwl = diag(ewl); 

[Lew, Lewl, n_right, n_left, status] = func_check_tolerance(eww, ewwl, status);
%guarantees, that number of eigenvalues (left and right) are identical

if any(real(eww) > 10^-8), disp('WARNING: positive Eigenvalues...'); status.positive_eigenvalues = 1; 
else
    status.positive_eigenvalues = 0;
end

status.no_solution = 0;
if n_left == 0
    
    Sig = [];
    
    disp('WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NO SOLUTION!!!!!!!!!!!!!!!!!!!!!!!!');
    disp(['minimal value: ', num2str(eww(abs(eww) == min(abs(eww))).')])
    status.minimal_value = eww(abs(eww) == min(abs(eww)));
    
    status.no_solution = 1;
    status.multiple_solutions = 0;
    status.error_text = ['no eigenvalues found with tolerance:', num2str(status.eigenvalue_tolerance), ', minimal value: ', num2str(eww(abs(eww) == min(abs(eww))).')];
    return
end


if n_right > 1, disp('WARNING: multiple solutions - NEW PROCEDURE'); 
    status.multiple_solutions = sum(Lew);
    
    %TODO: three methods, 
    % a) take minimal values (check if there are more singular values than others
    % b) check projection method
    % c) time evolution
   
    %a) analyse eigenvalues
    % deal with situation where one eigenvalue is 5 orders smaller:
    order_of_ew = log(abs(eww(Lew)))/log(10);
     
    order_of_ewl = log(abs(ewwl(Lewl)))/log(10);
    
    Lew_select =  max(order_of_ew)- order_of_ew   > 5;
    Lewl_select = max(order_of_ewl)- order_of_ewl   > 5;
    if (sum(Lew_select) < sum(Lew)) && (sum(Lewl_select) < sum(Lewl)) && (sum(Lew_select) >= 1) && (sum(Lewl_select) >= 1) && (sum(Lewl_select) == sum(Lew_select))
        disp('smaller eigenvalues taken (4 orders of magnitude smaller)')
        
        Lew(Lew) = Lew_select;
        Lewl(Lewl) = Lewl_select;
        
    end
        

    
    %b) projection method:
    disp('Apply projection + time evolution method!')
   
    
    er = ev(:,Lew);
    el = evl(:,Lewl).';
    for perm_ind = 1:sum(Lew)
        if abs((el(perm_ind,:) * er(:,1)) ) > 10^-12 %vectors are already orthogonal
            el([1,perm_ind],:) = el([perm_ind,1],:); % swaps some vectors - not necessary
            break;
        else
            warning('some eigenvectors are already orthogonal')
        end
    end
    el(1,:) = el(1,:) ./ (el(1,:)*er(:,1));
   
    

        
    for k = 2:sum(Lew)  %Algorithm to gain inverse eigenvectors in T^-1 (D = T^-1 * K * T) by mutual orthogonalisation
        
        
        el(k,:) = el(k,:) - sum(repmat((el(k,:) * er(:,1:k-1)).', 1, length(K)) .* el(1:k-1,:),1);
        if abs((el(k,:) * er(:,k)) ) > 10^-14
            el(k,:) = el(k,:) ./ (el(k,:)*er(:,k));
        else
            warning('some eigenvectors are already orthogonal this shouldnot happen - repair in get_sigma')
        end
        er(:,k) = er(:,k) - sum(repmat((el(1:k-1,:) * er(:,k)).' ,length(K), 1) .* er(:,1:k-1),2);
        
    end
    
    if mes_flag
        disp(['quality of eigenvectors |K*v_r|, |v_l*K: ', num2str([sqrt(sum((abs(K*er)).^2,1)), sqrt(sum((abs(el*K)).^2,2)).'])])
        disp(['orthogonality check: ', num2str(sum(sum(abs(el*er-eye(sum(Lew))))))])
    end
    
    A = er * el;
    % ---- creation of startvector, either random or equally distributed or groundstate like
    
    %TODO: make a parameter to set up
    %version a) random state
    %initial_sigma = Sigma(En, [], Tab_cut, Energy_cut, Energy_cut_2);
    %vec= initial_sigma.f_vector();
    
    %version b) maximum entangled state aka Linksvakuum
    max_ent_state = Tab_cut.f_vector(eye(Tab_cut.Numel_Datasets)/Tab_cut.Numel_Datasets);
    initial_sigma = Sigma([], max_ent_state , Tab_cut); % check if Energy_cut is needed
    vec = initial_sigma.f_vector();
    
    Sig_vec = A * vec;
    Sig_vec(abs(Sig_vec) < Tab_cut.Tolerance) = 0; %denoise
    
    
    
    
    %update Lew vector

else
        ew = eww(Lew);  %disp(ew);
        Sig_vec = ev(:,Lew);
        %denoise
        Sig_vec(abs(Sig_vec) < Tab_cut.Tolerance) = 0;
        status.multiple_solutions = 1;
end
       
       
% ------------------------- analyse the received Sigma vector -------------------------------
%Sig_vec_full = sparse(sigma_index,ones(size(sigma_index)), Sig_vec, sum(Tab.Block_Size.^2),1);
%Sig_diag = Sig_vec_full(Tab.Diagonal_Index());

Sig_diag = Sig_vec(diag_index_shortened);


Sig_trace = sum(Sig_diag);
if abs(Sig_trace) < 10^-16
    warning('Sigma trace is too small1!!')
end

%rough correction, a first try
%first correct for complex trace
%second compensate for absolute value

trace_argument = atan2(imag(Sig_trace), real(Sig_trace));
Sig_vec = Sig_vec * exp(-1i*trace_argument) /abs(Sig_trace); %correction by scalar
%denoise
Sig_vec(abs(Sig_vec) < Tab_cut.Tolerance) = 0;

Sig_diag = Sig_vec(diag_index_shortened);
Sig_trace = sum(Sig_diag);




% check argument of diagonal sigma (diagonal should not be 
status.not_same_argument = 0;
argument = atan2(imag(Sig_diag), real(Sig_diag));

argument(abs(Sig_diag) < 10^-12) = 0; % to ignore rounding errors of infinitesimal blocks
argument(abs(argument) < 10^-12) = 0; % to avoid +- eps around 0 / pi after mod function

argument = mod(argument, pi);
L_argument = abs(Sig_diag) > 10^-10;
if var(argument(L_argument)) > 10^-15
    disp('WARNING: additional rounding for argument of Sigma diagonal')
    %second step
    L_argument_2 = abs(Sig_diag) > 10^-8;
else
    L_argument_2 = L_argument;
end
if var(argument(L_argument_2)) > 10^-13
    disp(['WARNING: argument of Sigma not well defined - diagonal not real, var: ', num2str(var(mod(atan2(imag(Sig_diag), real(Sig_diag)), pi)))]);
    %TODO: Think about better way and criterion when to refine procedure, check status to identify problematic
    %matrices
    status.not_same_argument = 1;
    % procedure to repair that effect for L_blocks
    
    arg_block = Tab_cut.f_block(sparse(diag(argument)));
    Sig_diag_block = Tab_cut.f_block(sparse(diag(Sig_diag)));
    L_block = cellfun(@(x) trace(abs(x)) > Tab_cut.Tolerance, Sig_diag_block);
    
    %cellfun(@(x) vec(x), Tab.Diag_Index, 'Uniformoutput', false)
    if sum(cellfun(@(x) var(full(diag(x))), arg_block(L_block))) > 10^-8
        disp(['WARNING: individual blocks do not have same complex argument!!!: ', num2str(sum(cellfun(@(x) var(full(diag(x))), arg_block(L_block))))])
    end
    correction = cellfun(@(x) exp(-1i*mean(diag(x))), arg_block,'Uniformoutput', false);
    %TODO: implement a routine that does not run out of memory!
    block = cellfun(@(x,y) x*y, Tab_cut.f_block(Sig_vec), correction, 'Uniformoutput', false);
    
    %Sig = Sigma(En, block, Tab_cut, Energy_cut, Energy_cut_2);
    Sig = Sigma([],Tab_cut.f_block(Sig_vec./Sig_trace), Tab_cut, Energy_cut, Energy_cut_2);
else 
    Sig = Sigma([],Tab_cut.f_block(Sig_vec./Sig_trace), Tab_cut, Energy_cut, Energy_cut_2);
end

trace(Sig.f_matrix)
Sig.s_dom_blocks




%TODO: rethink analysis, and Blockstructure
%Sig.f_validate;
%Sig.Analysis
%Sig.Block_structure = En.Structure;
%Sig.s_plot

end




function [Lew, Lewl, n_right, n_left, status] = func_check_tolerance(eww, ewwl, status)
Lew = abs(eww) < status.eigenvalue_tolerance;
Lewl = abs(ewwl) < status.eigenvalue_tolerance;
n_right = sum(Lew);
n_left = sum(Lewl);

if n_right ~= n_left
    disp('# of LEFT and RIGHT EV not equal - refine eigenvalue_tolerance')
    status.eigenvalue_tolerance = status.eigenvalue_tolerance/100;
    [Lew, Lewl, n_right, n_left, status] = func_check_tolerance(eww, ewwl, status);
else
    if n_right > 0
        status.eigenvalue_tolerance = mean([abs(eww(Lew)); abs(ewwl(Lewl))]); 
    end
end
end
