function [Sig, status] = get_sigma(K, En, Tab_cut, status, mes_flag, Energy_cut, Energy_cut_2, sigma_tolerance)
if isempty(mes_flag), mes_flag = false; end

if isa(sigma_tolerance, 'Energy')
    error('replace energy_input in get_sigma call by sigma_tolerance')
end
    
diag_index_shortened = Tab_cut.Diagonal_Index();    

%procedure to extract steady state from K

%general idea: if some criteria is not fullified (not converged) increase Energy_cut

%first: extract eigenvector around zero

if  size(K,1) <= 4100
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
% TODO: the case not equal number of eigenvalues should be better monitored!!!

if any(real(eww) > 10^-8), warning(['WARNING: positive Eigenvalues... MAX: ', num2str(max(real(eww))), ', sum: ', num2str(sum(real(eww(real(eww)>10^-8))))]); 
    status.positive_eigenvalues = max(real(eww)); 
else
    status.positive_eigenvalues = 0;
end

status.no_solution = 0;
if n_left == 0 %is equal to n_right = 0 -> see func_check_tolerance 
    warning('in get_sigma: n_left = 0')
    Sig = [];
    
    
    status.minimal_value = min(eww(abs(eww) == min(abs(eww))));
    
    status.no_solution = 1;
    status.multiple_solutions = 0;
    status.error_text = ['no eigenvalues found with tolerance:', num2str(status.eigenvalue_tolerance), ', minimal value: ', num2str(eww(abs(eww) == min(abs(eww))).')];
    if mes_flag
        disp('WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NO SOLUTION!!!!!!!!!!!!!!!!!!!!!!!!');
        disp(status.error_text)
    end
    return
end


if n_right > 1
    warning('WARNING: multiple solutions - NEW PROCEDURE');
    status.multiple_solutions = sum(Lew);
    
    %TODO: three methods,
    % a) take minimal values (check if there are more singular values than others
    % b) check projection method
    % c) time evolution
    % d) check by corresponding energies!
    
    
    
    
    %=================================================================================
    %d) check average energy  (Problem with average energy if allowed gap is too
    %high!!!
    if false
        ew_index = find(Lew);
        if ismember('A', Tab_cut.Order)
            energies = cell2mat(arrayfun(@(x,y) repmat(x, y,1), Tab_cut.Average_Energies, Tab_cut.Block_Size.^2, 'uniformoutput', false));
        else
            energies = cell2mat(arrayfun(@(x,y) repmat(x, y,1), Tab_cut.Energies, Tab_cut.Block_Size.^2, 'uniformoutput', false));
        end
        E = zeros(numel(ew_index),1);
        for perm_ind = 1:numel(ew_index)
            er = ev(:,ew_index(perm_ind));
            E(perm_ind,1) = sum(abs(er).^2.*energies);
        end
        
        %d2) check energy
        ew_index = find(Lew);
        energies = En.Energies(Tab_cut.List_Index);
        E = zeros(numel(ew_index),1);
        for perm_ind = 1:numel(ew_index)
            er = ev(Tab_cut.Diagonal_Index,ew_index(perm_ind));
            E(perm_ind,1) = sum(abs(er).^2.*energies);
        end
        
        %NOTE: first attempt: exclude highest 10 % - now, take lowest energies + 20% tolerance!
        % Lew_reas_en = E < min(energies) + 0.9*(max(energies) - min(energies));
        Lew_reas_en = E < min(E) + 0.5*(max(max(energies), Energy_cut_2) - min(energies));
        
        
        if (sum(Lew_reas_en) < sum(Lew)) && (sum(Lew_reas_en) >= 1)
            Lew(Lew) = Lew_reas_en;
            Lewl(Lewl) = Lew_reas_en;
            disp(['eigenvectors with average eigenenergy of ', num2str(mean(E(~Lew_reas_en))), ' discarded!'])
        end
    end
    %=================================================================================
    
    %=================================================================================
    %a) analyse eigenvalues
    % deal with situation where one eigenvalue is 5 orders smaller:
    if false
        order_of_ew = log(abs(eww(Lew)))/log(10);
        
        order_of_ewl = log(abs(ewwl(Lewl)))/log(10);
        
        Lew_select =  max(order_of_ew)- order_of_ew   > 5;
        Lewl_select = max(order_of_ewl)- order_of_ewl   > 5;
        if (sum(Lew_select) < sum(Lew)) && (sum(Lewl_select) < sum(Lewl)) && (sum(Lew_select) >= 1) && (sum(Lewl_select) >= 1) && (sum(Lewl_select) == sum(Lew_select))
            disp('smaller eigenvalues taken (4 orders of magnitude smaller)')
            
            Lew(Lew) = Lew_select;
            Lewl(Lewl) = Lewl_select;
            
        end
        Sig_vec = ev(:,Lew);
    end
    %=================================================================================
    
    %=================================================================================
    %b) projection method:
    if false
        if sum(Lew) > 1
            
            if mes_flag
                disp('Apply projection + time evolution method!')
            end
            
            
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
                er(:,k) = er(:,k)./norm(er(:,k));
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
            Sig_vec(abs(Sig_vec) < sigma_tolerance) = 0; %denoise
        else
            Sig_vec = ev(:,Lew);
        end
    end
    
    %=================================================================================
    %c) time evolution
    beta_temp = 20;
    weighted_energies = exp(-beta_temp*(En.Energies(Tab_cut.Old_List_Index) - min(En.Energies(Tab_cut.Old_List_Index))));
    weighted_energies = denoise(weighted_energies/sum(weighted_energies), 10^-5);
    %rand_rho = Sigma([],[],Tab_cut); %
    rand_rho = Sigma([],sparse(diag(weighted_energies)),Tab_cut);
    Sig_vec = expm(K*10^4)*rand_rho.f_vector;
    converged_flag = true;
    if sum(abs(K*Sig_vec)) >= status.eigenvalue_tolerance || any(isnan(Sig_vec))
        warning('time evolution for t=10^4 not successful (no steady state)')
        Sig_vec = expm(K*10^8)*rand_rho.f_vector;
        if sum(abs(K*Sig_vec)) >= status.eigenvalue_tolerance || any(isnan(Sig_vec))
            warning('time evolution for t=10^8 not successful (no steady state)')
            Sig_vec = expm(K*10^12)*rand_rho.f_vector;
            if sum(abs(K*Sig_vec)) >= status.eigenvalue_tolerance || any(isnan(Sig_vec))
                warning('no convergend eigenstate with time evolution (t = 10^12) found!')
                converged_flag = false;
            end
        end
    end
    Sig_vec(abs(Sig_vec) < sigma_tolerance) = 0;
    
    if converged_flag == false
        %beta_temp = 20;
        %weighted_energies = exp(-beta_temp*En.Energies(Tab_cut.List_Index));
        %weighted_energies = denoise(weighted_energies/sum(weighted_energies), 10^-5);
        random_matrix = Sigma([],[], Tab_cut);
        rand_rho = Sigma([], 1/2*(random_matrix.f_matrix + Tab_cut.f_matrix(sum(ev(:,Lew),2)/sum(Lew))),Tab_cut);
        Sig_vec = expm(K*10^4)*rand_rho.f_vector;
        if sum(abs(K*Sig_vec)) >= status.eigenvalue_tolerance || any(isnan(Sig_vec))
            warning('time evolution for t=10^4 not successful (no steady state)')
            Sig_vec = expm(K*10^8)*rand_rho.f_vector;
            if sum(abs(K*Sig_vec)) >= status.eigenvalue_tolerance || any(isnan(Sig_vec))
                warning('time evolution for t=10^8 not successful (no steady state)')
                Sig_vec = expm(K*10^12)*rand_rho.f_vector;
                if sum(abs(K*Sig_vec)) >= status.eigenvalue_tolerance || any(isnan(Sig_vec))
                    error('no convergend eigenstate with time evolution (t = 10^12) found!')
                end
            end
        end
    end

    

    
    %update Lew vector
    
else
    ew = eww(Lew);  %disp(ew);
    Sig_vec = ev(:,Lew);
    %denoise
    Sig_vec(abs(Sig_vec) < sigma_tolerance) = 0;
    status.multiple_solutions = 1;
end
       
%alternative way to extract Hermitian density matrix
%NOTE: Matrix must always be Hermitian - otherwise big error 

Sig_vec_denoise = denoise(Sig_vec./norm(Sig_vec), sigma_tolerance);

Values = Tab_cut.f_block(Sig_vec_denoise);

block_weight = cellfun(@(x) sum(abs(full(x(:)))), Values);
L_ind = find(block_weight > sigma_tolerance);

weight = cell2mat(cellfun(@(x,y) repmat(sum(abs(full(x(:)))),y,1), Values(L_ind), num2cell(Tab_cut.Block_Size(L_ind)), 'uniformoutput', false));
argument = full(cell2mat(cellfun(@(x) atan2(imag(diag(x)), real(diag(x))), Values(L_ind), 'uniformoutput', false)));
argument = mod(argument, 2*pi);
if mean(argument) < pi % in order to avoid jumps!!!
    argument(argument > 5.5) = argument(argument > 5.5) - 2*pi;
else
    argument(argument < 1) = argument(argument < 1) + 2*pi;
end
weight = cell2mat( cellfun(@(x) diag(full(x)), Values(L_ind) ,'UniformOutput', false));
% if matrix is truely hermitian then only two values are possible!
relevant_argument = argument(weight > 10^-2);

[unique_phase,~,frequency_of_phase] = uniquetol([relevant_argument;7], 10^-1);
%NOTE: uniquetol uses a relative tolerance to the maximal value, since the
%interval is between [0,2pi] 7 is added and later removed to get a fixed
%absolute error tolerance
unique_phase = unique_phase(1:end-1); frequency_of_phase = frequency_of_phase(1:end-1);

phase_correction = mean(relevant_argument(frequency_of_phase == mode(frequency_of_phase)));

if numel(unique_phase) > 2%var(argument, weight)^2 > sigma_tolerance
    error('density matrix has not the same complex phase on all diagonal elements!')
    %NOTE: density matrix should be Hermitian, but could have negative
    %eigenvalues - therefore discrimination
end

if numel(unique_phase) == 2
    warning(['Density matrix is either negative or not Hermitian!!!'])
end

Sig_vec_denoise = denoise(Sig_vec_denoise * exp(-1i*phase_correction), sigma_tolerance);

S = Tab_cut.f_matrix(Sig_vec_denoise);
% NOTE: trace norm always 1 - not equal to trace(abs(rho))!!!
trace_abs = trace(sqrtm(full(S'*S))); 
Sig_vec_denoise = denoise(Sig_vec_denoise./trace_abs, sigma_tolerance);
Sig = Sigma([], Sig_vec_denoise, Tab_cut, Energy_cut, Energy_cut_2);
dom_blocks = Sig.s_dom_blocks;
if mes_flag
    disp(dom_blocks)
    disp(Sig.Analysis)
end

% A = cell2mat(cellfun(@(x) diag(x), Values, 'uniformoutput', false)); %NOTE: one could also use trace(Tab_cut.f_matrix(Values))
% total_trace = sum(A);
% Values = cellfun(@(x) denoise(1/total_trace*x, sigma_tolerance), Values, 'uniformoutput', false);
% 
% Sig = Sigma([],Values, Tab_cut, Energy_cut, Energy_cut_2);



%raw_trace = sum(cellfun(@(x) trace(x), Values(L_ind)));


if false
    show_warning_trace_once = true;
    warning_neg_ev_val = 0;
    %start to correct all dominant blocks individually!
    for k_ind = 1:numel(L_ind)
        k = L_ind(k_ind);
        Values{k} = 1/2 * full((Values{k} + Values{k}'));
        if any(diag(Values{k}) < -10^-10)
            
            if all(diag(Values{k}) < -10^-10)
                Values{k} = -(Values{k});
            else
                if show_warning_trace_once
                    warning('problem: some negative trace elements in sigma!')
                    show_warning_trace_once = false;
                end
            end
        end
        
        %check positivity, start with Gershgorin
        if any(diag(Values{k}) - sum(abs(Values{k} - diag(diag(Values{k}))),2) < -10^-10)
            [evec_block, eval_block] = eig(Values{k});
            
            if any(diag(eval_block) < -10^-10)
                eigval = diag(eval_block); %NOTE: do not overwrite eval - its a function in Matlab
                warning_neg_ev_val = warning_neg_ev_val + sum(abs(eigval(eigval <-10^-10)));
                
                Values{k} = evec_block * abs(eval_block) * evec_block';
            end
        end
        
    end
    
    if warning_neg_ev_val
        warning(['problem: negative eigenvalues in sigma, artificially repaired to valid density matrix! Sum: ', num2str(warning_neg_ev_val)])
    end
end



       
% % ------------------------- analyse the received Sigma vector if it is a correct density matrix -------------------------------
% %Aim: produce a correct density matrix, if not possible force it!
% %Sig_vec_full = sparse(sigma_index,ones(size(sigma_index)), Sig_vec, sum(Tab.Block_Size.^2),1);
% %Sig_diag = Sig_vec_full(Tab.Diagonal_Index());
% 
% Sig_diag = Sig_vec(diag_index_shortened);
% 
% 
% Sig_trace = sum(Sig_diag);
% if abs(Sig_trace) < 10^-16
%     warning('Sigma trace is too small1!!')
% end
% 
% %rough correction, a first try
% %first correct for complex trace
% %second compensate for absolute value
% 
% %define trace argument either by trace or by weighted mean value
% 
% argument = atan2(imag(Sig_diag), real(Sig_diag));
% 
% argument(abs(Sig_diag) < 10^-10) = 0; % to ignore rounding errors of infinitesimal blocks
% argument(abs(argument) < 10^-12) = 0; % to avoid +- eps around 0 / pi after mod function
% 
% argument = mod(argument, 2*pi);
% argument(argument > 5.5) = argument(argument > 5.5) - 2*pi;
% 
% weight = abs(Sig_diag);
% L_weight = weight > 10^-7;
% mean_arg = sum(argument(L_weight).*weight(L_weight))/sum(weight(L_weight));
% 
% 
% 
% trace_argument = atan2(imag(Sig_trace), real(Sig_trace));
% %if mean_arg and trace_argument differ a lot, there is a problem...
% if mes_flag
%     disp(['Difference in trace: ', num2str(mod(trace_argument, pi) - mod(mean_arg, pi))])
% end
% 
% Sig_vec = Sig_vec * exp(-1i*trace_argument) /abs(Sig_trace); %correction by scalar
% %denoise imaginary part:
% Sig_vec(abs(imag(Sig_vec)) < sigma_tolerance) = real(Sig_vec(abs(imag(Sig_vec)) < sigma_tolerance));
% %denoise real part:
% Sig_vec(abs(real(Sig_vec)) < sigma_tolerance) = imag(Sig_vec(abs(real(Sig_vec)) < sigma_tolerance));
% %denoise 
% Sig_vec(abs(Sig_vec) < sigma_tolerance) = 0;
% 
% 
% 
% 
% Sig_diag = Sig_vec(diag_index_shortened);
% Sig_trace = sum(Sig_diag);
% 
% 
% %check detailed:
% argument = atan2(imag(Sig_diag), real(Sig_diag));
% 
% argument(abs(Sig_diag) < 10^-12) = 0; % to ignore rounding errors of infinitesimal blocks
% argument(abs(argument) < 10^-12) = 0; % to avoid +- eps around 0 / pi after mod function
% 
% argument = mod(argument, 2*pi);
% argument(argument > 5.5) = argument(argument > 5.5) - 2*pi;
% 
% weight = abs(Sig_diag);
% mean_arg = sum(argument.*weight)/sum(weight);
% 
% variance_arg = 1/sum(weight) * sum( ( (argument - mean_arg).^2) .*weight);
% 
% % check argument of diagonal sigma (diagonal should not be 
% status.not_same_argument = 0;
% 
% L_argument = abs(Sig_diag) > 10^-10;
% 
% if variance_arg > 10^-5
%     
%     %second step
%     %L_argument_2 = abs(Sig_diag) > 10^-8;
%     %else
%     %    L_argument_2 = L_argument;
%     %end
%     %if var(argument(L_argument_2)) > 10^-13
%     if mes_flag
%         warning(['WARNING: argument of Sigma not well defined - diagonal not real, var: ', num2str(var(mod(atan2(imag(Sig_diag), real(Sig_diag)), pi)))]);
%     end
%     %TODO: Think about better way and criterion when to refine procedure, check status to identify problematic
%     %matrices
%     status.not_same_argument = 1;
%     % procedure to repair that effect for L_blocks
%     
%     %Tab_cut.f_block(Sig_vec)
%     
%     
%     arg_block = Tab_cut.f_block(sparse(diag(argument)));
%     Sig_diag_block = Tab_cut.f_block(sparse(diag(Sig_diag)));
%     L_block = cellfun(@(x) trace(abs(x)) > sigma_tolerance, Sig_diag_block);
%     
%     %cellfun(@(x) vec(x), Tab_cut.Diag_Index, 'Uniformoutput', false)
%     if sum(cellfun(@(x) var(full(diag(x))), arg_block(L_block))) > 10^-8
%         if mes_flag
%             warning(['WARNING: individual blocks do not have same complex argument!!!: ', num2str(sum(cellfun(@(x) var(full(diag(x))), arg_block(L_block))))])
%         end
%     end
%     
%     %TODO: implement a routine that does not run out of memory!
%     %TODO: check performance with large K!!!!
%     if numel(Sig_vec) < 2000
%         block = Tab_cut.f_block(Sig_vec);
%         correction = cellfun(@(x) exp(-1i*mean(diag(x))), arg_block,'Uniformoutput', false);
%         block(L_block) = cellfun(@(x,y) x*y, block(L_block), correction(L_block), 'Uniformoutput', false);
%         Sig_vec = Tab_cut.f_vector(block);
%         Sig_diag = Sig_vec(diag_index_shortened);
%         Sig_trace = sum(Sig_diag);
%     end
%     
%     Sig = Sigma([],Tab_cut.f_block(Sig_vec./Sig_trace), Tab_cut, Energy_cut, Energy_cut_2);
% else 
%     Sig = Sigma([],Tab_cut.f_block(Sig_vec./Sig_trace), Tab_cut, Energy_cut, Energy_cut_2);
% end





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
    warning('# of LEFT and RIGHT EV not equal - refine eigenvalue_tolerance')
    status.eigenvalue_tolerance = status.eigenvalue_tolerance/100;
    [Lew, Lewl, n_right, n_left, status] = func_check_tolerance(eww, ewwl, status);
else
    if n_right > 0
        status.minimal_value = min(eww(abs(eww) == min(abs(eww)))); 
    end
end
end
