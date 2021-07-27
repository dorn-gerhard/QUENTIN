function LL = filter_inp(input, category)
%elements in input should be bigger than one - check before calling this function

if nargin < 1 , error('no filter variables delivered'); end
% if input is empty, return all
if isempty(input)||  numel(input) < 1, input = {[]}; end

% vectors in category should all have same size

num_varargin = numel(input);


%input is the varargin argument containing the two possible schemes:
if ischar(input{1} )
    % ------------------------SELECTION MODE----------------------
    % 
    % user writes down which filters he sets ('SE', 1,{0,3})
    
    filter_sequence = input{1};
    is_selection_mode = 1;
else
    % ------------------------STRUCTURE MODE----------------------
    %
    % filter sequence given by structure, 
    % user just defines values of filter ([], 1, {0,3}) -> structure necessary ('NSE')
    is_selection_mode = 0;
    
    if numel(category) < num_varargin, error('Too many filter inputs for defined structure, try "NSE", 1,...'); end
    filter_sequence = cellfun(@(x) upper(x(1)), {category(1:num_varargin).Name}); % the number of filter arguments defines how much entries of structure will be used
end
    
    num_filter = numel(filter_sequence);
    if num_filter ~= num_varargin - is_selection_mode % given filter categories 'NE' corresponds to number of filters input - the one for selection mode
        error('number of filter categories does not match delivered filters');
    end
    
    %filter sequence
    % fixed categories, names should be unique, in filter sequence and in category labelling:
    % E ... Energy, include tolerance
    % S ... Spin
    % N ... Particle sector
    % D ... Double occupation / Degeneracy
    % TODO: also rethink this labelling (Categories just given by Category)
    
    %category should be table whos names identify filter categories
    FILTER_CAT_NAMES = {category(:).Name}; %full name
    FILTER_CAT_LABEL = cellfun(@(x) upper(x(1)), FILTER_CAT_NAMES); %first upper letter
    
    %FILTER = 'ESND'; %to be extended /// or match from categories itself (first letter 'upper')
    FILTER_CAT_TOL = [category(:).Tol]; % [tolerance, 0,0,0];
    
    
    filter_cat_ind = zeros(1, num_filter);
    L = cell(1,num_filter);
    LL = true(numel(category(1).Value), 1);
    
    
    for k = 1:num_filter
        [is_valid, filter_cat_ind(k)] = ismember(filter_sequence(k), FILTER_CAT_LABEL);
        if ~is_valid, error(['Filter category ', filter_sequence(k), ' does not match ', FILTER_CAT_LABEL]); end
        
        FILTER_CAT = category(filter_cat_ind(k)).Value;
        
        filter_temp = input{k+is_selection_mode};
        if iscell(filter_temp),
            if numel(filter_temp) ~= 2, error('for range filter use cell with two entries'); end
            L{k} = FILTER_CAT >= filter_temp{1} & FILTER_CAT <= filter_temp{2};
        elseif isscalar(filter_temp)
            L{k} = FILTER_CAT == filter_temp;
        elseif isempty(filter_temp),
            L{k} = true;
        elseif islogical(filter_temp) && all(size(filter_temp) == size(LL)), %if input is already a logical vector
            L{k} = filter_temp;
            
        else
            if FILTER_CAT_TOL(filter_cat_ind(k))
                %use tolerance (only works for newer Matlab > 2015a)
                vers = version;
                if str2double(vers(end-3:end-2)) >= 15
                    L{k} = ismembertol(FILTER_CAT, filter_temp, FILTER_CAT_TOL(filter_cat_ind(k)));
                else
                    if FILTER_CAT_TOL(filter_cat_ind(k))
                        filter_temp = round(filter_temp, -floor(log10(FILTER_CAT_TOL(filter_cat_ind(k)))));
                    end
                    L{k} = ismember(FILTER_CAT, filter_temp, FILTER_CAT_TOL(filter_cat_ind(k)));
                end
            else
                L{k} = ismember(FILTER_CAT, filter_temp);
            end
        end
        
        LL = LL & L{k};
    end
    
    
    
  