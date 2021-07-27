function str = plot_parameters(parameters, str, recursionlevel)
output = false; firstlevelmarker = '';
if nargin < 2 || isempty(str), output = true; str = []; firstlevelmarker = '\n';  end
if nargin < 3 || isempty(str), recursionlevel = 0; end

names = fieldnames(parameters);
for k = 1:numel(fieldnames(parameters))
    object = getfield(parameters, names{k});
    objectname = names{k};
    
    
    if isstruct(object)
        %recursion
        str = [str, firstlevelmarker, repmat('     ', 1, recursionlevel), names{k}, ':\n'];
        str = plot_parameters(object, str, recursionlevel + 1);
    elseif isobject(object)
        %disp('there is an object')
        str = [str, firstlevelmarker, repmat('     ', 1, recursionlevel), names{k}, ':\n'];
        str = plot_parameters(object, str, recursionlevel + 1);
    else
        if numel(object) <= 1
            
            
            if iscell(object)
                %warning('asdf')
                if isempty(object)
                    object = [];
                else
                    object = object{1};
                end
                if isobject(object)
                    %warning('class inside')
                    str = [str, firstlevelmarker, repmat('     ', 1, recursionlevel), names{k}, ':\n'];
                    str = plot_parameters(object, str, recursionlevel + 1);
                else
                    str = [str, firstlevelmarker, repmat('     ', 1, recursionlevel), objectname, ': ', repmat(' ', 1,20-numel(objectname)), num2str(object), '\n'];
                end
            else
                
                str = [str, firstlevelmarker, repmat('     ', 1, recursionlevel), objectname, ': ', repmat(' ', 1,20-numel(objectname)), num2str(object), '\n'];
            end
            
            
        else
            if ischar(object), min_length  = 30; else min_length = 8; end
            object_to_print = object(1:min(end,min_length));
            if iscell(object_to_print)
                if isobject(object_to_print{1})
                    disp('object')
                    str = [str, firstlevelmarker, repmat('     ', 1, recursionlevel), names{k}, ':\n'];
                    str = plot_parameters(object, str, recursionlevel + 1);
                    
                else
                    str = [str, firstlevelmarker, repmat('     ', 1, recursionlevel), objectname, ': ', repmat(' ', 1,20-numel(objectname)), 'cell' ,...
                        '  [', num2str(size(object,1)), ']x[', num2str(size(object,2)), '] = ', num2str(numel(object)) ' elements\n'];
                end
            else
                
                 str = [str, firstlevelmarker, repmat('     ', 1, recursionlevel), objectname, ': ', repmat(' ', 1,20-numel(objectname)), num2str(object_to_print(:).', 2),...
                     '  [', num2str(size(object,1)), ']x[', num2str(size(object,2)), '] = ', num2str(numel(object)) ' elements\n'];
            end
        end
    end
end



if output == true
    sprintf(str)
end
if recursionlevel == 0
    str = replace(str,'\n','\newline');
    str = replace(str,'_', '\_');
end


end