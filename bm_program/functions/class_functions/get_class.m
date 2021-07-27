function class = get_class(class_inp)
% manages creation / or loading of class depending on the input argument:
% three option:
% 1) input is object (class already)
% 2) input is a struct that has .class as field to identify which object of which class shall be created
%    the other fields of class_inp should be the parameters necessary to create object
% 3) input is a string - load object

if isobject(class_inp)
    class = class_inp;
elseif isa(class_inp, 'struct') 
    if isfield(class_inp, 'class')
        
        %param = struct2cell(class_inp);
        
        class = eval([class_inp.class, '(class_inp)']);
        %name valued Arguments
    else
        error('no "class" field in struct (class_inp) defined!')
    end
elseif isa(class_inp,'cell')
    %not implemented! - option if arguments of function have an order
    % e.g. class(var_1, var_2, var_3, ...)
    
elseif isa(class_inp, 'string')
    class = load(class_inp);
end

end