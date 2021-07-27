function struct=vars2struct(varargin)


for counter=1:numel(varargin)
    thisvar=evalin('caller', varargin{counter});
    THEWORKSPACE.(varargin{counter})=thisvar;
end

struct=THEWORKSPACE;

end
