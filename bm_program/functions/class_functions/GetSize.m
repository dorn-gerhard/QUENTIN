function totSize = GetSize(this) 
% get size (in Byte) of even classes / structs and variables...
props = properties(this); 
totSize = 0; 
for ii=1:length(props) 
    currentProperty = getfield(this, char(props(ii))); 
    s = whos('currentProperty'); 
    totSize = totSize + s.bytes; 
end
if isempty(props)
    s = whos('this');
    totSize = totSize + s.bytes; 
end
fprintf(1, '%d bytes\n', totSize);
