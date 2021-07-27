function [data_combined] = transform_data_combined(data_combined, form)

switch (form)
    case 'diagonal'
        names = fieldnames(data_combined);
        for k = 1:numel(names)
            
            temp_old = data_combined.(names{k});
            
            switch class(temp_old)
                case 'cell'
                    
                    siz_cell = size(temp_old);
                    if siz_cell(1) > 1 &&  siz_cell(2) > 1
                        temp_new = cell([siz_cell(1), 1, siz_cell(3:end)]);
                        
                        
                        for r = 1:siz_cell(2)
                            if length(siz_cell) > 2
                                for l = 1:prod(siz_cell(3:end))
                                    temp_new{r,1,l} = temp_old{r,r,l};
                                end
                            else
                                temp_new{r,1} = temp_old{r,r};
                            end
                        end
                    else
                        temp_new = temp_old;
                    end
                    
                    
                case 'double'
                    
                    if ~ismatrix(temp_old) % more dimensions than 2
                        siz3 = size(temp_old);
                        temp_new = zeros([siz3(1), 1, siz3(3:end)]);
                    
                        for l = 1:prod(siz3(3:end))
                            %note: this works also for multidimensional arrays > 3
                            temp_new(:,1,l) = diag(temp_old(:,:,l));
                        
                        end
                    else
                        temp_new = diag(temp_old);
                    end
                    
                case 'struct'
                    temp_new = temp_old;
                    
            end
            
            data_combined.(names{k}) = temp_new;
            
        end
            
            
            
    case 'line'
                %maybe define an additional argument to return the kth line of a matrix (horizontal or vertical)
        
    case 'mask'
        %TODO: use a mask to select relevant entries
                    
            
end


