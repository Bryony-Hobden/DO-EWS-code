function [h] = H(lambda, p, C)

    if (lambda< -1) 
        h = 50/4;
        
    
    elseif (lambda>=1.5) || C >= 0.6
        h = 5;

        %varying speed according to the global mean backgorund temperature
    else  
        
        if p.delta< -4.7933
 
            h = 5;       
            
        else
            h = sample(p.delta, p);
        end   
                
    end

end
