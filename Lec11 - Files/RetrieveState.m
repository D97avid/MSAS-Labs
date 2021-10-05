function [ output ] = outerfunction (x0, len, func, T , ep, e, f0, mu, tolr)

    % RHS of the outer ODE
    function fx = RHS( ~, xxx )
        r = xxx(1:end-1); % extract actual point
        rr = M(r, ep, e, f0, mu);
    
        ev=E(:,1)*((l2-l1)/(l2+l1))^2;

        % ensure directionality
        if dot((r-xlast),ev)<0
            ev = -ev;
        end
        
        fx = [ev; l2];
    end

    % save position in each step (also the initial one, ignoring flag)
    function status = callback( ~, x, flag )
        if isempty(flag) xlast = x(1:end-1,end); end
        status = 0 ;
    end
    xlast = x0;  % last integration point
    outeropt = odeset('OutputFcn',@callback, 'AbsTol', tolr, 'RelTol', tolr);
    [ ~, xr ] = ode113(@RHS,[0;len],[x0 ;0],outeropt);
    lavg = xr(end,end);
    xr = xr(1:end,1:end-1);
end