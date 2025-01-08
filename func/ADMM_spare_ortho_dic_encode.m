function [w,z] = ADMM_spare_ortho_dic_encode(x,phi,lambda,rho)


        [p,t]=size(phi);

        num_iters = 500;
        tol = 10^-6;
        
       %initialization
       I = eye(p);
       w = rand(1,p);
       z = w;
       Gamma=0;
       
       fit = 0;
       for i=1:num_iters
           
           
            % update w equation 12
             w = (x*phi'+rho*z+Gamma)./(1+rho);

            % update z equation 14 and 15 
             h =w-Gamma/rho;
             z = sign(h).*max(abs(h)-lambda/rho,0);

             % update Gamma equation 16
             Gamma = Gamma + rho*(z-w);

             old_fit = fit;
             fit = norm(x-w*phi)+ sum(abs(z),'all');
             
             fit_change =fit-old_fit;

%             disp(['obj-',num2str(i),'=',num2str(fit ),',','change=',num2str(fit_change)]);

             % if our fit hasn't changed much return
             if abs(fit_change)< tol
                
                 return

             end


            
       end



end

