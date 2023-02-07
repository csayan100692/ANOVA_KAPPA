function [W]=PBTKAPPA1(x,alpha)
p=3;
%B=200;
H=500;
n=size(x,1);
kappa_0=2;
%kappa_alt=1.5;
gamma_0=(5*kappa_0^3)/(5*kappa_0^2-1);
%mu=[1 1 1]/sqrt(3);
%statistic_value_array=zeros(1,B);
%statistic_value_array2=zeros(1,B);
%critical_value_array=zeros(1,B);
inner_loop_statistic_value_array=zeros(1,H);
%alpha=0.05;
%for b=1:B
  
    
   % x=Random_vMF_Wood(kappa_0,Mu=mu,n=n,d=3);
    
    
    y=[sum(x(:,1)) sum(x(:,2)) sum(x(:,3))];
    
    r_bar=sqrt(sum(y.^2))/n;
    r_bar1=sqrt(sum(y.^2));
    
    mu_hat_0=y./r_bar1;
    
    %k_hat_0{i}=matrix(c(mu_hat_0[[i]]),nrow = p-1,ncol = 1)
    %l=(r_bar*(p-r_bar^2))/(1-r_bar^2);
    %v1=1/(n*bes1(p,l));
    %statistic_value_array(b)=(l-kappa_0)/sqrt(v1);
    statistic_value=2*n*kappa_0*(1-r_bar);
    %statistic_value_array2(b)=2*n*gamma_0*(1-r_bar);
  
    
   for h = 1:H
   
     g=Random_vMF_Wood(kappa_0,Mu=mu_hat_0,n=n,d=3);
     %g1 <- x1[sample(nrow(x1),replace = TRUE),]
     z1=[sum(g(:,1)) sum(g(:,2)) sum(g(:,3))];
     q_bar=sqrt(sum(z1.^2))/n;
     q_bar1=sqrt(sum(z1.^2));
     
     %mu_hat1=z1/q_bar1;
     
     
     %l2=(q_bar*(p-q_bar^2))/(1-q_bar^2);
     %v2=1/(n*bes1(p,l2));
     %inner_loop_statistic_value_array(h)=(l2-kappa_0)/sqrt(v2);
     %inner_loop_statistic_value_array(h)=2*n*gamma_0*(1-q_bar);
     inner_loop_statistic_value_array(h)=2*n*kappa_0*(1-q_bar);
   end
   critical_value=quantile(inner_loop_statistic_value_array, 1-alpha);
   %end
%  for b=1:B
%   
%     
%     x=Random_vMF_Wood(kappa_alt,Mu=mu,n=n,d=3);
%     
%     
%     y=[sum(x(:,1)) sum(x(:,2)) sum(x(:,3))];
%     
%     r_bar=sqrt(sum(y.^2))/n;
%     r_bar1=sqrt(sum(y.^2));
%     
%     mu_hat_0=y./r_bar1;
%     
%     %k_hat_0{i}=matrix(c(mu_hat_0[[i]]),nrow = p-1,ncol = 1)
%     l=(r_bar*(p-r_bar^2))/(1-r_bar^2);
%     v1=1/(n*bes1(p,l));
%     %statistic_value_array2(b)=(l-kappa_0)/sqrt(v1);
%     statistic_value_array2(b)=2*n*gamma_0*(1-r_bar);
%  end
% count=0;
% for i = 1:B
% 
%   if statistic_value_array2(i)>critical_value_array(i)
%   %if statistic_value_array2(i)>chi2inv(1-alpha,(n-1)*(p-1))
%   %if statistic_value_array(i)>norminv(1-alpha)
%     count=count+1;
%   end
% end

W=[statistic_value, critical_value];
end