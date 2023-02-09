function [W]=PBTKAPPA3(x,alpha,kappa_0)
p=3;
H=500;
n=size(x,1);
inner_loop_statistic_value_array=zeros(1,H);
gamma_0=(5*kappa_0^3)/(5*kappa_0^2-1);
    y=[sum(x(:,1)) sum(x(:,2)) sum(x(:,3))];
    
    r_bar=sqrt(sum(y.^2))/n;
    r_bar1=sqrt(sum(y.^2));
    
    mu_hat_0=y./r_bar1;


    l=(r_bar*(p-r_bar^2))/(1-r_bar^2);
    v1=1/(n*bes1(p,l));
    
    statistic_value=(l-kappa_0)/sqrt(v1);
    
   for h = 1:H
   
     g=Random_vMF_Wood(kappa_0,Mu=mu_hat_0,n=n,d=3);
     %g1 <- x1[sample(nrow(x1),replace = TRUE),]
     z1=[sum(g(:,1)) sum(g(:,2)) sum(g(:,3))];
     q_bar=sqrt(sum(z1.^2))/n;
     q_bar1=sqrt(sum(z1.^2));
     l2=(q_bar*(p-q_bar^2))/(1-q_bar^2);
     v2=1/(n*bes1(p,l2));
     inner_loop_statistic_value_array(h)=(l2-kappa_0)/sqrt(v2);
     
     %inner_loop_statistic_value_array(h)=2*n*gamma_0*(1-q_bar);
     %inner_loop_statistic_value_array(h)=2*n*kappa_0*(1-q_bar);
   end
   critical_value=quantile(inner_loop_statistic_value_array, 1-alpha);
W=[statistic_value, critical_value];
end