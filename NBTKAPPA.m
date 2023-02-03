function [W]=NBTKAPPA(x,alpha)
p=3;
H=500;  
n=size(x,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inner_loop_statistic_value_array1=zeros(1,H);
inner_loop_statistic_value_array2=zeros(1,H);
inner_loop_statistic_value_array3=zeros(1,H);
inner_loop_statistic_value_array4=zeros(1,H);
inner_loop_statistic_value_array5=zeros(1,H);
inner_loop_statistic_value_array6=zeros(1,H);


   for h=1:H
   
   
     
       
       %g{i}=Random_vMF_Wood(kappa_hat_b(i),Mu=mu_hat_0{i},n=n,d=3);
       g=x(datasample((1:size(x,1)),size(x,1)),:);
       z=[sum(g(:,1)) sum(g(:,2)) sum(g(:,3))];
       q_bar1=sqrt(sum(z.^2));
       q_bar=q_bar1/size(x,1);
       q_bar1=sqrt(sum(z.^2));
      
       %mu_hat_hat_0=z./q_bar1;
       
       
        kappa_hat_hat_b=((q_bar*(p-q_bar^2))/(1-q_bar^2));
        
        
        p1=kappa_hat_hat_b;
        b0=bes(p,p1);
        

         b1=(-(p-1)*(b0/p1)+1-b0^2);
 
         b2=(p-1)*(b0/(p1)^2)-(p-1)*(b1/p1)-2*b0*b1;

         b3= -2*(p-1)*(b0/(p1)^3)+2*(p-1)*(b1/(p1)^2)-(p-1)*(b2/p1)-2*(b1^2)-2*b0*b2;
 
         b4=6*(p-1)*(b0/(p1)^4)-6*(p-1)*(b1/(p1)^3)+3*(p-1)*(b2/(p1)^2)-(p-1)*(b3/p1)-6*b1*b2-2*b0*b3;
 
         b5= -24*(p-1)*(b0/(p1)^5)+24*(p-1)*(b1/(p1)^4)-12*(p-1)*(b2/(p1)^3)+4*(p-1)*(b3/(p1)^2)-(p-1)*(b4/p1)-6*(b2^2)-8*b1*b3-2*b0*b4;
        
        
        
        delta2=((p+1)/(n*(p-1)));
        delta15_1=(1-(p+1)/(n*(p-1)));
        delta15_2=(p*(p-3))/(2*n*(p-1));
        delta16_1=1-(3/(n*(p+2)));
        delta16_2=(p*(p-1))/(2*n*kappa_hat_hat_b);
        delta17=((p-1)*b1-kappa_hat_hat_b*b2)/(2*n*kappa_hat_hat_b*(b1^2));
        delta18_1=1-(3*(p+1))/(2*n*(p+2));
        delta18_2=(p*(p-1))/(2*n*kappa_hat_hat_b);
%       
        l(1)=(kappa_hat_hat_b*(1)-0);
        l(2)=(kappa_hat_hat_b*(1-delta2)-0);
        l(3)=(kappa_hat_hat_b*(delta15_1)-delta15_2);
        l(4)=(kappa_hat_hat_b*(delta16_1)-delta16_2);
        l(5)=(kappa_hat_hat_b*(1)-delta17);
        l(6)=(kappa_hat_hat_b*(delta18_1)-delta18_2);
    

  inner_loop_statistic_value_array1(h)= l(1);
  inner_loop_statistic_value_array2(h)=l(2);
  inner_loop_statistic_value_array3(h)=l(3);
  inner_loop_statistic_value_array4(h)=l(4);
  inner_loop_statistic_value_array5(h)=l(5);
  inner_loop_statistic_value_array6(h)=l(6);
   end
  
lower1= quantile(inner_loop_statistic_value_array1,  alpha/2);
lower2= quantile(inner_loop_statistic_value_array2,  alpha/2);
lower3= quantile(inner_loop_statistic_value_array3,  alpha/2);
lower4= quantile(inner_loop_statistic_value_array4,  alpha/2);
lower5= quantile(inner_loop_statistic_value_array5,  alpha/2);
lower6= quantile(inner_loop_statistic_value_array6,  alpha/2);




    upper1= quantile(inner_loop_statistic_value_array1, 1 - alpha/2);
    upper2= quantile(inner_loop_statistic_value_array2, 1 - alpha/2);
    upper3= quantile(inner_loop_statistic_value_array3, 1 - alpha/2);
    upper4= quantile(inner_loop_statistic_value_array4, 1 - alpha/2);
    upper5= quantile(inner_loop_statistic_value_array5, 1 - alpha/2);
    upper6= quantile(inner_loop_statistic_value_array6, 1 - alpha/2);
    
    

    l1=upper1-lower1;
    l2=upper2-lower2;
    l3=upper3-lower3;
    l4=upper4-lower4;
    l5=upper5-lower5;
    l6=upper6-lower6;
    
  
W=[ l1 l2  l3 l4  l5 l6];
end