clc
clear all
load('data6.mat')
data=data6;

%% 各积分端点B样条及其导数值
% h=0.005;
% tb=0:0.001:0.02;
% 
% for i=1:21
%     B(i)=BSPLINE( 0:0.005:0.02,tb(i),h);
%     DB(i)=DBSPLINE( 0:0.005:0.02,tb(i),h);
%     DDB(i)=DDBSPLINE( 0:0.005:0.02,tb(i),h);
% end

h=0.00125;
tb=0:2.5e-04:0.005;

for i=1:21
    B(i)=BSPLINE( 0:h:0.005,tb(i),h);
    DB(i)=DBSPLINE( 0:h:0.005,tb(i),h);
    DDB(i)=DDBSPLINE( 0:h:0.005,tb(i),h);
end

%% Romberg 求积公式
% m=0.135;
zeta=0.05; %阻尼因子
% h_span=h*4;

x=data(:,1);
h_span=data(2,1)-data(1,1);
num_use_data=(length(x)-1)/5-20;
num_it=1000;

% d0=-0.0071;
% d3=0.0242;
% 
% d1=0.0195;
% d2=0.02;

d0=-10;
d3=25;

d1=18;
d2=20;

theta_hat=zeros(20,1);
delete=0;df_delete=0;DF_delete=[];

for i=1:num_it

for j=1:num_use_data
        temp_m=0;
        temp_c=0;
        
        temp_V1_5=0;%第一段五次方系数
        temp_V1_4=0;%第一段四次方系数，以此类推
        temp_V1_3=0;
        temp_V1_2=0;
        temp_V1_1=0;
        temp_V1_0=0;
        
        temp_V2_5=0;
        temp_V2_4=0;
        temp_V2_3=0;
        temp_V2_2=0;
        temp_V2_1=0;
        temp_V2_0=0;
        
        temp_V3_5=0;
        temp_V3_4=0;
        temp_V3_3=0;
        temp_V3_2=0;
        temp_V3_1=0;
        temp_V3_0=0; 
        
        temp_F=0;
        
        t_temp=data((j-1)*5+1:(j-1)*5+1+40,1); 
        x_temp=data((j-1)*5+1:(j-1)*5+1+40,2);  %% (x(t))中的x值
        dx_temp=data((j-1)*5+1:(j-1)*5+1+40,3);  %% (dx(t))中的x值
        
        for k=1:21
            temp_m=temp_m+h_span*DDB(k)*x_temp(k);
            temp_c=temp_c+h_span*DB(k)*x_temp(k);
            
            temp_V1_5=temp_V1_5+h_span*B(k)*((x_temp(k)).^5)*phi1(x_temp(k),d0,d1,d2,d3);
            temp_V1_4=temp_V1_4+h_span*B(k)*((x_temp(k)).^4)*phi1(x_temp(k),d0,d1,d2,d3);
            temp_V1_3=temp_V1_3+h_span*B(k)*((x_temp(k)).^3)*phi1(x_temp(k),d0,d1,d2,d3);
            temp_V1_2=temp_V1_2+h_span*B(k)*((x_temp(k)).^2)*phi1(x_temp(k),d0,d1,d2,d3);          
            temp_V1_1=temp_V1_1+h_span*B(k)*((x_temp(k)).^1)*phi1(x_temp(k),d0,d1,d2,d3);
            temp_V1_0=temp_V1_0+h_span*B(k)*phi1(x_temp(k),d0,d1,d2,d3);
                     
            temp_V2_1=temp_V2_1+h_span*B(k)*((x_temp(k)).^1)*phi2(x_temp(k),d0,d1,d2,d3);
            temp_V2_0=temp_V2_0+h_span*B(k)*phi2(x_temp(k),d0,d1,d2,d3);
            
            temp_V3_3=temp_V3_3+h_span*B(k)*((x_temp(k)).^3)*phi3(x_temp(k),d0,d1,d2,d3);
            temp_V3_2=temp_V3_2+h_span*B(k)*((x_temp(k)).^2)*phi3(x_temp(k),d0,d1,d2,d3);          
            temp_V3_1=temp_V3_1+h_span*B(k)*((x_temp(k)).^1)*phi3(x_temp(k),d0,d1,d2,d3);
            temp_V3_0=temp_V3_0+h_span*B(k)*phi3(x_temp(k),d0,d1,d2,d3);
            
            %temp_F=temp_F+h_span*B(k)*5*pi^2*sin(10*pi*t_temp(k));
            %temp_F=temp_F+h_span*B(k)*0.045*pi^2*sin(30*pi*t_temp(k));
            %temp_F=temp_F+h_span*B(k)*2880*pi^2*sin(24*pi*t_temp(k));
            temp_F=temp_F+h_span*B(k)*42000*pi^2*sin(20*pi*t_temp(k));
        end
        F(j)=temp_F;
        %X theta=F
       DF(j,:)=[temp_m -temp_c temp_V1_5 temp_V1_4 temp_V1_3 temp_V1_2 temp_V1_1 temp_V1_0...
                   temp_V2_1 temp_V2_0...
                   temp_V3_3 temp_V3_2 temp_V3_1 temp_V3_0];
end

X=DF;
if delete~=0
    for d=1:delete
        index=DF_delete(d);
        DF(:,index)=0;
    end
    X=DF;
    X=X(:,~all(X==0,1));
    widthX=width(X);
end

    theta=((X'*X)^(-1)*X'*F');
    DF_del=sort(DF_delete);
    if df_delete~=0
        for d=1:delete
            index=DF_del(d);
            if index<=length(theta)
                theta=[theta(1:index-1);0;theta(index:end)];
            else
                theta=[theta(1:index-1);0];
            end
        end
    end
if length(theta_hat)==length(theta)
    sum(abs(theta_hat-theta))
    if sum(abs(theta_hat-theta))<1e-3
        fprintf('theta ok');
        break
    end
end

theta_hat=theta
  
%% d1值更新   
p1_5=theta_hat(3);
p1_4=theta_hat(4);
p1_3=theta_hat(5);
p1_2=theta_hat(6);
p1_1=theta_hat(7);
p1_0=theta_hat(8);

p2_1=theta_hat(9);
p2_0=theta_hat(10);

p3_3=theta_hat(11);
p3_2=theta_hat(12);
p3_1=theta_hat(13);
p3_0=theta_hat(14);

temp_d1=0;
temp_d2=0;
 for jjj=1:num_use_data
        t_temp=data((jjj-1)*5+1:(jjj-1)*5+1+40,1); 
        x_temp=data((jjj-1)*5+1:(jjj-1)*5+1+40,2);
        dx_temp=data((jjj-1)*5+1:(jjj-1)*5+1+40,3);
        
        for k=1:21
            temp_d1=temp_d1+h_span*B(k)*(...
                d1_phi1(x_temp(k),d0,d1,d2,d3)*...
                (p1_5*x_temp(k)^5+p1_4*x_temp(k)^4+p1_3*x_temp(k)^3+...
                p1_2*x_temp(k)^2+p1_1*x_temp(k)^1+p1_0)+...
                d1_phi2(x_temp(k),d0,d1,d2,d3)*...
                (p2_1*x_temp(k)^1+p2_0)+...     
                d1_phi3(x_temp(k),d0,d1,d2,d3)*...
                (p3_3*x_temp(k)^3+...
                p3_2*x_temp(k)^2+p3_1*x_temp(k)^1+p3_0));
            
            temp_d2=temp_d2+h_span*B(k)*(...
                d2_phi1(x_temp(k),d0,d1,d2,d3)*...
                (p1_5*x_temp(k)^5+p1_4*x_temp(k)^4+p1_3*x_temp(k)^3+...
                p1_2*x_temp(k)^2+p1_1*x_temp(k)^1+p1_0)+...
                d2_phi2(x_temp(k),d0,d1,d2,d3)*...
                (p2_1*x_temp(k)^1+p2_0)+...     
                d2_phi3(x_temp(k),d0,d1,d2,d3)*...
                (p3_3*x_temp(k)^3+...
                p3_2*x_temp(k)^2+p3_1*x_temp(k)^1+p3_0));
        end     
 end
        dd1=d1-0.000001*temp_d1;
        dd2=d2+0.0000001*temp_d2;
        
%         if abs(dd1-d1)+abs(dd2-d2)<1e-2
%             theta2=abs(theta_hat(3:end));
%             theta2(theta2==0)=NaN;
%             fprintf('delete\n');
%             df_delete=find(theta2==min(theta2))+2;
%             if abs(theta(df_delete))<0.05
%                 delete=delete+1
%                 DF_delete=[DF_delete,df_delete];%删除的项数
%             end
%         else
%             delete=delete;
%         end
        if abs(dd1-d1)+abs(dd2-d2)<1e-6
            fprintf('d1d2 ok');
            d1=dd1
            d2=dd2
            break
        end
        d1=dd1
        d2=dd2
end