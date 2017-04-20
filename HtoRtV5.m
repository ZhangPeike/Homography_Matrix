function [ Rt,isPureRotation ] = HtoRtV5( H,K1,K2,matches )
% Pose estiamting considering planar strucutre
% Reference: Faugeras 1988
% Zhang Peike 2016
% 2016年7月31日20:51:21
% MVG P
% V3 2017年1月17日17:07:42
% V5 2017年3月2日08:02:23
% Approximate R correctly.
% Modified the solution set,
% Treat w2 as to be 1.
A=K2\H*K1;
%A=H;
[~,w,~]=svd(A);
d2=w(2,2);
fprintf('Diag Matrix in SVD of A:%.3f %.3f %.3f\n',w(1,1),w(2,2),w(3,3));
%Enforce the property 1
disp(A);
A=A/d2;
A=A*det(A);
disp(A);
[U,w,V]=svd(A);
% 2016年10月9日16:12:30
%U=det(U)*U;
%Vt=det(Vt)*Vt;
% s=det(U)*det(V);
% fprintf('Variable s sign:%.4f\n',s);
% d2 must be equal to 1
%d2=w(2,2);
fprintf('Diag Matrix in SVD of G:%.5f %.5f %.5f\n',w(1,1),w(2,2),w(3,3));
d1=w(1,1);
d2=w(2,2);
d3=w(3,3);
if abs(d1/d3)<=1.02
    fprintf('Pure Rotation is found!\n');
   
    isPureRotation=true;
    Rt=[U*V,zeros(3,1)];
    fprintf('Pure Rotation Check (det sign):%.1f\n',det(Rt(:,1:3)));
    return;
else
    isPureRotation=false;
end
Rp=zeros(3,3,4);
tp=zeros(3,4);
R=zeros(3,3,4);
t=zeros(3,4);
% if s>0
%     disp('s>0');
aux1=sqrt((d1^2-d2^2)/(d1^2-d3^2));
aux3=sqrt((d2^2-d3^2)/(d1^2-d3^2));
x1=[aux1;aux1;-aux1;-aux1];
x3=[aux3;-aux3;aux3;-aux3];
aux_sin_Theta=sqrt((d1^2-d2^2)*(d2^2-d3^2))/((d1+d3)*d2);
cos_Theta=(d2^2+d1*d3)/((d1+d3)*d2);
sin_Theta=[aux_sin_Theta;-aux_sin_Theta;-aux_sin_Theta;aux_sin_Theta];
for i=1:4
    %(13)
    Rp(:,:,i)=[cos_Theta,0,-sin_Theta(i);
        0,        1,         0;
        sin_Theta(i),0,cos_Theta];
    % Faugeras (8)
    R(:,:,i)=s*U*Rp(:,:,i)*V';
    %Eula(i,:)=dcm2eul(R(:,:,i))*180/pi;
    %(14)
    tp(:,i)=(d1-d3)*[x1(i);0;-x3(i)];
    t(:,i)=U*tp(:,i);
end
%end
%{
if s<0
    disp('s<0');
    aux_sin_Phi=sqrt((d1^2-d2^2)*(d2^2-d3^2))/((d1-d3)*d2);
    cos_Phi=(d1*d3-d2^2)/((d1-d3)*d2);
    sin_Phi=[aux_sin_Phi;-aux_sin_Phi;-aux_sin_Phi;aux_sin_Phi];
    for i=1:4
        Rp(:,:,i)=[cos_Phi,   0,sin_Phi(i);
                    0,       -1,         0;
                    sin_Phi(i),0,-cos_Phi];
        R(:,:,i)=s*U*Rp(:,:,i)*V';
        %Eula(i,:)=dcm2eul(R(:,:,i))*180/pi;
        tp(:,i)=(d1+d3)*[x1(i);0;x3(i)];
        t(:,i)=U*tp(:,i);
    end
end
%}
P1=K1*[eye(3) zeros(3,1)];
Vote=zeros(1,4);
for i=1:4
    P2=K2*[R(:,:,i),t(:,i)];
    Vote(i)=0;
    % Triangulate
    X=Triangulation( P1,P2,matches );
    % Turn into inhomogenerous
    X = X(1:3,:) ./ X([4 4 4],:);
    mmcos=zeros(1,size(X,2));
    for m=1:size(X,2)
        %齐次坐标转化为非齐次
        %X(:,m)=X(1:3,m)./X([4 4 4],m);
        %判断从第二个相机焦点到X点的矢量(X(:,m)-Rt(:,4,i))与第二个相机的光轴的点积
        %R矩阵的第三行Rt(3,1:3,i)即Z轴方向的单位矢量在世界坐标系的坐标
        %四个解的判断方法与SFMedu的比较理解
        %2016年5月25日17:00:04，更正，第二个相机的角点坐标C=-R'*t
        %可增加重构点与两个光轴的夹角
        mmcos(m)=R(3,:,i)*(X(:,m)+R(:,:,i)'*t(:,i));        
        Vote(i)=(X(3,m)>0&&mmcos(m)>0)+Vote(i);
    end    
end
%disp(dcm2eul(R{1})*180/pi);
%disp(dcm2eul(R{2})*180/pi);
[~,index]=max(Vote);
R=R(:,:,index);
t=t(:,index);
Rt=[R t];

disp('HtoRtV5 4 Candidate Poses');
disp(Vote);
fprintf('Unique R:%.4f %.4f %.4f\n',dcm2eul(R)*180/pi);
fprintf('Unique t:%.4f %.4f %.4f\n',t);

end
