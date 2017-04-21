function [ Rt ] = H2Rt_Fa( G )
% Pose estiamting considering planar strucutre
% imsize [w,h]
% Reference: Faugeras 1988
% 张培科 2016
% 2016年7月31日20:51:21
% DiagEq=0.001;
% MVG P
% V3 2017年1月17日17:07:42
% The input homography must be normalized
    [~,w,~]=svd(G);
    d2=w(2,2);
    G=G/d2;
    Rt=zeros(3,4,8);
    Rp=zeros(3,3,8);
    tp=zeros(3,8);
    [U,w,V]=svd(G);
    % 2016年10月9日16:12:30
    %U=det(U)*U;
    %Vt=det(Vt)*Vt;
    s=det(U)*det(V);
    d1=w(1,1);
    d2=w(2,2);
    d3=w(3,3);
    %fprintf('Diag Matrix in SVD of A：%.5f %.5f %.5f\n',d1,d2,d3);
    %if abs(d1-d2)<=DiagEq&&abs(d2-d3)<=DiagEq%e-3
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
        Rt(:,1:3,i)=s*U*Rp(:,:,i)*V';
        %Eula(i,:)=dcm2eul(R(:,:,i))*180/pi;
        %(14)
        tp(:,i)=(d1-d3)*[x1(i);0;-x3(i)];
        Rt(:,4,i)=U*tp(:,i);
    end
    aux_sin_Phi=sqrt((d1^2-d2^2)*(d2^2-d3^2))/((d1-d3)*d2);
    cos_Phi=(d1*d3-d2^2)/((d1-d3)*d2);
    sin_Phi=[aux_sin_Phi;-aux_sin_Phi;-aux_sin_Phi;aux_sin_Phi];
    for i=5:8
        Rp(:,:,i)=[cos_Phi,   0,sin_Phi(i-4);
                    0,       -1,         0;
                    sin_Phi(i-4),0,-cos_Phi];
        Rt(:,1:3,i)=s*U*Rp(:,:,i)*V';
        %Eula(i,:)=dcm2eul(R(:,:,i))*180/pi;
        tp(:,i)=(d1+d3)*[x1(i-4);0;x3(i-4)];
        Rt(:,4,i)=U*tp(:,i);
    end    
end
