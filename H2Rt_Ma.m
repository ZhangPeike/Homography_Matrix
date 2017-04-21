function [ Rt ] = H2Rt_Ma( G )
%Reference MVG&MaYi
%20160809
%2017-04-02 10:00
%ZhangPeike
%V2
%ZhangPeike finding
%V3 only putout the candidate pose
%The input homography must be normalized
    Rt=zeros(3,4,4);
    if det(G)<0
        G=-G;
    end
    [~,w,~]=svd(G);
    d2=w(2,2);
    %Ma Yi P135 Lemma 5.18
    G=G/d2;
    %Reference
    AtA=G'*G;
    %[V,D]=eig(AtA);
    % sigmaSqua Ascending order
    %[V,sigmaSqua]=schur(AtA);
    [~,sigmaSqua,V]=svd(AtA);
    %Ma book P136 V belong to SO3
    V=det(V)*V;
    u1=(sqrt(1-sigmaSqua(3,3))*V(:,1)+sqrt(sigmaSqua(1,1)-1)*V(:,3))/sqrt(sigmaSqua(1,1)-sigmaSqua(3,3));
    u2=(sqrt(1-sigmaSqua(3,3))*V(:,1)-sqrt(sigmaSqua(1,1)-1)*V(:,3))/sqrt(sigmaSqua(1,1)-sigmaSqua(3,3));
    U1=[V(:,2),u1,vec2skew_symme(V(:,2))*u1];
    W1=[G*V(:,2),G*u1,vec2skew_symme(G*V(:,2))*G*u1];
    U2=[V(:,2),u2,vec2skew_symme(V(:,2))*u2];
    W2=[G*V(:,2),G*u2,vec2skew_symme(G*V(:,2))*G*u2];
    %MaYi P137 table 5.1
    %t is the ratio between t and d the plane to the second camera center distance.
    Rt(:,:,1)=[W1*U1',(G-W1*U1')*U1(:,3)];
    Rt(:,:,2)=[W2*U2',(G-W2*U2')*U2(:,3)];
    Rt(:,:,3)=[Rt(:,1:3,1),-Rt(:,4,1)];
    Rt(:,:,4)=[Rt(:,1:3,2),-Rt(:,4,2)];
end