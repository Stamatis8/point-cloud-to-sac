function [ srf ] = KCSsimModeler( p )
%nurbs_ship_model accepts a vector (p) of parameters
%and returns a bicubic B-Spline surface representing 
%half of a ship hull
srf = [];
X=1;
Y=2;
Z=3;
if length(p)==29
    Lwl = p(1); %Length at waterline
    B = p(2)*3*Lwl/20+Lwl/10; %Breadth: assume B between L/10 and L/4
    T = p(3)*Lwl/15+Lwl/30; %Draft: assume T between L/30 and L/10
    Mid_L = p(4) * Lwl * 0.8; %Length of midship part [0 to 80%] of Lwl
    Mid_X = p(5) * (0.8 * Lwl - Mid_L) + (0.1 * Lwl + Mid_L / 2); % position of parallel midship part's mid section
    tmp = min(B / 2, T);
    RiseOfFloor = p(6) * tmp/20;
    tmp = tmp-RiseOfFloor;
    Bilge_R = p(7) * 0.98*tmp+0.01*tmp; %Bilge radius
    FoS_Fwd_X = (1-p(8)) * (Mid_X - Mid_L / 2);
	FoS_Aft_X = Mid_X + Mid_L/2 + p(9) * (Lwl - Mid_X - Mid_L / 2);
    FoS_trans_Fwd = p(10);
    FoS_trans_Aft = p(11);
	FoB_Fwd_X = (1-p(12)) * (Mid_X - Mid_L / 2);
	FoB_Aft_X = Mid_X + Mid_L/2 + p(13) * (Lwl - Mid_X - Mid_L / 2);
    FoB_trans_Fwd = p(14);
    FoB_trans_Aft = p(15);
    BatFP = p(16)*0.45*B+0.05*B; %max is B/2 and min is 0.1B/2
	Bulb_B = p(17) * (B / 5);
    Bulb_L = p(18) * (Lwl / 25);
    Bulb_D = p(19) * (2 * T / 3);
    Bulb_tip = p(20) * Bulb_B/2;
    Fwd_rise = p(21) * (Bulb_D / 4);
    Bow_x = p(22); %Bow and Bulb intesection position
    BatAP = p(23)*0.45*B+0.05*B; %max is B/2 and min is 0.1B/2
    BatAPLow = p(24)*0.45*B+0.05*B; %max is B/2 and min is 0.1B/2
    Stern_D = p(25)*T/3;
    Stern_B = p(26)*BatAP;
    Stern_L = p(27) * (Lwl - max(FoS_Aft_X,FoB_Aft_X));
    tmp=2;
    if T/3<2
        tmp = T/3;
    end 
    Prop_H = p(28)*T/6+tmp/2;
    Prop_W = p(29)*tmp;
    %array of surface control points
    srf_pnts = zeros(3,9,12); %DIM,V=vertical,U=longitudinal
    crv_pnts = zeros(3,9);
    %create first section of midship part
    crv_pnts(:,1)=[Mid_X-Mid_L/2,B/2,0];
    crv_pnts(:,5)=crv_pnts(:,1);
    crv_pnts(Z,5)=crv_pnts(Z,1)-T+RiseOfFloor;
    crv_pnts(:,4)=crv_pnts(:,1);
    crv_pnts(Z,4)=crv_pnts(Z,5)+Bilge_R;
    crv_pnts(:,3)=crv_pnts(:,1);
    crv_pnts(Z,3)=crv_pnts(Z,4)*2/3;
    crv_pnts(:,2)=crv_pnts(:,1);
    crv_pnts(Z,2)=crv_pnts(Z,4)*1/3;
    crv_pnts(:,9)=crv_pnts(:,1);
    crv_pnts(Y,9)=crv_pnts(Y,1)-B/2;
    crv_pnts(Z,9)=crv_pnts(Z,1)-T;
    crv_pnts(:,6)=crv_pnts(:,9);
    crv_pnts(Y,6)=crv_pnts(Y,5)-Bilge_R;
    crv_pnts(Z,6)=crv_pnts(Z,9)+(RiseOfFloor*(crv_pnts(Y,6)-crv_pnts(Y,9))/(crv_pnts(Y,5)-crv_pnts(Y,9)));
    crv_pnts(:,7)=crv_pnts(:,9);
    crv_pnts(Y,7)=crv_pnts(Y,9)+(crv_pnts(Y,6)-crv_pnts(Y,9))/2;
    crv_pnts(Z,7)=crv_pnts(Z,9)+(RiseOfFloor*(crv_pnts(Y,7)-crv_pnts(Y,9))/(crv_pnts(Y,5)-crv_pnts(Y,9)));
    crv_pnts(:,8)=crv_pnts(:,9);
    crv_pnts(Y,8)=crv_pnts(Y,9)+(crv_pnts(Y,7)-crv_pnts(Y,9))/10; %for satisfying tangential continuity at hull's bottom
    MSS1 = 5; %id of first midship section
    srf_pnts(:,:,MSS1)=crv_pnts;
    %create three remaining curves by translation
    for i=1:3
        crv_pnts(X,:) = crv_pnts(X,:)+Mid_L/3;
        srf_pnts(:,:,MSS1+i)=crv_pnts;
    end
    %create the Fore FoS-FoB boundary
    crv_pnts(:,:)=srf_pnts(:,:,MSS1);
    crv_pnts(X,1)=FoS_Fwd_X;
    crv_pnts(X,4) = crv_pnts(X,4)-(crv_pnts(X,4)-FoS_Fwd_X)/20; %internal parameter
    crv_pnts(X,3) = crv_pnts(X,4)-FoS_trans_Fwd*(crv_pnts(X,4)-FoS_Fwd_X);
    crv_pnts(Y,2)=crv_pnts(Y,1);
    crv_pnts(X,2)=(3*crv_pnts(X,1)+crv_pnts(X,3))/4;
    crv_pnts(Z,2)=crv_pnts(Z,1)-abs(crv_pnts(Z,4)-crv_pnts(Z,1))/20; %internal parameter
    crv_pnts(X,5) = crv_pnts(X,5)-(crv_pnts(X,5)-max(FoB_Fwd_X,FoS_Fwd_X))/20; %internal parameter
    crv_pnts(X,6) = crv_pnts(X,6)-(crv_pnts(X,6)-FoS_Fwd_X)/20; %internal parameter
    crv_pnts(X,7) = crv_pnts(X,6)-FoB_trans_Fwd*(crv_pnts(X,6)-FoB_Fwd_X);
    crv_pnts(X,9)=FoB_Fwd_X;
    crv_pnts(X,8)=(3*crv_pnts(X,9)+crv_pnts(X,8))/4;
    srf_pnts(:,:,MSS1-1)=crv_pnts;
    %create bow profile
    crv_pnts(:,1)=[0,0,0];
    crv_pnts(X,9)=.95*FoB_Fwd_X;%internal parameter
    crv_pnts(:,8)=crv_pnts(:,9);
    crv_pnts(X,8)=crv_pnts(X,9)*0.95;
    crv_pnts(:,7)=[0,0,Fwd_rise*crv_pnts(X,9)/(crv_pnts(X,9)+Bulb_L)-T];
    crv_pnts(:,6)=[-Bulb_L,0,Fwd_rise-T];
    crv_pnts(:,5)=[-Bulb_L,0,Fwd_rise+Bulb_D+Bulb_tip-T];
    crv_pnts(:,4)=[0,0,Fwd_rise+Bulb_D-T];
    crv_pnts(:,3)=crv_pnts(:,4);
    crv_pnts(X,3)=min(FoS_Fwd_X,FoB_Fwd_X)*Bow_x;
    crv_pnts(:,2)=[crv_pnts(X,3)/2,0,crv_pnts(Z,3)/2];

    srf_pnts(:,:,1)=crv_pnts;

    % create transition curve between FoS-FoB boundary & bow profile
    crv_pnts(:,1)=[crv_pnts(X,3),BatFP,0];
    crv_pnts(X,3)=(3*crv_pnts(X,3)+srf_pnts(X,3,MSS1-1))/4;
    crv_pnts(Y,3) = BatFP/10;
    crv_pnts(:,2) = (crv_pnts(:,3)+crv_pnts(:,1))/2;
    crv_pnts(Y,4:8) = Bulb_B;
    Z1 = (crv_pnts(:,5)*3+crv_pnts(:,6))/4;
    Z2 = (crv_pnts(:,6)*3+crv_pnts(:,5))/4;
    crv_pnts(:,5) = Z1;
    crv_pnts(:,6)= Z2;
    Z1 = (crv_pnts(:,4)*3+crv_pnts(:,7))/4;
    Z2 = (crv_pnts(:,7)*3+crv_pnts(:,7))/4;
    crv_pnts(:,4) = Z1;
    crv_pnts(:,7)= Z2;
    crv_pnts(X,9)=(srf_pnts(X,9,1)+2*srf_pnts(X,9,MSS1-1))/3;
    srf_pnts(:,:,3)=crv_pnts(:,:);

    %secure bow profile tangency 
    srf_pnts(:,:,2)=srf_pnts(:,:,1);
    srf_pnts(Y,1:8,2)=(4*srf_pnts(Y,1:8,2)+srf_pnts(Y,1:8,3))/5;
    srf_pnts(X,9,2)=(2*srf_pnts(X,9,1)+srf_pnts(X,9,MSS1-1))/3;


    %create the Aft FoS-FoB boundary
    crv_pnts(:,:)=srf_pnts(:,:,MSS1+3);
    crv_pnts(X,1)=FoS_Aft_X;
    crv_pnts(Y,2)=crv_pnts(Y,1);
    crv_pnts(Z,2)=crv_pnts(Z,1)-abs(crv_pnts(Z,4)-crv_pnts(Z,1))/20;
    crv_pnts(X,9)=FoB_Aft_X;
    crv_pnts(X,4) = crv_pnts(X,4)+(FoS_Aft_X-crv_pnts(X,4))/20;
    crv_pnts(X,3) = crv_pnts(X,4)+FoS_trans_Aft*(FoS_Aft_X-crv_pnts(X,4));
    crv_pnts(X,2)=(crv_pnts(X,1)*3+crv_pnts(X,3))/4;
    crv_pnts(X,5) = crv_pnts(X,5)+(min(FoB_Aft_X,FoS_Aft_X)-crv_pnts(X,5))/20;
    crv_pnts(X,6) = crv_pnts(X,6)+(FoS_Aft_X-crv_pnts(X,6))/20;
    crv_pnts(X,7) = crv_pnts(X,6)+FoB_trans_Aft*(FoB_Aft_X-crv_pnts(X,6));
    crv_pnts(X,8)=(crv_pnts(X,9)*3+crv_pnts(X,7))/4;
    srf_pnts(:,:,MSS1+4)=crv_pnts;

    %create stern profile
    crv_pnts(:,1)=[Lwl,0,0];
    crv_pnts(:,2)=[Lwl,0,-Stern_D];
    crv_pnts(:,3)=[Lwl-Stern_L,0,-Stern_D];
    
    crv_pnts(:,5)=[Lwl-Stern_L,0,-T+Prop_H+Prop_W/2];
    crv_pnts(:,6)=[Lwl-Stern_L,0,-T+Prop_H];
    crv_pnts(:,7)=[Lwl-Stern_L,0,-T+Prop_H-Prop_W/2];
    crv_pnts(:,4)=[(Lwl-Stern_L+max(FoB_Aft_X,FoS_Aft_X))/2,0,(crv_pnts(Z,5)+crv_pnts(Z,3))/2];
    crv_pnts(:,8)=[FoB_Aft_X+0.1*Stern_L,0,-T];
    crv_pnts(:,9)=[FoB_Aft_X+0.05*Stern_L,0,-T];
    srf_pnts(:,:,end)=crv_pnts;
    
    %secure G1 continutiy for stern profile
    srf_pnts(:,:,end-1)=srf_pnts(:,:,end);
    srf_pnts(Y,1:8,end-1)=0.1*Stern_B; %internal parameter
    srf_pnts(X,9,end-1)=(2*srf_pnts(X,9,end)+FoB_Aft_X)/3;
    
    %create transition between FoB-FoS curve and stern profile
    srf_pnts(:,1,end-2) = [Lwl-Stern_L,BatAP,0];
    srf_pnts(:,2,end-2) = [Lwl-Stern_L,BatAP,-Stern_D];
    srf_pnts(:,3,end-2) = [Lwl-Stern_L,BatAP,-Stern_D];
    srf_pnts(:,4,end-2) = [((Lwl-Stern_L)*2+max(FoB_Aft_X,FoS_Aft_X))/3,BatAP/2,(srf_pnts(Z,5,end)+srf_pnts(Z,3,end))/2];
    srf_pnts(X,5:7,end-2) = srf_pnts(X,4,end-2);
    srf_pnts(Y,5:7,end-2) = BatAPLow;
    srf_pnts(Z,5,end-2) = srf_pnts(Z,5,end);
    srf_pnts(Z,6,end-2) = srf_pnts(Z,6,end);
    srf_pnts(Z,7,end-2) = srf_pnts(Z,7,end);
    srf_pnts(:,9,end-2)=[(srf_pnts(X,9,end)+2*FoB_Aft_X)/3,0,-T];
    srf_pnts(:,8,end-2)=srf_pnts(:,9,end-2);
    srf_pnts(Y,8,end-2)=0.1*Stern_B; %internal parameter
    
    srf_pnts(X,:,:)=srf_pnts(X,:,:)*-1;
    %create surface
    srf_size = size(srf_pnts);
    srf = nrbmak(srf_pnts,{UniformKnotVector(3,srf_size(2)),UniformKnotVector(3,srf_size(3))});
end
    function u = UniformKnotVector(d, nP)
        n_spans = nP-d;
        u = zeros(1,n_spans+1);
        for j=1:n_spans+1
            u(j)=(j-1)/(n_spans);
        end
        u = [zeros(1,d), u, ones(1,d)];
    end
end