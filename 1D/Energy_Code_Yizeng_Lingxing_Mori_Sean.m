%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Two phase Cell Migration (1D Code)
%  Contact: Yizeng Li (liyizeng52@hotmail.com or yli54@kennesaw.edu)
%           Lingxing Yao (lyao@uakron.edu)
%           Yoichiro Mori (y1mori@math.upenn.edu)
%           Sean Sun (ssun@jhu.edu)
%  Updated on: October 11, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Steady-state solution
% Given: theta_c/n, Jactin, pressure in the network
% Unknowns: Pc, vc, vn, thetan, v0
% In this version, the velocity is in nm/s to reduce numerical error,
% the corresponding coefficients have been rescaled with nm/s.
% Use a matrix index system to plot all the contours
% Use a nonlinear focal adhesion vs actin velocity relation.
% If we let c2 be infinite, then we go back to
% the original approximation that the focal adhesion friction is a linear
% function in actin velocity.
% Include electro-netural ions
% Focal traction: -eta_st*theta_n*v_n/(1+(vn/c2)^m)

% 1: Jactin
% 2: nust
% 3: xim
% 4: c0f
% 5: dg
% 6: Phin
% 7: Thetan
% 8: eta
% 9: etan
% 10: zetac
% 11: zetan
% 12: fextf
% 13: fextb
% 14: Jwater (from Jcactivef)
% 15: nust0
% 17: c0b
% 18: Jcactivef
% 19: p0b
% 21: alpha
% 22: sigmaa0

clear
clc

RankCheck = 0;

Parameter1 = 1;
Parameter2 = 2;

Iter = 21;
b = 3.d-6;          % (m) cell width
w = 10.d-6;         % (m) cell depth
L = 50.d-6;         % (m) cell length
h = 0.5d-6;         % (m) membrane thickness
p0f = 0*1d5;            % (Pa) external pressure at the front
p0b = 0*1d5;            % (Pa) external pressure at the back
c0f = 340;            % (mol/m^3 = mM) external ion concentration at the front
c0b = 340;            % (mol/m^3 = mM) external ion concentration at the back
cinb = 340.35;        % expected value, not prescribed
cinf = 340.87;        % expected value, not prescribed
Dcommon = 1.d-3;    % (nm m/s) diffusion constant for ions
vc0 = 0;           % (nm) expected average cytosol velocity

R = 8.31451;        % (J/mol K) Ideal gas constant
T = 300;            % (K) absolute temperature
kB = 1.38d-23;      % (J/K) Boltzmann constant
NA = 6.02d23;       % (1/mol) Avogadios number

AF = 25.d-18;       % (m^2) cross sectional area of one actin filament
delta = 3;          % (nm) effective length of each G-actin
GATP = 25*kB*T;     % (J) energy of each ATP
NC = 1;             % each pumped ion takes this number of ATP

alphaf = 1.d-1;     % (nm/Pa/s) coefficient of water permeation
alphab = 1.d-1;     % (nm/Pa/s) coefficient of water permeation
gf = 5d4;         % (nm/s) passive channel coefficient at the front
gb = 5d4;         % (nm/s) passive channel coefficient at the back
% Jcactivef = 1.d-6;      % (mM nm/s) actin flux at the front
% Jcactiveb = 1*1d-2*(gf*(cinf-c0f)+Dcommon/L*(cinf-cinb)-vc0*(cinb+cinf)/2);%gb*(cinb-c0b)-Dcommon/L*(cinf-cinb)-vc0*(cinb+cinf)/2;      % (mM m/s) actin flux at the back
% Jcactivef = 1*1d2*Jcactiveb;
% Jcactiveb = Jcactivef;
Jcactiveb = 0*(gb*(cinb-c0b)-Dcommon/L*(cinf-cinb)+vc0*cinb);
Jcactivef = 0*(gf*(cinf-c0f)+Dcommon/L*(cinf-cinb)-vc0*cinf);
% Jcactiveb = 0.01*(gf*(cinf-c0f)+Dcommon/L*(cinf-cinb)-vc0*(cinb+cinf)/2)*1.4368;


fextf = 0d2;  % (Pa) external force per unit area at the front of the cell
fextb = 0d2;  % (Pa) external force per unit area at the back of the cell

eta  = 1.d3;         % (Pa s/m nm) drag coefficient between two phases
etan = 12/9*1.d3;         % (Pa s/m nm) effective viscous coefficient for the actin phase
xim  = 1.d-3;         % (Pa s/nm) drag coefficient between the membrane and the substrate
dg   = 0*1d3;         % (Pa s/nm) drag coefficient from the displaced water

Thetan = 0.02;
theta_0 = Thetan/sqrt(2);
theta_star = Thetan;
Jactinf = 0*50*Thetan;    % (nm/s) actin flux at the front the cell

N = 21;
dx = L/(N-1);
x = linspace(0,L,N);

Nucleus = round(N/3):round(2*N/3);
Lzeta = (length(Nucleus)-1)*dx;
rho0 = 3.4;             % (mM) concentration of the actin network
sigmap0 = R*T*rho0;     % swelling pressure the actin network
sigmaa0 = 0;         % active contractile pressure from myosin
sigmaa = -sigmaa0*linspace(1,0,N)';     
nust = 1.d7*ones(N,1);     % (Pa s/nm m) drag coefficient between the actin and the substrate
c2 = 10*ones(N,1)*1d3;      % (nm/s)
m = 4;          % exponent in focal traction
nust0 = 1.d5;
Zetac = 0.d4;
Zetan = 0.d5;
zetac = zeros(N,1);
zetan = zeros(N,1);
D = Dcommon*ones(N,1);   % (nm m/s) diffusion constant for ions
zetac(Nucleus) = Zetac;
zetan(Nucleus) = Zetan;
Dzeta = Dcommon;         % (nm m/s) diffusion constant for ions at the nucleus
D(Nucleus) = Dzeta;


if Parameter1 == 0 
    N1 = 1;
else
    N1 = 31;
end
if Parameter2 == 0
    N2 = 1;
else
    N2 = 33;
end

if Parameter1 == 1
    JACTINF = linspace(1,50,N1);
    Ym = JACTINF;
elseif Parameter1 == 2
    NUST = logspace(3,7,N1);
    Ym = NUST*1d9/1d12;
elseif Parameter1 == 3
    XIM = logspace(-2,0,N1);
    Ym = XIM*1d9/1d6;
elseif Parameter1 == 4
    C0F = linspace(200,340,N1);
    Ym = C0F;
elseif Parameter1 == 5
    DG = logspace(-1,3,N1);
    Ym = DG*1d9/1d6;
elseif Parameter1 == 6
    PHIN = logspace(1,4,N1);
    Ym = PHIN;
elseif Parameter1 == 7
    THETAN = linspace(0.01,0.05,N1);
    Ym = THETAN;
elseif Parameter1 == 8
    ETA = logspace(3,7,N1);
    Ym = ETA*1d9/1d12;
elseif Parameter1 == 9
    ETAN = logspace(1,7,N1);
    Ym = ETAN*1d9;
elseif Parameter1 == 10
    ZETAC = logspace(1,6,N1);
    Ym = ZETAC*1d9;
elseif Parameter1 == 11
    ZETAN = logspace(2,7,N1);
    Ym = ZETAN*1d9;
elseif Parameter1 == 12
    FEXTF = linspace(-1d2,1d2,N1);
    Ym = FEXTF;
elseif Parameter1 == 13
    FEXTB = linspace(-1d2,1d2,N1);
    Ym = FEXTB;
elseif Parameter1 == 17
    C0B = linspace(0.6,1,N1);
    Ym = C0B;
elseif Parameter1 == 18 || Parameter1 == 14
    JACTIVEF = (gf*(cinf-c0f)+Dcommon/L*(cinf-cinb)-vc0*(cinb+cinf)/2)*linspace(0.01,1,N1)*1.4368;
    Ym = JACTIVEF;
elseif Parameter1 == 19 
    P0B = p0b + linspace(0,500,N1);
    Ym = P0B;
elseif Parameter1 == 21 
    ALPHA = logspace(-3,1,N1);
    Ym = ALPHA*1d9/1d6;
elseif Parameter1 == 22 
    SIGMAA0 = logspace(0,3,N1);
    Ym = SIGMAA0;
end

if Parameter2 == 1
    JACTINF = linspace(0,50,N2);
    Xm = JACTINF;
elseif Parameter2 == 2
    NUST = logspace(3,7,N2);
    Xm = NUST*1d9/1d12;
elseif Parameter2 == 3
    XIM = logspace(-4,-1,N2);
    Xm = XIM*1d9/1d6;
elseif Parameter2 == 4
    C0F = linspace(200,340,N2);
    Xm = C0F;
elseif Parameter2 == 5
    DG = logspace(-1,3,N2);
    Xm = DG*1d9/1d6;
elseif Parameter2 == 6
    PHIN = logspace(1,4,N2);
    Xm = PHIN;
elseif Parameter2 == 7
    THETAN = linspace(0.01,0.05,N2);
    Xm = THETAN;
elseif Parameter2 == 8
    ETA = logspace(1,7,N2);
    Xm = ETA*1d9/1d12;
elseif Parameter2 == 9
    ETAN = logspace(1,7,N2);
    Xm = ETAN*1d9;
elseif Parameter2 == 10
    ZETAC = logspace(1,6,N2);
    Xm = ZETAC*1d9;
elseif Parameter2 == 11
    ZETAN = logspace(2,7,N2);
    Xm = ZETAN*1d9;
elseif Parameter2 == 12
    FEXTF = linspace(-2d3,0,N2);
    Xm = FEXTF;
elseif Parameter2 == 13
    FEXTB = linspace(0,200,N2);
    Xm = FEXTB;
elseif Parameter2 == 17
    C0B = linspace(0.5,1,N2);
    Xm = C0B;
elseif Parameter2 == 18 || Parameter2 == 14
    JACTIVEF = (gf*(cinf-c0f)+Dcommon/L*(cinf-cinb)-vc0*(cinb+cinf)/2)*linspace(0.01,1,N2);
    Xm = JACTIVEF;
elseif Parameter2 == 19
    P0B = p0b + linspace(0,500,N2);
    Xm = P0B;
elseif Parameter2 == 21 
    ALPHA = logspace(-3,1,N2);
    Xm = ALPHA*1d9/1d6;
elseif Parameter2 == 22 
    SIGMAA0 = logspace(0,3,N2);
    Xm = SIGMAA0;
end

VN = zeros(N1,N2);
VC = zeros(N1,N2);
V0 = zeros(N1,N2);
JWATERF = zeros(N1,N2);
TAUF = zeros(N1,N2);
TAUB = zeros(N1,N2);
PSTARF = zeros(N1,N2); 
PSTARB = zeros(N1,N2); 
PCF = zeros(N1,N2);
PCB = zeros(N1,N2);
V0_L = zeros(N1,N2);
KC = zeros(N1,N2);


ICELL = zeros(N1,N2);
IEXT  = zeros(N1,N2);
DFRIC = zeros(N1,N2);
DFLOW = zeros(N1,N2);
DSOL  = zeros(N1,N2);
E_SUBSTRATE = zeros(N1,N2);
E_MEMBRANE  = zeros(N1,N2);
E_HYDRAULIC = zeros(N1,N2);
E_EXTERNAL  = zeros(N1,N2);
E_PUMP      = zeros(N1,N2);
E_ACTIN     = zeros(N1,N2);
E_MYOSIN    = zeros(N1,N2);
E_INTERFACE = zeros(N1,N2);
E_NUCLEUS   = zeros(N1,N2);
E_FLOW      = zeros(N1,N2);
ATP_ACTIN   = zeros(N1,N2);
ATP_ION     = zeros(N1,N2);


for loop2 = 1:N2
    
    loop2
    
    if Parameter2 == 1
        Jactinf = JACTINF(loop2)*Thetan;
    elseif Parameter2 == 2
        nust = NUST(loop2)*ones(N,1);
    elseif Parameter2 == 3
        xim = XIM(loop2);
    elseif Parameter2 == 4
        c0f = C0F(loop2);
    elseif Parameter2 == 5
        dg = DG(loop2);
    elseif Parameter2 == 6
        sigmap0 = PHIN(loop2);
    elseif Parameter2 == 7
        Thetan = THETAN(loop2);
    elseif Parameter2 == 8
        eta = ETA(loop2);
    elseif Parameter2 == 9
        etan = ETAN(loop2);
    elseif Parameter2 == 10
        zetac(Nucleus) = ZETAC(loop2);
    elseif Parameter2 == 11
        zetan(Nucleus) = ZETAN(loop2);
    elseif Parameter2 == 12
        fextf = FEXTF(loop2);
    elseif Parameter2 == 13
        fextb = FEXTB(loop2);
    elseif Parameter2 == 17
        c0b = C0B(loop2);
    elseif Parameter2 == 18 || Parameter2 == 14
        Jcactivef = JACTIVEF(loop2);
    elseif Parameter2 == 19
        p0b = P0B(loop2);
    elseif Parameter2 == 21
        alphaf = ALPHA(loop2);
        alphab = ALPHA(loop2);
    elseif Parameter2 == 22
        sigmaa0 = SIGMAA0(loop2);
   end
    
    for loop1 = 1:N1
        
        loop1; 
        
        if Parameter1 == 1
            Jactinf = JACTINF(loop1)*Thetan;
        elseif Parameter1 == 2
            nust = NUST(loop1)*ones(N,1);
        elseif Parameter1 == 3
            xim = XIM(loop1);
        elseif Parameter1 == 4
            c0f = C0F(loop1);
        elseif Parameter1 == 5
            dg = DG(loop1);
        elseif Parameter1 == 6
            sigmap0 = PHIN(loop1)*ones(N,1).*linspace(1,5,N)';
        elseif Parameter1 == 7
            Thetan = THETAN(loop1);
        elseif Parameter1 == 8
            eta = ETA(loop1);
        elseif Parameter1 == 9
            etan = ETAN(loop1);
        elseif Parameter1 == 10
            zetac(Nucleus) = ZETAC(loop1);
        elseif Parameter1 == 11
            zetan(Nucleus) = ZETAN(loop1);
        elseif Parameter1 == 12
            fextf = FEXTF(loop1);
        elseif Parameter1 == 13
            fextb = FEXTB(loop1);
        elseif Parameter1 == 17
            c0b = C0B(loop1);
        elseif Parameter1 == 18 || Parameter1 == 14
            Jcactivef = JACTIVEF(loop1);
        elseif Parameter1 == 19
            p0b = P0B(loop1);
        elseif Parameter1 == 21
            alphaf = ALPHA(loop1);
            alphab = ALPHA(loop1);
        elseif Parameter1 == 22
            sigmaa0 = SIGMAA0(loop1);
        end
        
        dgf = dg/2;
        dgb = dg/2;
%         Jcactiveb = -Jcactivef;
        sigmaa = -sigmaa0*linspace(1,0,N)';  
       
        DF = zeros(5*N+1,5*N+1);
        Fn = zeros(5*N+1,1);
        
        % initial guess
        if loop1 == 1 && loop2 == 1
            pc = 1.5d3*ones(N,1);
            vc = -3.d1*ones(N,1);
            vn = -10*ones(N,1);
            thetan = Thetan*ones(N,1);
            cin = linspace(cinb,cinf,N)';
            v0 = 3;
        elseif loop1 == 1 && loop2 > 1
            pc = temp_pc;
            vc = temp_vc;
            vn = temp_vn;
            thetan = temp_thetan;
            cin = temp_cin;
            v0 = temp_v0;
        end
        X = [pc; vc; vn; thetan; cin; v0];
        
        iter = 0;
        ITER = true;
        while ITER
            iter = iter + 1;
            
            sigma = -sigmap0*(theta_0./thetan - thetan/2/theta_0) + sigmaa;
            dsigma = sigmap0*(theta_0./thetan.^2 + 1/2/theta_0);
            
            Fn(1) = pc(2) + sigma(2) - pc(1) - sigma(1) ...
                + dx*etan*thetan(1)*vn(1) + dx*zetac(1)*vc(1) + dx*zetan(1)*thetan(1)*vn(1) ...
                + dx*thetan(1)*nust(1)*(vn(1)+v0)/(1+((vn(1)+v0)/c2(1))^m);
            Fn(2:N-1) = pc(3:N) + sigma(3:N) - pc(1:N-2) - sigma(1:N-2) ...
                + 2*dx*etan*thetan(2:N-1).*vn(2:N-1) ...
                + 2*dx*zetac(2:N-1).*vc(2:N-1) + 2*dx*zetan(2:N-1).*thetan(2:N-1).*vn(2:N-1) ...
                + 2*dx*thetan(2:N-1).*nust(2:N-1).*(vn(2:N-1)+v0)./(1+((vn(2:N-1)+v0)./c2(2:N-1)).^m);
            Fn(N) = -alphab/alphaf/(1+alphab*dgb)*(pc(1) - p0b - R*T*cin(1) + R*T*c0b)...
                -1/(1+alphaf*dgf)*(pc(N) - p0f - R*T*cin(N) + R*T*c0f) ...
                + (dgf-alphab/alphaf*dgb)/(1+alphab*dgb)/(1+alphaf*dgf)*v0;
            Fn(N+1) = -pc(2)+pc(1) + dx*eta*thetan(1)*(vn(1)-vc(1)) - dx*zetac(1)*vc(1);
            Fn(N+2:2*N-1) = -pc(3:N)+pc(1:N-2) + 2*dx*eta*thetan(2:N-1).*(vn(2:N-1)-vc(2:N-1)) - 2*dx*zetac(2:N-1).*vc(2:N-1);
            Fn(2*N) = vc(N) - alphaf/(1+alphaf*dgf)*(pc(N)-p0f-R*T*cin(N)+R*T*c0f) ...
                + alphaf*dgf/(1+alphaf*dgf)*v0;
            Fn(2*N+1:3*N-1) = thetan(2:N).*vn(2:N) - thetan(1:N-1).*vn(1:N-1);
            Fn(3*N) = sum(thetan(1:N-1)+thetan(2:N))*dx/2/L - Thetan;
            Fn(3*N+1:4*N-1) = thetan(2:N).*vn(2:N) + vc(2:N) ...
                - thetan(1:N-1).*vn(1:N-1) - vc(1:N-1);
            Fn(4*N) = thetan(N)*vn(N) + Jactinf;
            Fn(4*N+1) = -gb*(cin(1)-c0b) + Jcactiveb +(D(1)+D(2))/2*(cin(2)-cin(1))/dx ...
                - (vc(1)+vc(2))/2*(cin(1)+cin(2))/2;
            Fn(4*N+2:5*N-1) = -(D(2:N-1)+D(1:N-2))/2/dx.*(cin(2:N-1)-cin(1:N-2)) ...
                + 1/4*(vc(2:N-1)+vc(1:N-2)).*(cin(2:N-1)+cin(1:N-2)) ...
                + (D(2:N-1)+D(3:N))/2/dx.*(cin(3:N)-cin(2:N-1)) ...
                - 1/4*(vc(2:N-1)+vc(3:N)).*(cin(2:N-1)+cin(3:N));
            Fn(5*N) = -gf*(cin(N)-c0f) + Jcactivef - (D(N)+D(N-1))/2*(cin(N)-cin(N-1))/dx ...
                + (vc(N)+vc(N-1))/2*(cin(N)+cin(N-1))/2;
            Fn(5*N+1) = (1+(Zetac*Lzeta - dg)*alphaf/(1+alphaf*dgf))*(pc(N)-p0f) - (pc(1)-p0b)...
                + sigma(N) - sigma(1) ...
                - (Zetac*Lzeta - dg)*alphaf/(1+alphaf*dgf)*R*T*(cin(N)-c0f)...
                + (etan*L*Thetan - dg + (dg - Zetac*Lzeta)*alphaf*dgf/(1+alphaf*dgf))*v0 ...
                - (etan*L*Thetan + 2*xim*L/b)*v0 - etan*L*Jactinf - Zetan*Lzeta*Jactinf  + (fextf-fextb);
            
            DF(1,[1, 2]) = [-1, 1];
            DF(1,N+1) = dx*zetac(1);
            DF(1,2*N+1) = dx*(etan + zetan(1))*thetan(1) + dx*nust(1)*thetan(1)*(1+(1-m)*((vn(1)+v0)/c2(1))^m)/(1+((vn(1)+v0)/c2(1))^m)^2;
            DF(1,3*N+1) = -dsigma(1) + dx*etan*vn(1) ...
                + dx*zetan(1)*vn(1) + dx*nust(1)*(vn(1)+v0)/(1+((vn(1)+v0)/c2(1))^m);
            DF(1,3*N+2) = dsigma(2);
            DF(1,5*N+1) = dx*nust(1)*thetan(1)*(1+(1-m)*((vn(1)+v0)/c2(1))^m)/(1+((vn(1)+v0)/c2(1))^m)^2;
            DF(N,[1, N]) = [-alphab/alphaf/(1+alphab*dgb), -1/(1+alphaf*dgf)];
            DF(N,4*N+1) = alphab/alphaf/(1+alphab*dgb)*R*T;
            DF(N,5*N) = 1/(1+alphaf*dgf)*R*T;
            DF(N,5*N+1) = (dgf-alphab/alphaf*dgb)/(1+alphab*dgb)/(1+alphaf*dgf);
            DF(N+1,[1,2]) = [1, -1];
            DF(N+1,[N+1,2*N+1]) = dx*[-eta*thetan(1)-zetac(1), eta*thetan(1)];
            DF(N+1,3*N+1) = dx*eta*(vn(1)-vc(1));
            DF(2*N,[N,2*N]) = [-alphaf/(1+alphaf*dgf), 1];
            DF(2*N,5*N) = alphaf/(1+alphaf*dgf)*R*T;
            DF(2*N,5*N+1) = alphaf*dgf/(1+alphaf*dgf);
            DF(2*N+1,[2*N+1,2*N+2]) = [-thetan(1), thetan(2)];
            DF(2*N+1,[3*N+1,3*N+2]) = [-vn(1), vn(2)];
            DF(3*N,3*N+1:4*N) = dx/2/L*[1,2*ones(1,N-2),1];
            DF(3*N+1,[N+1,N+2]) = [-1, 1];
            DF(3*N+1,[2*N+1,2*N+2]) = [-thetan(1), thetan(2)];
            DF(3*N+1,[3*N+1,3*N+2]) = [-vn(1), vn(2)];
            DF(4*N,[3*N,4*N]) = [thetan(N), vn(N)];
            DF(4*N+1,[N+1,N+2]) = -1/4*(cin(1)+cin(2));
            DF(4*N+1,4*N+1) = -gb - 1/2/dx*(D(1)+D(2)) - 1/4*(vc(1)+vc(2));
            DF(4*N+1,4*N+2) = 1/2/dx*(D(1)+D(2)) - 1/4*(vc(1)+vc(2));
            DF(5*N,[2*N-1,2*N]) = 1/4*(cin(N)+cin(N-1));
            DF(5*N,5*N-1) = 1/2/dx*(D(N)+D(N-1)) + 1/4*(vc(N)+vc(N-1));
            DF(5*N,5*N) = -gf - 1/2/dx*(D(N)+D(N-1)) + 1/4*(vc(N)+vc(N-1));
            DF(5*N+1,[1,N]) = [-1, 1-(dg - Zetac*Lzeta)*alphaf/(1+alphaf*dgf)];
            DF(5*N+1,[3*N+1,4*N]) = [-dsigma(1), dsigma(N)];
            DF(5*N+1,5*N) = (dg - Zetac*Lzeta)*alphaf/(1+alphaf*dgf)*R*T;
            DF(5*N+1,5*N+1) = - dg + (dg - Zetac*Lzeta)*alphaf*dgf/(1+alphaf*dgf) - 2*xim*L/b;
            
            
            
            for i = 2:N-1
                DF(i,[i-1, i+1]) = [-1, 1];
                DF(i,N+i) = 2*dx*zetac(i);
                DF(i,2*N+i) = 2*dx*(etan + zetan(i))*thetan(i) + 2*dx*nust(i)*thetan(i)*(1+(1-m)*((vn(i)+v0)/c2(i))^m)/(1+((vn(i)+v0)/c2(i))^m)^2;
                DF(i,[3*N+i-1, 3*N+i+1]) = [-dsigma(i-1), dsigma(i+1)];
                DF(i,3*N+i) = 2*dx*etan*vn(i) + 2*dx*zetan(i)*vn(i) ...
                    + 2*dx*nust(i)*(vn(i)+v0)/(1+((vn(i)+v0)/c2(i))^m);
                DF(i,5*N+1) =  2*dx*nust(i)*thetan(i)*(1+(1-m)*((vn(i)+v0)/c2(i))^m)/(1+((vn(i)+v0)/c2(i))^m)^2;
                DF(N+i,[i-1,i+1]) = [1, -1];
                DF(N+i,[N+i,2*N+i]) = 2*dx*[-eta*thetan(i)-zetac(i), eta*thetan(i)];
                DF(N+i,3*N+i) = 2*dx*eta*(vn(i)-vc(i));
                DF(2*N+i,[2*N+i,2*N+i+1]) = [-thetan(i), thetan(i+1)];
                DF(2*N+i,[3*N+i,3*N+i+1]) = [-vn(i), vn(i+1)];
                DF(3*N+i,[N+i,N+i+1]) = [-1, 1];
                DF(3*N+i,[2*N+i,2*N+i+1]) = [-thetan(i), thetan(i+1)];
                DF(3*N+i,[3*N+i,3*N+i+1]) = [-vn(i), vn(i+1)];
                DF(4*N+i,N+i-1) = 1/4*(cin(i)+cin(i-1));
                DF(4*N+i,N+i) = 1/4*(cin(i)+cin(i-1)) - 1/4*(cin(i)+cin(i+1));
                DF(4*N+i,N+i+1) = - 1/4*(cin(i)+cin(i+1));
                DF(4*N+i,4*N+i-1) = 1/2/dx*(D(i)+D(i-1)) + 1/4*(vc(i)+vc(i-1));
                DF(4*N+i,4*N+i) = - 1/2/dx*(D(i)+D(i-1)) + 1/4*(vc(i)+vc(i-1)) ...
                    - 1/2/dx*(D(i)+D(i+1)) - 1/4*(vc(i)+vc(i+1));
                DF(4*N+i,4*N+i+1) = 1/2/dx*(D(i)+D(i+1)) - 1/4*(vc(i)+vc(i+1));
            end
            
            DF = sparse(DF);
            X = X - DF\Fn;
            
            pc = X(1:N);
            vc = X(N+1:2*N);
            vn = X(2*N+1:3*N);
            thetan = X(3*N+1:4*N);
            cin = X(4*N+1:5*N);
            v0 = X(5*N+1);
            
            if iter > 1
                error = abs((X-temp_X)./X);
                error = sum(error)/(5*N+1);
                if error < 1d-4 || iter == Iter
                    ITER = false;
                end
            end
            
            temp_X = X;
            
            if Parameter1 == 0 && Parameter2 == 0
                if iter == Iter
                    figure(21);
                    subplot(2,1,1)
                    [hAx,hLine1,hLine2] = plotyy(x*1d6,pc,x*1d6,vc);
                    xlabel('x ({\mu}m)','fontsize',13)
                    ylabel(hAx(1),'p_c (Pa)','fontsize',13) % left y-axis
                    ylabel(hAx(2),'v_c\prime (nm/s)','fontsize',13) % right y-axis
                    set(gca,'fontsize',13); set(hAx(2),'fontsize',13);
                    set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
                    vaxis = axis; vaxis(2) = L*1d6; axis(vaxis);
                    subplot(2,1,2)
                    [hAx,hLine1,hLine2] = plotyy(x*1d6,vn,x*1d6,thetan);
                    xlabel('x ({\mu}m)','fontsize',13)
                    ylabel(hAx(1),'v_n\prime (nm/s)','fontsize',13) % left y-axis
                    ylabel(hAx(2),'\theta_n','fontsize',13) % right y-axis
                    set(gca,'fontsize',13); set(hAx(2),'fontsize',13);
                    set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
                    vaxis = axis; vaxis(2) = L*1d6; axis(vaxis);
                    figure(22)
                    plot(x*1d6,cin,'-k','linewidth',2)
                    set(gca,'fontsize',13)
                    xlabel('x ({\mu}m)','fontsize',13)
                    ylabel('c (mM)','fontsize',13);
                    vaxis = axis; vaxis(2) = L*1d6; axis(vaxis);
                end
            end
        end
        
        if loop1 == 1
            temp_pc = pc;
            temp_vc = vc;
            temp_vn = vn;
            temp_thetan = thetan;
            temp_cin = cin;
            temp_v0 = v0;
        end
                    
        Jwaterf = -alphaf/(1+alphaf*dgf)*(pc(N)-p0f-R*T*cin(N)+R*T*c0f) + alphaf*dgf/(1+alphaf*dgf)*v0;
        Jwaterb = -alphab/(1+alphab*dgb)*(pc(1)-p0b-R*T*cin(1)+R*T*c0b) - alphab*dgb/(1+alphab*dgb)*v0;
        tauf = b/2*(pc(N)+sigma(N) - p0f)-b*dgf/2*(v0-Jwaterf)+b/2*fextf;
        taub = b/2*(pc(1)+sigma(1) - p0b)+b*dgb/2*(v0-Jwaterf)+b/2*fextb;
        pstarf = p0f + dgf*(v0-Jwaterf);
        pstarb = p0b - dgb*(v0-Jwaterf);
        
        NF = b*w*thetan(N)/AF;
        ATP_actin = NF/delta*Jactinf/thetan(N)*GATP;
        ATP_ion_f = b*w*NA*abs(Jcactivef)*NC*GATP/1d9;
        ATP_ion_b = b*w*NA*abs(Jcactiveb)*NC*GATP/1d9;
        ATP_ion = ATP_ion_f + ATP_ion_b;
        
        ATP_ACTIN(loop1,loop2) = ATP_actin;
        ATP_ION(loop1,loop2) = ATP_ion;
        
        % The following calculates the energy
        v0_origin = v0/1d9;
        vn_origin = vn/1d9;
        vc_origin = vc/1d9;
        eta_origin = eta*1d9;
        zetac_origin = zetac*1d9;
        zetan_origin = zetan*1d9;
        nust_origin = nust*1d9;
        xim_origin = xim*1d9;
        dg_origin = dg*1d9;
        dgf_origin = dgf*1d9;
        dgb_origin = dgb*1d9;
        alphaf_origin = alphaf/1d9;
        alphab_origin = alphab/1d9;
        Jwaterf_origin = Jwaterf/1d9;
        Jwaterb_origin = Jwaterb/1d9;
        Jactinf_origin = Jactinf/1d9;
        Jactinb_origin = -Jactinf_origin;
        gf_origin = gf/1d9;
        gb_origin = gb/1d9;
        D_origin = D/1d9;
        Jcactivef_origin = Jcactivef/1d9;
        Jcactiveb_origin = Jcactiveb/1d9;
        
        dpf = pc(N) - p0f - dgf_origin*(v0_origin-Jwaterf_origin);
        dpb = pc(1) - p0b + dgb_origin*(v0_origin+Jwaterb_origin);
        iwaterf = alphaf_origin*R*T*(cin(N)-c0f);
        iwaterb = alphab_origin*R*T*(cin(1)-c0b);
        dphif = dpf - R*T*(cin(N)-c0f);
        dphib = dpb - R*T*(cin(1)-c0b);
        dmuf = R*T*log(cin(N)/c0f);
        dmub = R*T*log(cin(1)/c0b);
        dpsi0 = (p0f-p0b) - R*T*(c0f-c0b);
        dmu0 = R*T*log(c0f/c0b);
        Jcf = -gf_origin*(cin(N)-c0f) + Jcactivef_origin;
        Jcb = -gb_origin*(cin(1)-c0b) + Jcactiveb_origin;
        
        dvn = zeros(N,1);
        dvn(1) = (vn_origin(2)-vn_origin(1))/dx;
        dvn(N) = (vn_origin(N)-vn_origin(N-1))/dx;
        dvn(2:N-1) = (vn_origin(3:N)-vn_origin(1:N-2))/2/dx;
        Intg1 = sigmaa.*dvn;
        Intg2 = eta_origin*thetan.*(vc_origin-vn_origin).^2 + zetac_origin.*vc_origin.^2 ...
            + zetan_origin.*thetan.*vn_origin.^2 + nust_origin.*thetan.*(vn_origin+v0_origin).^2;
        logc = log(cin);
        dlogc = zeros(N,1);
        dlogc(1) = (logc(2)-logc(1))/dx;
        dlogc(N) = (logc(N)-logc(N-1))/dx;
        dlogc(2:N-1) = (logc(3:N)-logc(1:N-2))/2/dx;
        Intg3 = D_origin.*cin.*dlogc.^2;
        
        cbar = sigmap0/2*theta_0/theta_star^2 - sigmap0/2/theta_0*log(theta_star) - sigmap0/2/theta_0;
        en = sigmap0/2*theta_0./thetan + sigmap0/2*thetan/theta_0.*log(thetan) + cbar*thetan;
        den = -sigmap0/2*theta_0./thetan.^2 + sigmap0/2/theta_0*log(thetan) + sigmap0/2/theta_0 + cbar;
        e_pump = b*(dmuf*Jcactivef_origin + dmub*Jcactiveb_origin);
        e_actin = b*(pc(N)*1+den(N)+sigmaa(N)/thetan(N) - pc(1)*1-den(1)-sigmaa(1)/thetan(1))*Jactinf_origin;
        e_myosin = b*sum((Intg1(2:N)+Intg1(1:N-1))/2)*dx;
        
        e_interface = eta_origin*thetan.*(vc_origin-vn_origin).^2;
        e_interface = b*sum((e_interface(2:N)+e_interface(1:N-1))/2*dx);
        e_nucleus = zetac_origin.*vc_origin.^2 + zetan_origin.*thetan.*vn_origin.^2;
        e_nucleus = b*sum((e_nucleus(2:N)+e_nucleus(1:N-1))/2*dx);
        e_substrate = nust_origin.*thetan.*(vn_origin+v0_origin).^2;
        e_substrate = b*sum((e_substrate(2:N)+e_substrate(1:N-1))/2*dx);
        e_membrane = 2*xim_origin*L*v0_origin^2;
        
        e_flow = b*(alphaf_origin*dphif^2 + alphab_origin*dphib^2);
        e_hydraulic = b*dg_origin*(v0_origin-Jwaterf_origin)^2;
        
        Iext = b*((fextf-fextb)*v0_origin - (p0f-p0b)*v0_origin + dpsi0*Jwaterf_origin + dmu0*Jcf);
        Icell = e_pump + e_actin + e_myosin;
        Dfric = e_interface + e_nucleus + e_substrate + e_membrane;
        Dflow = e_flow + e_hydraulic;
        Dsol = b*(gf_origin*(cin(N)-c0f)*dmuf + gb_origin*(cin(1)-c0b)*dmub + R*T*sum((Intg3(2:N)+Intg3(1:N-1))/2)*dx);
        % End of energy calculation

        
        if RankCheck == 1
            DF = full(DF);
            RankCheck = (5*N+2) - rank(DF);
            if RankCheck ~= 0
                fprintf('RankCheck = %d\n',RankCheck);
            end
        end
        if Parameter1 == 0 && Parameter2 == 0
            fprintf('Jwaterf = %g nm/s\n',Jwaterf);
            fprintf('tauf = %g Pa m\n',tauf);
            fprintf('taub = %g Pa m\n',taub);
            fprintf('v0 = %g nm/s\n',v0);
        end
        
        VN(loop1,loop2) = vn(N);
        VC(loop1,loop2) = vc(N);
        V0(loop1,loop2) = v0;
        JWATERF(loop1,loop2) = Jwaterf;
        TAUF(loop1,loop2) = tauf;
        TAUB(loop1,loop2) = taub;
        PSTARF(loop1,loop2) = pstarf;
        PSTARB(loop1,loop2) = pstarb;
        PCF(loop1,loop2)  = pc(N);
        PCB(loop1,loop2)  = pc(1);
        IEXT(loop1,loop2) = Iext;
        ICELL(loop1,loop2) = Icell;
        DFRIC(loop1,loop2) = Dfric;
        DFLOW(loop1,loop2) = Dflow;
        DSOL(loop1,loop2)  = Dsol;
        E_SUBSTRATE(loop1,loop2) = e_substrate;
        E_MEMBRANE(loop1,loop2)  = e_membrane;
        E_HYDRAULIC(loop1,loop2) = e_hydraulic;
        E_PUMP(loop1,loop2) = e_pump;
        E_ACTIN(loop1,loop2) = e_actin;
        E_MYOSIN(loop1,loop2) = e_myosin;
        E_FLOW(loop1,loop2) = e_flow;
        E_INTERFACE(loop1,loop2) = e_interface;
        E_NUCLEUS(loop1,loop2) = e_nucleus;
        
        
        t_xim = xim/b;
        t_dg = dg/L;
        t_zeta = R*T/(Dcommon*(1-Thetan*0) + gf*L/2);
        t_alpha = alphaf*L;
        t_eta = (eta*Thetan+Zetac)/(1-Thetan*0);
        
        Kc = ((Thetan*nust(1)+2*t_xim+t_dg)*(2+t_alpha*(t_dg+t_eta+t_zeta*c0f)) - t_alpha*t_dg^2)...
            /(2+t_alpha*(t_dg+t_eta+t_zeta*c0f));
        fact1 = t_alpha*t_dg/(2+t_alpha*(t_dg+t_eta+t_zeta*c0f));
        v0_L = 1/Kc*((nust(1)+eta*fact1)*Jactinf + t_zeta*fact1*Jcactivef);
        V0_L(loop1,loop2) = v0_L;
        KC(loop1,loop2) = Kc;
        
    end
end

if Parameter1 == 14
    Ym = JWATERF(:,1)';
end
if Parameter2 == 14
    Xm = JWATERF(1,:);
end

%%

if Parameter1 > 0 && Parameter2 > 0
    [Xmesh,Ymesh] = meshgrid(Xm,Ym);

    figure(1)
    [~,hd] = contour(Xmesh,Ymesh,V0,'fill','on'); colorbar; hold on
    set(hd,'LevelList',linspace(min(min(V0)),max(max(V0)),35));
    [C,hd] = contour(Xmesh,Ymesh,V0,'k','linewidth',2); colorbar;  hold off
    set(hd,'LevelList',[40]);
    AxesProperties20170905
    set(gca,'fontsize',15)

end


