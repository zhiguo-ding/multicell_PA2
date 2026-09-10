%% 
clear
%figure
%close all
f = 28e9; % 28 GHz
c = 3e8; % speed of light
eta = (c/4/pi/f)^2;
lambda = 0.001;%intensity
d = 3; %waveguide height
eps = 0.5; %factional power parameter
%P = 0.1; % transmit power
sigmanoise = 10^(-12); %noise power
beta = 0.5; % LoS probability parameter
ct = 50000;
%the following parameter needs to be changed based on others
Ptxmean = 1.5170;%for rc=30, beta=0.5,eps=0.5
% %1.5539;%for rc=40, beta=0.5,eps=0.5
%2.0636;%for rc=40, beta=0.5,eps=1

%1.3620; %for rc=10, beta=0.5,eps=0.5
% 
%1.5214;%for rc=30, beta=0.5,eps=1 
% % 0.5142; %for rc=30, beta=0.5,eps=0.5
%Ptxmean = mean(mean(Ptrcontest)) / mean(mean(Ptrtest))


D_all = 1000; %overall area D_all X D_all
r_c = 30; % size of each small cluster

snrdbm = [0: 10 : 30];

for isnr = 1: length(snrdbm)
    P = 10^((snrdbm(isnr)-30)/10);
    sum_test=0;
    rate_conv_ict = zeros(ct,1);
    rate_pa_ict = zeros(ct,1);
    for ict = 1 : ct
        number_BS = poissrnd(lambda * D_all * D_all); %number of BSs
        z_BS = sign(randn(number_BS,2)).*rand(number_BS, 2) * D_all/2; %locations of BSs
                
        z_BS = [0 0; z_BS]; %add the tpical base station
        number_BS = number_BS+1;
        id_ty = 1; %the index of typical BS (user, antenna)
        
        %create poisson cluster process to find users' locations
        theta_temp = 2*pi*rand(number_BS,1);   % angle: random between 0 and 2pi
        r_temp = r_c*sqrt(rand(number_BS,1)); % random radius,sqrt to fit the pdf
        z_temp = [r_temp.*cos(theta_temp) r_temp.*sin(theta_temp)]; %generate within one circle
        z_user = z_BS + z_temp; %each user is located with a shift by its BS
        %z_user(id_ty,:) = z_BS(id_ty,:);%typical user is right underneath of its base station
        typical_user = z_user(id_ty,:);
        
        %for pinching antenna case
        %need to move the locations of BSs
        angle_temp = 2*pi*rand(number_BS,1); %random tilded angles of waveguides
        v = [cos(angle_temp), sin(angle_temp)]; %unit direction of waveguides
        t = sum((z_user-z_BS) .* v, 2); %projection scalar
        t = max(min(t, r_c), -r_c); %wrap it within the diameter
        z_PA = z_BS + t.*v; %the locations of actual transmiters
        %z_PA(id_ty,:) = z_BS(id_ty,:); %typical BS/user do not move        

        %conventional antennas%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find the PA transmit powers using fractional power
        R_user_conv = d^2 + sum((z_BS - z_user).^2,2); %the squired distance between the user and its own BS
        P_LoS_R_conv = exp(-beta*sqrt(R_user_conv));
        isLoS_R_conv = rand(number_BS,1) < P_LoS_R_conv; % the index of LoS
        isLoS_R_conv(id_ty) = 1; %typical user has LoS link
        Ptr_conv = zeros(number_BS,1);        
        Ptr_conv(isLoS_R_conv) = P * R_user_conv(isLoS_R_conv).^(2*eps/2); % the actual transmit powers LoS %P*ones(number_BS,1);%
        Ptr_conv(~isLoS_R_conv) = P * R_user_conv(~isLoS_R_conv).^(4*eps/2); % the actual transmit powers NLoS %P*ones(number_BS,1);%
        Ptr_inf_conv = Ptr_conv; Ptr_inf_conv(id_ty) = []; %contain interference lines only
        
        %find the channel gains (including the legitmate and interference channel
        %gains)
        D_user_conv = d^2 + sum((z_BS - typical_user).^2,2); %the squired distance between the typical user and all PAs
        P_LoS_conv = exp(-beta*sqrt(D_user_conv));
        isLoS_conv = rand(number_BS,1) < P_LoS_conv; % the index of LoS
        isLoS_conv(id_ty) = 1; %typical user has LoS link
        g_conv = zeros(number_BS,1);
        g_conv(isLoS_conv)  = eta./D_user_conv(isLoS_conv); 
        fading_temp_conv = exprnd(1,number_BS,1);
        g_conv(~isLoS_conv) = eta* fading_temp_conv(~isLoS_conv) ./ D_user_conv(~isLoS_conv).^2;
        g_inf_conv = g_conv; g_inf_conv(id_ty) = [];%this exclude the user's channel gain
        
        %rate_conv_ict(ict) = log(1 + Ptr_conv(id_ty)*g_conv(id_ty)/(Ptr_inf_conv'*g_inf_conv + sigmanoise));
        rate_conv_ict(ict) = log(1 + Ptr_conv(id_ty)/Ptxmean*g_conv(id_ty)/(Ptr_inf_conv'/Ptxmean*g_inf_conv + sigmanoise));
        %rate_conv_ict(ict) = log(1 + Ptr_conv(id_ty)*g_conv(id_ty)/(Ptr_inf_conv'/Ptxmean*g_inf_conv + sigmanoise));
        %rate_cons_ict(ict) = log(1 + P*r_c^(2*eps)*g_conv(id_ty)/(P*r_c^(2*eps)*ones(number_BS-1,1)'*g_inf_conv + sigmanoise)); % this uses the constant transmit power

        %pinching antenna part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %find the PA transmit powers using fractional power
        R_user = d^2 + sum((z_PA - z_user).^2,2); %the squired distance between the user and its own PA
        
        %decide the power based the existence of LoS for Ri
        P_LoS_R = exp(-beta*sqrt(R_user));
        isLoS_R = rand(number_BS,1) < P_LoS_R; % the index of LoS
        isLoS_R(id_ty) = 1; %typical user has LoS link
        Ptr = zeros(number_BS,1);
        Ptr(isLoS_R) = P * R_user(isLoS_R).^(2*eps/2); % the actual transmit powers for LoS
        Ptr(~isLoS_R) = P * R_user(~isLoS_R).^(4*eps/2); % the actual transmit powers for NLoS
        Ptr_inf = Ptr; Ptr_inf(id_ty) = []; %contain interference lines only
        
        %using mean(mean(Ptrtest)) and mean(mean(Ptrcontest))
        Ptrtest(isnr,ict) = Ptr(id_ty)/P;
        Ptrcontest(isnr,ict) = Ptr_conv(id_ty)/P;

        %find the channel gains (including the legitmate and interference channel
        %gains)
        D_user = d^2 + sum((z_PA - typical_user).^2,2); %the squired distance between the typical user and all PAs
        Dxx_conv = D_user_conv; Dxx_conv(id_ty)=[];
        Dxx = D_user; Dxx(id_ty)=[];
        zz = [sort(Dxx) sort(Dxx_conv)];
        P_LoS = exp(-beta*sqrt(D_user));
        isLoS = rand(number_BS,1) < P_LoS; % the index of LoS
        isLoS(id_ty) = 1; % typical user has LoS
        g = zeros(number_BS,1);
        g(isLoS)  = eta./D_user(isLoS); 
        fading_temp = exprnd(1,number_BS,1);
        g(~isLoS) = eta* fading_temp(~isLoS) ./ D_user(~isLoS).^2;
        g_inf = g; g_inf(id_ty) = [];%this exclude the user's channel gain
        
        rate_pa_ict(ict) = log(1 + Ptr(id_ty)*g(id_ty)/(Ptr_inf'*g_inf + sigmanoise));



    end

    %analytical results
    stepy = 2*r_c/100;   rho = P/sigmanoise;
    xupp = Inf;stepx=0.1;stept=stepx;
    sum1 = 0; sum2 = 0;

    %%%%%%
    a = rho*eta*d^(2*(eps-1));

    %outer_fun = @(t) outer_int(t,a,lambda,r_c,beta,eta,rho,d,eps,xupp);

    %ana(isnr) = integral(outer_fun, 0, Inf,'RelTol', 1e-6, 'AbsTol', 1e-9);     
    
    rate_pa(isnr) = mean(rate_pa_ict);
    rate_conv(isnr) = mean(rate_conv_ict);
    %rate_cons(isnr) = mean(rate_cons_ict);
 
end

plot( snrdbm,rate_conv,snrdbm,rate_pa )
%plot( snrdbm,rate_conv,snrdbm,rate_pa, snrdbm,ana)
 
 
 
function val = outer_int(t,a,lambda,r_c,beta,eta,rho,d,eps,xupp)

    val = zeros(size(t));

    %avoid potential singularity issue
    idx0 = (t == 0);
    idx1 = (t ~= 0);
    val(idx0) = a; % if L(0)=1, the limit is a

    % here t is a vector, so use arrayfun.
    val(idx1) = arrayfun(@(tt) ...
        ((1-exp(-a*tt))/tt) * ...
        lapl_integral2(tt,lambda,r_c,beta,eta,rho,d,eps,xupp) * ...
        exp(-tt), t(idx1));
end
 
function L = lapl_integral2(s,lambda,r_c,beta,eta,rho,d,eps,xupp)

    integrand = @(x,y) lapl_xy_integrand(x,y,s,r_c,beta,eta,rho,d,eps);
    %x for interference distance and y for own distance 

    sum1 = integral2(integrand, 0, xupp, -r_c, r_c,'RelTol', 1e-5, 'AbsTol', 1e-8);

    L = exp(-2*pi*lambda*sum1);
end

function out = lapl_xy_integrand(x,y,s,r_c,beta,eta,rho,d,eps)

    R2 = y.^2 + d^2; %squared distance for the typical user's own link 
    D2 = x.^2 + d^2; %squared distance for interference links
    R = sqrt(R2); %own link
    D = sqrt(D2); %interference link

    pR = exp(-beta*R); % LoS probabilities for own link 
    pD = exp(-beta*D); %LoS probabilities for interference  link 

    weight_y = 2*sqrt(max(r_c^2 - y.^2, 0))/(pi*r_c^2); %pdf of y 

    %there are four cases due to lo
    % Case 1: both los, TX power exponent: R2^eps, 
    % interferencechannel gain: eta/D2
    term_LL = pR.*pD .* exp( ...
        -s*rho*eta*R2.^eps ./ D2 );

    % Case 2: user NLoS, interference LoS, TX power exponent: R2^(2eps)
    % interference channel gain: eta/D2
    term_NL = (1-pR).*pD .* exp( ...
        -s*rho*eta*R2.^(2*eps) ./ D2 );

    % Case 3: USER LoS, interference NLoS, TX power exponent: R2^eps
    % interference channel gain: eta*h/D2^2, h~exp(1)
    term_LN = pR.*(1-pD) ./ ...
        (1 + s*rho*eta*R2.^eps ./ D2.^2);

    % Case 4: user NLoS, interference NLoS, TX power exponent: R2^(2eps)
    % interference channel: eta*h/D2^2, h~exp(1)
    term_NN = (1-pR).*(1-pD) ./ ...
        (1 + s*rho*eta*R2.^(2*eps) ./ D2.^2);

    qxy = term_LL + term_NL + term_LN + term_NN;

    out = (1 - qxy).*x.*weight_y;
end