%% 
clear
%figure
%close all
f = 28e9; % 28 GHz
c = 3e8; % speed of light
eta = (c/4/pi/f)^2;
lambda = 0.001;%intensity
d = 3; %waveguide height
eps = 1; %factional power parameter
%P = 0.1; % transmit power
sigmanoise = 10^(-12); %noise power
beta = 0.5; % LoS probability parameter
ct = 50000;
Rtr = 2; %target data rate
 
D_all = 50000; %overall area D_all X D_all
r_c = 30; % size of each small cluster
[Ptxmean,Ptxmean_LoS] = power_norm(ct,r_c,beta, d,eps);% need to do power normalization

snrdbm = [0: 10 : 30];
lambda_vector = [ 0.02 0.01 0.009 0.006 0.003 0.001];

for isnr = 1: length(lambda_vector)
    lambda = lambda_vector(isnr);
    P = 10^((snrdbm(1)-30)/10);
    sum_test=0;
    rate_conv_ict = zeros(ct,1);
    rate_pa_ict = zeros(ct,1);
    sum1 = 0; sum2=0; sum3=0;
    for ict = 1 : ct
        number_BS = poissrnd(lambda * D_all); %number of BSs
        z_BS = rand(number_BS, 1) * D_all - D_all/2; %locations of BSs
        z_BS = [0; z_BS];
        number_BS = number_BS + 1;

        id_ty =  1; %typical BS close to origin
        z_temp = z_BS; z_temp(id_ty) = -D_all; %need to avoid that position to be selected
        [~, id_next] = min(abs(z_temp - z_BS(id_ty))); %find the closest BS 

        %create poisson cluster process to find users' locations
        z_temp = sign(randn(number_BS,1)).*rand(number_BS, 1) * r_c; %[-rc rc]
        z_user = z_BS + z_temp; %each user is located with a shift by its BS
        %z_user(id_ty,:) = z_BS(id_ty,:);%typical user is right underneath of its base station
        typical_user = z_user(id_ty); %typical user
        interf_user = z_user(id_next); %interference  
        

        %for pinching antenna case
        %need to move the locations of BSs: for 1D, simple
        z_PA = z_user;  
        

        %conventional antennas%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find the PA transmit powers using fractional power
        R_user_conv = d^2 + sum((z_BS - z_user).^2,2); %the squired distance between the user and its own BS
        P_LoS_R_conv = exp(-beta*sqrt(R_user_conv));
        isLoS_R_conv = rand(number_BS,1) < P_LoS_R_conv; % the index of LoS
        %israndom = 1; % typical user's link is random
        isLoS_R_conv(id_ty) = 1;israndom=0; %typical user has LoS link
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
        
        Pt_ty = Ptxmean*(israndom==1) + Ptxmean_LoS*(israndom==0);
        rate_conv_ict(ict) = log(1 + Ptr_conv(id_ty)/Pt_ty*g_conv(id_ty)...
            /(Ptr_inf_conv'/Ptxmean*g_inf_conv + sigmanoise));
        P0 = P * R_user_conv(id_ty).^(2*eps/2);%power to typical user, always LoS
        PI = P * R_user_conv(id_next).^(2*eps/2);%power to interference user, always LoS
        rate_conv_case_ict(ict) = log(1 + P0/Pt_ty*eta./D_user_conv(id_ty)...
            /(PI/Ptxmean_LoS*eta./D_user_conv(id_next) + sigmanoise));
   
        %pinching antenna part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %find the PA transmit powers using fractional power
        R_user = d^2 + sum((z_PA - z_user).^2,2); %the squired distance between the user and its own PA
        
        %decide the power based the existence of LoS for Ri
        P_LoS_R = exp(-beta*sqrt(R_user));
        isLoS_R = rand(number_BS,1) < P_LoS_R; % the index of LoS
        %isLoS_R(id_ty) = 1; %typical user has LoS link
        Ptr = zeros(number_BS,1);
        Ptr(isLoS_R) = P * R_user(isLoS_R).^(2*eps/2); % the actual transmit powers for LoS
        Ptr(~isLoS_R) = P * R_user(~isLoS_R).^(4*eps/2); % the actual transmit powers for NLoS
        Ptr_inf = Ptr; Ptr_inf(id_ty) = []; %contain interference lines only

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
        
        rate_pa_ict(ict) = log(1 + Ptr(id_ty)*g(id_ty)...
            /(Ptr_inf'*g_inf + sigmanoise));
        P0 = P * R_user(id_ty).^(2*eps/2);%power to typical user, always LoS
        PI = P * R_user(id_next).^(2*eps/2);%power to interference user, always LoS
        rate_pa_case_ict(ict) = log(1 + P0*eta./D_user(id_ty)...
            /(PI*eta./D_user(id_next) + sigmanoise));

        %outage probability
        if rate_pa_case_ict(ict)<Rtr
            sum1 = sum1 + 1; 
        end 
        if rate_conv_case_ict(ict)<Rtr
            sum3 = sum3 + 1; 
        end

        % %test for analysis
        % z0 = sign(randn)*rand*r_c;
        % zI = sign(randn)*rand*r_c;
        % r = abs(z_BS(id_ty) - z_BS(id_next));
        % a = sqrt(1/(1/d^2/(exp(Rtr)-1) - sigmanoise/eta/P/d^(2*eps))-d^2);
        % if z0-zI-r <= a & z0-zI-r >= -a
        %     sum2 = sum2 + 1;
        % end
        % 
        % zz=80;
        % if r<zz
        %     sum3 = sum3 + 1;
        % end
    end
    % sum3/ct
    % 1-exp(-2*lambda*zz)
    %%%%%%
    
    rate_pa(isnr) = mean(rate_pa_ict);
    rate_conv(isnr) = mean(rate_conv_ict);
    rate_pa_case(isnr) = mean(rate_pa_case_ict);
    rate_conv_case(isnr) = mean(rate_conv_case_ict);
    out_pa(isnr) = sum1/ct; 
    out_conv(isnr) = sum3/ct;
    out_pa_test(isnr) = sum2/ct; 

    %analytical results 
    a = sqrt(1/(1/d^2/(exp(Rtr)-1) - sigmanoise/eta/P/d^(2*eps))-d^2);
    if 0<=a & a<2*r_c
        out_ana(isnr) = a/r_c - a^2/4/r_c^2 + ...
            (2*exp(-2*lambda*a)-2+exp(-2*lambda*(2*r_c-a)) -exp(-2*lambda*(2*r_c+a)))...
            /16/lambda^2/r_c^2;
    else
        out_ana(isnr) = 1 - (exp(-2*lambda*(a-2*r_c)) - 2*exp(-2*lambda*a) + exp(-2*lambda*(a+2*r_c)))...
            /16/lambda^2/r_c^2;
    end
end

%plot( snrdbm,rate_conv,snrdbm,rate_conv_case,snrdbm,rate_pa,snrdbm,rate_pa_case )
semilogy(lambda_vector,out_conv,lambda_vector,out_pa,lambda_vector,out_ana)%,lambda_vector,out_pa_test) 
 
 
function [Ptxmean,Ptxmean_LoS] = power_norm(ct,r_c,beta,d,eps)

    for i = 1 : ct 
        Ri = rand(1, 1) * r_c;
        Ptx_conv(i) = exp(-beta*Ri)*(Ri^2+d^2)^(2*eps/2)+...
            (1-exp(-beta*Ri))*(Ri^2+d^2)^(4*eps/2);

        Ptx_conv_LoS(i) = (Ri^2+d^2)^(2*eps/2);
    end
    Ptrmean_conv = mean(Ptx_conv);
    Ptrmean_conv_LoS = mean(Ptx_conv_LoS);
    Ptxmean = Ptrmean_conv/d^(2*eps);%the distance for PA is just d
    Ptxmean_LoS = Ptrmean_conv_LoS/d^(2*eps);%the distance for PA is just d
end  