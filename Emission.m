%% new instance of TMat Script
close all
clc
clear


[Em_data]=xlsread('Em_data_pce.xls');   %load in emission data file
wEm=size(Em_data,2);     %amount of columns = amount of different thicknesses +1
lEm=size(Em_data,1);     %amount of rows = amount of different wavelengths +1
thickness_range=150:25:150;  %Em_data(1,2:wEm)
lambda_range=760; %(Em_data(2:50:lEm,1))'

for lambda_var=lambda_range
for thickness_var=thickness_range  

lambda=535;                 %incident light wavelength
lambda_emitted=round(lambda_var);
spectrum=400:1:1000;            %examined spectrum, in nm
%% Solar Panel Model
% standard model = 5 layers: glass, ITO, PEDOT/PSS, the active layer (PCPDTBT/PC60BM) and the aluminum 
% (PEDOT/PSS is poly(3,4-ethylenedioxythiophene)/poly(styrene sulfonate) and ITO is indium tin oxide)

% initiate layer parameters
layers = {'Air' 'Glass' 'PCD3' 'Air'}; % Names of layers of materials starting from side light is incident from
                                                                           % Names must match the strings in the excel n,k source file
thicknesses = [7000000 thickness_var];                                            % thickness of each corresponding layer in nm
theta_in=0*pi/180;                                                                                                                                         
active_layer=2; %in reference to thicknesses

% Load in index of refraction for each material
ntotal = zeros(size(layers,2),size(spectrum,2));
for index = 1:size(layers,2)
    ntotal(index,:) = LoadRefrIndex(layers{index},spectrum);
end

% Constants
h = 6.62606957e-34; 	% Js Planck's constant
c = 2.99792458e8;	% m/s speed of light
q = 1.60217657e-19;	% C electric charge 

%% Calculate the R,T for the smaller stack

% initiate layer parameters

layers2 = {'Glass' layers(3) 'Air'}; % Names of layers of materials starting from side light is incident from
                                                                           % Names must match the strings in the excel n,k source file
thicknesses2 = [thicknesses(active_layer)];                                            % thickness of each corresponding layer in nm
theta_in2=0*pi/180;                                                                                                                                         
active_layer2=1; %in reference to thicknesses

% Load in index of refraction for each material
ntotal2 = zeros(size(layers2,2),size(spectrum,2));
for index = 1:size(layers2,2)
    ntotal2(index,:) = LoadRefrIndex(layers2{index},spectrum);
end

posLam=find(spectrum==lambda_emitted);
n_extra=ntotal(:,posLam:posLam);

t=thicknesses;
t2 = thicknesses2;
ls=length(spectrum);

[~,T_stack] = Tmat(lambda_emitted,n_extra,t2,theta_in2);

%% Absorbance and Generation
% absorption coefficient a gives the fraction of incident
% radiant energy absorbed per unit thickness (usually cm-1 but here taken m-1)
% The absorbed intensity in one layer is therefore int(a*I(x)dx)
%I=some constant*n*E2 but we use n*E2(x) instead of I(x) since absorbance is a
%relative value (I/I0)

% % Load in 1sun AM 1.5 solar global tilt in mW/cm2
% AM=xlsread('AM15.xls');                                 

i=0;
for lam=lambda:lambda 
    i=lam-spectrum(1)+1;
    n = ntotal(:,i:i).';
    [R,T,E2_lam,Eps,Ems,Epp,Emp] = Tmat(lam,n,t,theta_in);  %Calculate new E2 for every lambda
    EppLam(lam)=Epp;
    EmpLam(lam)=Emp;
    EpsLam(lam)=Eps;
    EmsLam(lam)=Ems;
%   note AM15 is left out at this point

        a=4*pi*imag(n(active_layer+1))/lam/1e-9; %in metres (usually in cm)
        for z=0:t(active_layer)
            alpha(z+1)=a*real(n(active_layer+1))/real(n(1))*E2_lam(z+1);  %calculate Qj for a certain position/layer
                                                                %note that alpha=Iabs/I0 means we can leave out constants, 
                                                                %but we still need to include n of I0. |E+0|=1
            G(z+1,lam)=a*real(n(active_layer+1))*E2_lam(z+1) *lam*1e-9/h/c; %# of generated pairs in the active layer per m2, at zi
                                                                                  %G=a*E2*nl /(hf) = a*E2*nl /(hc)*lam
                                                                                  %to see why I=AM15*n*E2, see notes p.13
        end
end
 

%% Current density
%Calculate J(x) from G(x)

% Model boundaries for layers and get the total generation rate over all lambda
for z=0:t(active_layer)
    Gtot(z+1)=sum(G(z+1,:))*1E-9; %*1E-9 because per nm
end

xJ=0:t(active_layer);
% figure(4)
% hold on
% plot(xJ,G(:,lambda),'LineWidth',2)
% ylabel('G(z) [#pairs]')
% xlabel('z')
% title('Generated pairs/(s.nm.m^2) in active layer','FontSize',12)

Jsc=sum(Gtot(:))*1e-9*q * 0.1; %in mA/cm^2 (q = 1.60217657e-19;	% C electric charge)


%% Emission Model

%Implement the Gtot (all wavelengths) f(z) profile to radiative
%recombination. For now assuming arbitrary chance of radiative recombination. For the
%position, calculate an expected average lifetime, multiply by mobility to
%find average travelled distance. Maybe for this you need the dipole model.

%%Diffusion equations solver for R(z)
Dn=8.5/10*10^14; %in nm^2/s 
t_RNR=2*10^-12; %in s, this is for PM6, 2 picoseconds

Diff_length=sqrt(Dn*t_RNR) %in nm
%note in this case diffusion length is rather long, typically it is 5-20nm
%in organics

%inititate some differential equation variables
Elps=EpsLam(lambda);
Elms=EmsLam(lambda);
Elpp=EppLam(lambda);
Elmp=EmpLam(lambda);
posLam=find(spectrum==lambda);
k=imag(ntotal(active_layer+1,posLam:posLam));
n=real(ntotal(active_layer+1,posLam:posLam));
B=a*n*lambda*1e-9/h/c*1E-9; %factor beta before E2 in G(x), twice *1E-9 because of photon energy and because per nm

% figure(7)
% for z=0:t(active_layer)
%     Gtest(z+1)=1/2*(B*1E9*(abs(Elps)^2*exp(-4*pi*k*z/lambda)+abs(Elms)^2*exp(4*pi*k*z/lambda)...
%     +2*real(Elps*conj(Elms))*cos(4*pi*n*z/lambda)-2*imag(Elps*conj(Elms))...
%     *sin(4*pi*n*z/lambda))+B*1E9*(abs(Elpp)^2*exp(-4*pi*k*z/lambda)+abs(Elmp)^2*exp(4*pi*k*z/lambda)...
%     +2*real(Elpp*conj(Elmp))*cos(4*pi*n*z/lambda)-2*imag(Elpp*conj(Elmp))...
%     *sin(4*pi*n*z/lambda)));
% end
% zz=0:t(active_layer);
% plot(zz,Gtest)

syms y(z)
ode = -Dn*diff(y,z,2)+y/t_RNR == 1/2*(B*1E9*(abs(Elps)^2*exp(-4*pi*k*z/lambda)+abs(Elms)^2*exp(4*pi*k*z/lambda)...
    +2*real(Elps*conj(Elms))*cos(4*pi*n*z/lambda)-2*imag(Elps*conj(Elms))...
    *sin(4*pi*n*z/lambda))+B*1E9*(abs(Elpp)^2*exp(-4*pi*k*z/lambda)+abs(Elmp)^2*exp(4*pi*k*z/lambda)...
    +2*real(Elpp*conj(Elmp))*cos(4*pi*n*z/lambda)-2*imag(Elpp*conj(Elmp))...
    *sin(4*pi*n*z/lambda)));
Dy = diff(y,z);
DDy=diff(Dy,z);
cond1 = Dy(0) == 0;
cond2 = Dy(t(active_layer)) == 0;
conds = [cond1 cond2];
ySol(z) = dsolve(ode,conds);
R=matlabFunction(ySol);

% figure(7)
% fplot(ySol/t_RNR,[0 t(active_layer)])
figure(5)
plot(xJ,R(xJ)/t_RNR,'LineWidth',2)
ylabel('R(z) [excitons/(nm.m^2)/tau]')
xlabel('z')
title('Recombination per nm.m^2 under equilibrium','FontSize',11)
drawnow

for i=1:thickness_var
    RR(i,1)=i;
    RR(i,2)=R(thickness_var-i+1);
end

%%Plug in the Recombination profile in Emission model

%new lambda, to be modified by Stokes shift
lambda_old=lambda;
lambda=lambda_emitted;
posLam=find(spectrum==lambda);
Emission_z=zeros(1,t(active_layer));
n=ntotal(:,posLam:posLam);

 theta_det=8/180*pi;
% %theta calculation (Snell's law), where theta_det is the angle of detector
% %switch glass/x
ntemp=n(2);
n(2)=real(n(3));
n(3)=ntemp;

theta=zeros(1,length(n));
theta(1)=theta_det;
for i=2:length(n)
    theta(i)=asin(n(i-1)/n(i)*sin(theta(i-1)));
end

t(active_layer)=t(active_layer)/abs(cos(theta(2)));

omeg=pi; %the solid angle (placeholder)

%Fresnel coefficients calculation
rxg_s=(-n(3)*cos(theta(3))+n(2)*cos(theta(2)))/(n(3)*cos(theta(3))+n(2)*cos(theta(2)));
rxg_p=(-n(3)*cos(theta(2))+n(2)*cos(theta(3)))/(n(3)*cos(theta(2))+n(2)*cos(theta(3)));

Rga_s=abs((n(2)*cos(theta(2))-n(1)*cos(theta(1)))/(n(2)*cos(theta(2))+n(1)*cos(theta(1))))^2;
Rga_p=abs((n(2)*cos(theta(1))-n(1)*cos(theta(2)))/(n(2)*cos(theta(1))+n(1)*cos(theta(2))))^2;

Txg_s=1/real(n(3))*abs(2*n(2)*cos(theta(2))/(n(3)*cos(theta(3))+n(2)*cos(theta(2))))^2;
Txg_p=1/real(n(3))*abs(2*n(2)*cos(theta(2))/(n(3)*cos(theta(2))+n(2)*cos(theta(3))))^2;

rxa_s=(n(2)*cos(theta(2))-n(1)*cos(theta(1)))/(n(2)*cos(theta(2))+n(1)*cos(theta(1)));
rxa_p=(n(2)*cos(theta(1))-n(1)*cos(theta(2)))/(n(1)*cos(theta(2))+n(2)*cos(theta(1)));

txa_s=2*n(2)*cos(theta(2))/(n(1)*cos(theta(1))+n(2)*cos(theta(2)));
txa_p=2*n(2)*cos(theta(2))/(n(1)*cos(theta(2))+n(2)*cos(theta(1)));

%run over all positions, mapping R(z) on the updated thickness
for z=1:t(active_layer)

    t0=z;
    t1=t(active_layer)-z;
    tF=t0+t1;

    E_fwd_s=1/2/omeg*sqrt(cos(theta(2))*R(round(z*cos(theta(2)))))*(1+rxg_s*exp(1i*4*pi*n(2)/lambda*t0))/(1-rxg_s*rxa_s*exp(1i*4*pi*n(2)/lambda*tF));
    E_end_s(z)=E_fwd_s*exp(1i*2*pi*n(2)/lambda*t1);
    
    E_fwd_p=1/2/omeg*sqrt(cos(theta(2))*R(round(z*cos(theta(2)))))*(1+rxg_p*exp(1i*4*pi*n(2)/lambda*t0))/(1-rxg_p*rxa_p*exp(1i*4*pi*n(2)/lambda*tF));
    E_end_p(z)=E_fwd_p*exp(1i*2*pi*n(2)/lambda*t1);
    
    E_rev_s=1/2/omeg*sqrt(cos(theta(2))*R(round(z*cos(theta(2)))))*(1+rxa_s*exp(1i*4*pi*n(2)/lambda*t1))/(1-rxa_s*rxg_s*exp(1i*4*pi*n(2)/lambda*tF));
    
    E_rev_p=1/2/omeg*sqrt(cos(theta(2))*R(round(z*cos(theta(2)))))*(1+rxa_p*exp(1i*4*pi*n(2)/lambda*t1))/(1-rxa_p*rxg_p*exp(1i*4*pi*n(2)/lambda*tF));
    
    I_incoh(z)=1/4*(abs(E_rev_s)^2*Txg_s*Rga_s+abs(E_rev_p)^2*Txg_p*Rga_p); %*c*eps*n(4)/2
end

Emission_total_s=sum(E_end_s);
Emission_total_p=sum(E_end_p);
I_incoh_total=sum(I_incoh)*T_stack; %T_stack is wrong here but negligible anyway
I_total=(I_incoh_total+abs(1/2*Emission_total_s*txa_s)^2+abs(1/2*Emission_total_p*txa_p)^2)*cos(theta(1)) %*c*eps*n(4)/2, final cos = lambertian

% figure(6)
% hold on
% plot(1:t(active_layer),(real(E_end_s)+real(E_end_p))/2,1:t(active_layer),(imag(E_end_s)+imag(E_end_p))/2,'LineWidth',2)
% ylabel('End field contribution [a.u]')
% xlabel('Layer position [nm] (Glass interface = 0nm)')
% title('E Field at end from each point (NOT total Emission E field)') % that's right... you are not counting interference of all the E fields! 
% legend('real','imag')

% pos_t=find(Em_data(1,:)==thickness_var)
% pos_l=find(Em_data(:,1)==lambda_var)

% data=Em_data(pos_l,pos_t)


%Emission_adjusted(lambda_emitted,pos_t)=data/I_total;

clear('I_incoh')
clear('E_end_s')
clear('E_end_p')
end
end

%This is why you need a dipole model or at the very least an electric field
%model (which means restoring the M{} in TmatEM which you cut for speed.
%Note however that you only need to count E fields from first emission,
%since the chance for Rad recomb is so low that recycling is much lower,
%and the big incident field also cannot cancel the small emitted field,
%which is what matters. So you need to calculat E+ and E- again and follow
%the same procedure for the actual Inner field. 
%Another reason for the dipole model is the higher field when close to
%surface! This might cause discrepancy with Setfos even after interference
%implementation

%======================================================================================================================================================
%======================================================================================================================================================

%Other functions, used above

%======================================================================================================================================================
%======================================================================================================================================================
%% Tmat function code
% basic Transformation Matrix function where the input is a wavelength and two vectors
% (refractive indices and lengths of the respective layers)
% Returns Reflectance, Transmittance, and the D matrix and t vector necessary for calculating the E-field
function [R,T,Mstore,tstore,Efact] = Tmat_old(lambda,n,L) 
lay=length(L);                              % lay is the number of layers
for i=1:lay+3                                   %convert n to k vectors, note that surrounding medium indices are included as first and last
    k(i)=2*pi*n(i)/lambda;
end

%incoherent first layer: glass layer R and T
tGLASS=2*n(1)/(n(1)+n(2));      
tGLASSREV=2*n(2)/(n(1)+n(2));
rGLASS=(n(1)-n(2))/(n(1)+n(2));
T_glass=real(n(2)/n(1))*abs(tGLASS)^2;           %Note the capitalization difference between intensity R,T and E-field r,t
T_glassrev=real(n(1)/n(2))*abs(tGLASSREV)^2;        %also the index ratios to account for intensity in different media
R_glass=abs(rGLASS)^2;
% T_glass=1;                            %turn off incoherent layer
% T_glassrev=1;
% R_glass=0;

Mtot=eye(2);                                %Mtot will be the transfer matrix of the entire stack
Mstore{lay+2}=eye(2);                           %Mstore stores the transfer matrix at each layer for E-field distribution (Dj in paper)

for i=1:lay+1                             %calculate the total transfer matrix to be returned, Abeles formalism, p-polarized, thetaj=0
    if(i==1)
        d(1)=0;
    else
        d(i)=k(i+1)*L(i-1);                 %pay attention to indices, k(i) and n(i) have extra entries for surrounding layers   
    end
    rp(i)=(n(i+1)-n(i+2))/(n(i+1)+n(i+2));
    tp(i)=2*n(i+1)/(n(i+1)+n(i+2));
    M{i}=Mcalc(rp(i),d(i));            
    Mtot= (Mtot*M{i});                    %note that Mtot does not include the tp(i) yet!
    if(i==1)
        tTOTa=tp(i);
    else
        tTOTa=tTOTa*tp(i);
    end
    tstore(i)=tTOTa;                        %stores ti for E-field, elements are t1, t1*t2, t1*t2*t3, ...
end
for i=1:lay+1
    Mstore{lay+2-i}=M{lay+2-i}*Mstore{lay+3-i};     %creates D matrices for Efield
end

%Calculate T and R
tSTACK=tTOTa/Mtot(1,1);
T_stack=real(n(lay+3)/n(2))*abs(tSTACK)^2; %prefix of n/n because comparing intensities in bottom(air) VS in top(glass)

rSTACK=Mtot(2,1)/Mtot(1,1);
R_stack=abs(rSTACK)^2;

R=R_stack/(1-R_stack*R_glass)*T_glass*T_glassrev+R_glass;   %series expansion of incoherent layer internal reflections
T=T_stack*T_glass/(1-R_stack*R_glass);

Efact=tGLASS/sqrt(1-R_stack*R_glass);  %Efield normalization coefficient after incoherent layer (allowed because incoherent)                             
end

%% Transfer matrix function
%calculate one transfer matrix (1 layer)
%Abeles formalism is used
function M = Mcalc(r,d)                 
                                        
M=ones(2);
M(1,1)=exp(-j*d);
M(1,2)=r*exp(-j*d);
M(2,1)=r*exp(j*d);
M(2,2)=exp(j*d);
end

%% Efield function code
function [Epz,Emz,Ep,Em] = Efield(Mstore,tstore,lambda,n,L,Efact)
lay=length(L);                      % lay is the number of layers
for i=1:lay+2                                   %inputs get partitioned to seperate n,k vectors
    k(i)=2*pi*n(i)/lambda;
end

%calculate the total distance of layers
Ltot=0;
Lsum(1)=0;
for i=1:lay
    Ltot=Ltot+L(i);
    Lsum(i+1)=Ltot;           %Lsum makes sure z runs correctly
end
%initiate E-fields
Epz=linspace(0,0,Ltot+1);
Emz=linspace(0,0,Ltot+1);
for i=1:lay
     Ep(i)=Efact*tstore(i)*Mstore{i+1}(1,1)/Mstore{1}(1,1); %Mstore{1} is simply Mtot, E+0 is Efact(E+ at the glass/ITO) assuming Eincident is 1 in air
     Em(i)=Efact*tstore(i)*Mstore{i+1}(2,1)/Mstore{1}(1,1);                                          
     for z=Lsum(i):Lsum(i+1)
         Epz(z+1)=Ep(i)*exp(1i*k(i+2)*(z-Lsum(i)));           %again pay attention to index of k, we need the k in the layer                                                     
         Emz(z+1)=Em(i)*exp(-1i*k(i+2)*(z-Lsum(i)));
     end
end
end

%% Function LoadRefrIndex 
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library2.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function ntotal = LoadRefrIndex(name,wavelengths)

%Data in IndRefr, Column names in IndRefr_names
[IndRefr,IndRefr_names]=xlsread('Optical_Constants.xls');

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr(:,strmatch(strcat(name,'_lambda'),IndRefr_names));

n=IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
k=IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));   

Nan=find(isnan(file_wavelengths));
if(~isempty(Nan))
    file_wavelengths=file_wavelengths(1:Nan(1)-1);
    n=n(1:Nan(1)-1);
    k=k(1:Nan(1)-1);
end

% Interpolate/Extrapolate data linear*Ly to desired wavelengths
n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
ntotal = n_interp+1i*k_interp; 
end

function [EQEe] = LoadEQE(name,wavelengths)

[IndRefr,IndRefr_names]=xlsread('EQE.xls');
file_wavelengths=IndRefr(:,strmatch(strcat(name,'_lambda'),IndRefr_names));

p=IndRefr(:,strmatch(strcat(name,'_p'),IndRefr_names));

Nan=find(isnan(file_wavelengths));
if(~isempty(Nan))
file_wavelengths=file_wavelengths(1:Nan(1)-1);
p=p(1:Nan(1)-1);
end

% Interpolate/Extrapolate data linear*Ly to desired wavelengths
p_interp=interp1(file_wavelengths, p, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
EQEe = p_interp.';
end