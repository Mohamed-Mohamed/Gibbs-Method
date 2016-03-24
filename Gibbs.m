function [ h, mag_h, i, omega, e_vector, mag_e, w, theta2, rp, zp ] = Gibbs ( r1, r2, r3, muo, R )
% this function is about Gibbs Method of determination the orbit by 3 position vectors
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% INPUTS:
% r1   :  first position vector (1x3)
% r2   :  second position vector (1x3)
% r3   :  third position vector (1x3)
% muo :  Gravitational Parameter
% R    : raduis of plant
%% OUTPUTS:
% h              : specific angular momentum vector
% mag_h     : specific angular momentum magnitude
% i               : inclination angle in degree
% omega     : right ascension of the ascending node in degree
% e              : eccentricity vector
% mag_e     : eccentricity magnitude
% w             : argument of perigee in degree
% theta2       : true anomaly of r2 in degree
% rp             : radius of perigee
% zp             : height of perigee
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
mag_r1=norm(r1);
mag_r2=norm(r2);
mag_r3=norm(r3);
C12=cross(r1,r2);
C23=cross(r2,r3);
C31=cross(r3,r1);
check=dot(r1/mag_r1,C23/norm(C23));
display(['r1_hat.C23_hat (must be = 0 ) = ' num2str(check)]);
N_vector=mag_r1*C23+mag_r2*C31+mag_r3*C12;
N_mag=norm(N_vector);
D_vector=C12+C23+C31;
D_mag=norm(D_vector);
S_vector=r1*(mag_r2-mag_r3)+r2*(mag_r3-mag_r1)+r3*(mag_r1-mag_r2);
v2_vector=sqrt(muo/N_mag/D_mag)*(cross(D_vector,r2)/mag_r2+S_vector);
v2_mag=norm(v2_vector);
% orbital element procdure
vr2=dot(r2,v2_vector)/mag_r2;
h=cross(r2,v2_vector);
mag_h=norm(h);
i=acosd(h(3)/mag_h);
N=cross([0,0,1],h);
mag_N=norm(N);
if N(2) >= 0
    omega=acosd(N(1)/mag_N);
elseif N(2) < 0
    omega=360-acosd(N(1)/mag_N);
end
e_vector=((v2_mag^2-muo/mag_r2)*r2-mag_r2*vr2*v2_vector)/muo;
mag_e=norm(e_vector);
if e_vector(3) >= 0
    w=acosd(dot(N,e_vector)/mag_N/mag_e);
elseif e_vector(3) < 0
    w=360-acosd(dot(N,e_vector)/mag_N/mag_e);
end
if vr2 >= 0
    theta2=acosd(dot(e_vector,r2)/mag_r2/mag_e);
elseif vr2 < 0
    theta2=360-acosd(dot(e_vector,r2)/mag_r2/mag_e);
end
rp=mag_h^2/muo/(1+mag_e*cosd(0));
zp=rp-R;
end