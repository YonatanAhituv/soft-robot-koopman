function out1 = vf_koopman(in1,in2)
%VF_KOOPMAN
%    OUT1 = VF_KOOPMAN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    23-Oct-2024 13:07:26

u1 = in2(1,:);
u2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
t2 = u1.^2;
t3 = u1.^3;
t4 = u2.^2;
t5 = u2.^3;
t6 = x1.^2;
t7 = x1.^3;
t8 = x2.^2;
t9 = x2.^3;
et1 = t2.*(-2.20730292686887e-2)+t3.*2.320252337774905e-2-t4.*1.44058559580748e-2+t5.*7.817174470143146e-3+t6.*2.997785103448483e+2;
et2 = t7.*(-2.317735907534977e+2)-t8.*2.99739512154453e+2+t9.*2.31599993485944e+2+u1.*1.029037176678734e-1+u2.*1.026759801514546e-1;
et3 = x1.*(-2.907980088087964e+2)+x2.*2.911889925252844e+2+t2.*u2.*1.737484641333233e-2-t4.*u1.*1.361466671048331e-2+t6.*u1.*8.958424517718754;
et4 = t6.*u2.*8.927407796912714-t8.*u1.*8.609706874484884-t8.*u2.*8.609734886178233-u1.*u2.*2.978370123784272e-2-t2.*x1.*2.104086687480277e-1;
et5 = t2.*x2.*2.581181176788271e-1-t4.*x1.*2.102785344734557e-1+t4.*x2.*2.890011437768317e-1-t6.*x2.*7.737957807845797e+1+t8.*x1.*7.689258891551921e+1;
et6 = u1.*x1.*(-8.796751514477034)+u1.*x2.*8.254195485398752-u2.*x1.*8.796745395056867+u2.*x2.*8.269713314807806+x1.*x2.*1.950223557347971e-2;
et7 = u1.*u2.*x1.*(-2.104304647738505e-1)+u1.*u2.*x2.*2.88872528925434e-1+u1.*x1.*x2.*1.432616603063253e-1+u2.*x1.*x2.*1.743021321753162e-1+1.093158795354993e-2;
et8 = t2.*(-2.210213925496916e-2)-t3.*5.812719304384699e-3-t4.*2.976917840595957e-2+t5.*9.572838707749278e-3-t6.*2.997395101136605e+2;
et9 = t7.*2.314432986946653e+2+t8.*2.997785123896618e+2-t9.*2.319302863377638e+2+u1.*1.026759871323454e-1+u2.*1.029037106871543e-1;
et10 = x1.*2.911889925311807e+2-x2.*2.907980088146929e+2+t2.*u2.*1.515861255372912e-5+t4.*u1.*3.100459953499408e-2-t6.*u1.*8.630439895437824;
et11 = t6.*u2.*(-8.59942278482233)+t8.*u1.*8.937691774559177+t8.*u2.*8.937719597687838-u1.*u2.*1.439126880414027e-2+t2.*x1.*2.787001451405579e-1;
et12 = t2.*x2.*(-1.898267267368646e-1)+t4.*x1.*2.785697081828603e-1-t4.*x2.*2.207097493321984e-1+t6.*x2.*7.704928518542067e+1-t8.*x1.*7.722288397281527e+1;
et13 = u1.*x1.*8.261957513861031-u1.*x2.*8.788989509416091+u2.*x1.*8.261951315195756-u2.*x2.*8.804507428967698+x1.*x2.*1.949814996781169e-2;
et14 = u1.*u2.*x1.*2.787219342694194e-1-u1.*u2.*x2.*2.205811891373925e-1+u1.*x1.*x2.*1.84723044939398e-1+u2.*x1.*x2.*1.536826095239618e-1+1.093158795377279e-2;
out1 = [et1+et2+et3+et4+et5+et6+et7;et8+et9+et10+et11+et12+et13+et14];
end
