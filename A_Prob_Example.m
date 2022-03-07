%clear all;
%clc;
%dx -> histogram bins
function [N1,N2,N3,N4,N5,N6,N7,N8,N9,N10]=A_Morb(dx, decimate, Moments)

  Dimension = 176; %dimmesion matrix -> days
  
  [R1_A_HOSP,X1,N1]          = Density_Analysis(Dimension, 240.15/decimate, dx);
  [R2_A_HOSP,X2,N2]          = Density_Analysis(Dimension, 186.55/decimate, dx);
  [R3_A_HOSP,X3,N3]          = Density_Analysis(Dimension, 268.20/decimate, dx);
  [R4_A_HOSP,X4,N4]          = Density_Analysis(Dimension, 360.77/decimate, dx);
  [R5_A_HOSP,X5,N5]          = Density_Analysis(Dimension, 293.93/decimate, dx);
  [R1_A_AMB,X6,N6]           = Density_Analysis(Dimension, 824.69/decimate, dx);
  [R2_A_AMB,X7,N7]           = Density_Analysis(Dimension, 1019.51/decimate,dx);
  [R3_A_AMB,X8,N8]           = Density_Analysis(Dimension, 992.95/decimate, dx);
  [R4_A_AMB,X9,N9]           = Density_Analysis(Dimension, 1600.37/decimate,dx);
  [R5_A_AMB,X10,N10]         = Density_Analysis(Dimension, 1191.83/decimate,dx);
  
  v = [length(R1_A_HOSP),length(R2_A_HOSP),length(R3_A_HOSP),length(R4_A_HOSP),length(R5_A_HOSP),length(R1_A_AMB),length(R2_A_AMB),length(R3_A_AMB),length(R4_A_AMB),length(R5_A_AMB)];
  L = max(v);
 
  L1 = Moments*max(R1_A_HOSP);
  L2 = Moments*max(R2_A_HOSP);
  L3 = Moments*max(R3_A_HOSP);
  L4 = Moments*max(R4_A_HOSP);
  L5 = Moments*max(R5_A_HOSP);
  L6 = Moments*max(R1_A_AMB);
  L7 = Moments*max(R2_A_AMB);
  L8 = Moments*max(R3_A_AMB);
  L9 = Moments*max(R4_A_AMB);
  L10 = Moments*max(R5_A_AMB);
  
  R1 = max(X1);
  R2 = max(X2);
  R3 = max(X3);
  R4 = max(X4);
  R5 = max(X5);
  R6 = max(X6);
  R7 = max(X7);
  R8 = max(X8);
  R9 = max(X9);
  R10 = max(X10);


  figure(1);
  hold on;
  plot(X1, Moments*R1_A_HOSP,'p-');
  plot(X2, Moments*R2_A_HOSP,'p-');
  plot(X3, Moments*R3_A_HOSP,'p-');
  plot(X4, Moments*R4_A_HOSP,'p-');
  plot(X5, Moments*R5_A_HOSP,'p-');
  plot(X6, Moments*R1_A_AMB);
  plot(X7, Moments*R2_A_AMB);
  plot(X8, Moments*R3_A_AMB);
  plot(X9, Moments*R4_A_AMB);
  plot(X10, Moments*R5_A_AMB);
  legend({'R1-HOSP','R2-HOSP','R3-HOSP','R4-HOSP','R5-HOSP','R1-AMB','R2-AMB','R3-AMB','R4-AMB','R5-AMB'});
  xlabel('Scaled largest Eigenvalues');
  ylabel('Probability');
  title('ASTHMA Comorbidity Distribution');
  set(gca,"fontsize",15,"linewidth",2);
  hold off;

  pathoutRadio    = "C:/DrAlberto/ProbMorbilidadesRMT/Paso2/Ver2/A_Prob.csv"
  Radio           = [[L1,R1]',[L2,R2]',[L3,R3]',[L4,R4]',[L5,R5]',[L6,R6]',[L7,R7]',[L8,R8]',[L9,R9]',[L10,R10]'];
  csvwrite (pathoutRadio,Radio);
  
endfunction