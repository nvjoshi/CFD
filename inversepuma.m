function[o]=inverse_puma(x_c,y_c,z_c)
r11=1;r12=0;r13=0;
r21=0;r22=01;r23=0;
r31=0;r32=0;r33=1;
theta1=atan2(-x_c,y_c)*180/pi;
theta2=(atan2(-x_c,y_c))*180/pi;
D=(x_c^2+y_c^2+(z_c-26.45)^2-17.5^2-17^2)/(2*17*17.5);
psi1=atan2(sqrt(1-D^2),D)*180/pi;
psi2=atan2(-1*sqrt(1-D^2),D)*180/pi;
phi1=(atan2(z_c-26.45,(sqrt(x_c^2+y_c^2)))-atan2(17.5*sin(psi1*pi/180),17+17.5*cos(psi1*pi/180)))*180/pi;
phi2=(atan2(z_c-26.45,(sqrt(x_c^2+y_c^2)))-atan2(17.5*sin(psi2*pi/180),17+17.5*cos(psi2*pi/180)))*180/pi;

b1=(atan2(sqrt(1-(-1*sin(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r12+cos(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r22-sin((phi1*pi/180+psi1*pi/180))*r32)^2),-1*sin(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r12+cos(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r22-sin((phi1*pi/180+psi1*pi/180))*r32))*180/pi;
b2=(atan2(sqrt(1-(-1*sin(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r12+cos(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r22-sin((phi1*pi/180+psi1*pi/180))*r32)^2),-1*sin(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r12+cos(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r22-sin((phi1*pi/180+psi1*pi/180))*r32))*180/pi;
b3=(atan2(sqrt(1-(-1*sin(theta2*pi/180)*cos(phi2*pi/180-psi2*pi/180)*r12+cos(theta2*pi/180)*cos(phi2*pi/180-psi2*pi/180)*r22-sin((phi2*pi/180+psi2*pi/180))*r32)^2),-1*sin(theta2*pi/180)*cos(phi2*pi/180-psi2*pi/180)*r12+cos(theta2*pi/180)*cos(phi2*pi/180-psi2*pi/180)*r22-sin((phi2*pi/180+psi2*pi/180))*r32))*180/pi;
b4=(atan2(sqrt(1-(-1*sin(theta2*pi/180)*cos(phi2*pi/180-psi2*pi/180)*r12+cos(theta2*pi/180)*cos(phi2*pi/180-psi2*pi/180)*r22-sin((phi2*pi/180+psi2*pi/180))*r32)^2),-1*sin(theta2*pi/180)*cos(phi2*pi/180-psi2*pi/180)*r12+cos(theta2*pi/180)*cos(phi2*pi/180-psi2*pi/180)*r22-sin((phi2*pi/180+psi2*pi/180))*r32))*180/pi;

d1=(atan2(cos(theta1*pi/180)*r12+sin(theta1*pi/180)*r22,-1*sin(theta1*pi/180)*sin(psi1*pi/180+phi1*pi/180)*r12-cos(theta1*pi/180)*sin(psi1*pi/180+phi1*pi/180)*r22+cos(psi1*pi/180-phi1*pi/180)*r32))*180/pi;
d2=(atan2(cos(theta1*pi/180)*r12+sin(theta1*pi/180)*r22,-1*sin(theta1*pi/180)*sin(psi2*pi/180+phi2*pi/180)*r12-cos(theta1*pi/180)*sin(psi2*pi/180+phi2*pi/180)*r22+cos(psi2*pi/180-phi2*pi/180)*r32))*180/pi;
d3=(atan2(cos(theta1*pi/180)*r12+sin(theta1*pi/180)*r22,-1*sin(theta1*pi/180)*sin(psi1*pi/180+phi1*pi/180)*r12-cos(theta1*pi/180)*sin(psi1*pi/180+phi1*pi/180)*r22+cos(psi1*pi/180-phi1*pi/180)*r32))*180/pi;
d4=(atan2(cos(theta1*pi/180)*r12+sin(theta1*pi/180)*r22,-1*sin(theta1*pi/180)*sin(psi2*pi/180+phi2*pi/180)*r12-cos(theta1*pi/180)*sin(psi2*pi/180+phi2*pi/180)*r22+cos(psi2*pi/180-phi2*pi/180)*r32))*180/pi;

a1=(atan2(-1*sin(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r11+cos(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r21-sin((phi1*pi/180+psi1*pi/180))*r31,-1*sin(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r13+cos(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r23-sin((phi1*pi/180+psi1*pi/180))*r33))*180/pi;
a2=(atan2(-1*sin(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r11+cos(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r21-sin((phi1*pi/180+psi1*pi/180))*r31,-1*sin(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r13+cos(theta1*pi/180)*cos(phi1*pi/180-psi1*pi/180)*r23-sin((phi1*pi/180+psi1*pi/180))*r33))*180/pi;

o=[theta1,phi1,psi1,d1,b1,a1;theta1,phi1,psi1,d1,b2,a1;theta1,phi2,psi2,d2,b3,a1;theta1,phi2,psi2,d2,b4,a1;theta2,phi1,psi1,d3,b1,a2;theta2,phi1,psi1,d3,b2,a2;theta2,phi2,psi2,d4,b3,a2;theta2,phi2,psi2,d4,b4,a2];

end
