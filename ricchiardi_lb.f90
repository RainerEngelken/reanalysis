! -*- coding: utf-8 -*-
!
! ricchardi_lb.py
! 
! The efficent  implementation of ricchardi integral for the
! rate of LIF neurons.
!
! Copyright (C) 2015 Carl van Vreeswijk and Farzad Farkhooi
!
! This code is provided based on GNU GPLv3 license.
! Please concult with farzad@bccn-berlin.de in case problem running it.



function ricciardi_lb(mu,sigma,taum,taur,VT,Vr,Vl)

implicit none
real*8, parameter ::sqpi=1.7724538509055159d0
real*8 :: ricciardi_lb
real*8 :: mu,sigma,taum,taur,VT,Vr,Vl
real*8 :: xm,xp,xl,gm,gp,hl,rat
real*8 :: f_riccilb,g_riccilb,h_riccilb

if(sigma==0.d0)then
 if(mu<VT)then
  ricciardi_lb=0.d0
 else
  ricciardi_lb=1.d0/(taur+taum*log((mu-Vr)/(mu-VT)))
 endif
 return
endif

xp=(mu-Vr)/sigma
xm=(mu-VT)/sigma
xl=(mu-Vl)/sigma
ricciardi_lb=f_riccilb(xp)-f_riccilb(xm)
gm=g_riccilb(xm)
gp=g_riccilb(xp)
hl=h_riccilb(xl)
if(xm>0.d0)then
  rat=xm*exp(gm-gp)/xp
  ricciardi_lb=ricciardi_lb-sqpi*xp*exp(gp+hl+log(1.d0-rat))
elseif(xp>0.d0)then
  ricciardi_lb=ricciardi_lb-sqpi*xm*(2.d0-exp(hl))*exp(gm)
  ricciardi_lb=ricciardi_lb-sqpi*xp*exp(gp+hl)
elseif(xl>0.d0)then
  rat=xp*exp(gp-gm)/xm
  ricciardi_lb=ricciardi_lb-sqpi*xm*(2.d0-exp(hl))*exp(gm+log(1.d0-rat))
else
  rat=xp*exp(gp-gm)/xm
  ricciardi_lb=ricciardi_lb-sqpi*xm*exp(gm+hl+log(1.d0-rat))
endif
ricciardi_lb=1.d0/(taur+taum*ricciardi_lb)

return
end function ricciardi_lb

function f_riccilb(x)

implicit none
real*8 :: f_riccilb,x
real*8 :: x1,z

x1=abs(x)
z=x1/(1+x1)
f_riccilb=log(2*x1+1)+z*(-2.2757881388024176d-1+z*( &
  & 7.7373949685442023d-1+z*(-3.2056016125642045d-1+z*( &
  & 3.2171431660633076d-1+z*(-6.2718906618071668d-1+z*( &
  & 9.3524391761244940d-1+z*(-1.0616084849547165d0 +z*( &
  & 6.4290613877355551d-1+z*(-1.4805913578876898d-1)))))))))

return
end function f_riccilb

function g_riccilb(x)
implicit none
real*8 :: x,z,g_riccilb

z=x*x/(x*x+4.5d0)
g_riccilb=x*x-log(2*x*x+1)+z*(1-z)*( &
  &  5.9990862761060271+z*(-27.478249116288669 +z*( &
  &  136.63730568542917+z*(-734.98742766587850 +z*( &
  &  3370.5293232494396+z*(-12111.995350989970 +z*( &
  &  31233.658292386182+z*(-50400.802949023360 +z*( &
  &  25835.064338634256+z*( 93030.171040421192 +z*( &
  & -273478.58569046418+z*( 377392.98460630426 +z*( &
  & -313425.63518440089+z*( 158938.36390696367 +z*( &
  & -45315.549863323235+z*5551.8485546204729)))))))))))))))
return
end function g_riccilb

function h_riccilb(x)

real*8 :: x,t,z,h_riccilb

z=abs(x)
t=1./(1.d0+5.d-1*z)
h_riccilb=-z*z+log(t)-1.26551223+t*(1.00002368+t*(.37409196+&
  & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+&
  & t*(1.48851587+t*(-.82215223+t*.17087277))))))))

return
end function h_riccilb
