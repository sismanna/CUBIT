!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The program celle subdivide the space of catalogue in cell containing each one neve events !!!!!!!!!
!!!! 		   It looks for the largest earthquake j not yet assigned to a cell             !!!!!!!!!
!!!!	   Then looks for all the events around the j_th event within a circle of radius d      !!!!!!!!!
!!!!        If the number of events is larger than neve, d is reduced by a quantity delta       !!!!!!!!!
!!!!      If the number of events is smaller than neve, d is increased by a quantity delta      !!!!!!!!!
!!!!  delta is adjusted on the base of the average interdistance between the events in the cell !!!!!!!!!
!!!!      the research is performed for ntry times						!!!!!!!!!
!!!!                     Each cell is charcterized by its index ke			        !!!!!!!!!
!!!!      Each event in the catalogue is marked with the cell index through the vector iflag    !!!!!!!!!
!!!!  The procedure continues up to when all the events in the catalogue are assignd to a cell  !!!!!!!!!
!!!!The procedure ends when the number of events not yet assigned to a cell is smaller than neve!!!!!!!!!
!!!!  The exit of the program is a new catalogue of events with their coordinates 		!!!!!!!!!
!!!!             the cell index, the completeness magnitude and the b value     		!!!!!!!!!
!!!!                   attributed to each earthquake						!!!!!!!!!
!!!! the completeness magnitude is estimated evaluating the cv and uses a thresold value in the !!!!!!!!!
!!!!                                 setup file                                                 !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



parameter(k=3000000,ncelmax=1000)
dimension x(k),y(k),z(k),q(k),iflag(k),n(ncelmax),qdummy(k)
double precision d,dr,x,y,z,p,dr1,dr2,dd,df,xr,yr,zr,xx,yy,zz,dist,xs,ys,zs,xm,ym,zm
double precision to(k),tdummy(k),xdummy(k),ydummy(k),zdummy(k)
character*30 inp,outp,outc

p=3.1415926535879/180. !this is the pi value
p=p*1.d+0

!reads from setup the name of the input, the name of the output and the parameters neve, ntol, delta,ntry,the initial radius

open(50,file='setupmore')

read(50,*)inp  !inpunt file name
read(50,*)outp  !output file name for all the catalogue events
read(50,*)outc  !output file name for the individuated cells
read(50,*)neve  !number of events per cell
read(50,*)ntol  !tolerance of the number of events per cell
read(50,*)delta  !spatial increment distance
read(50,*)ntry  !maximum number of times for looking the cell
read(50,*)dini  !initial distance should be order of magnitude of the analysed area (km)
read(50,*)qrange  !minimum magniotude range for the b value estimation
 close(50)
 
open(80,file=inp)
open(90,file=outp)
open(70,file=outc)

!!!!  reads the data from the opened file 
!!!!  sets iflag(i) and find the maximum event in the catalogue 

do i=1,3000000
  read(80,*,end=99)to(i),x(i),y(i),z(i),q(i)
  q(i)=q(i)+0.1*(rand()-.5)
!  q(i)=.23+.86*q(i)
  iflag(i)=0
end do
99 npt=i-1
 close(80)

call magminmax(q,npt,iflag,qmintot,qmax,imax)

if(qmintot>0.)qmintot=0.

do i=1,npt
  q(i)=q(i)-qmintot !shift magnitude if minimum magnitude is negative for estimating the completeness magnitude
end do

ke=0
ntot=0
do ic=1,100000000   !loop for looking all the possible cells in the catalogue 
  ieve=imax
  ke=ke+1  !number of cell
  iflag(ieve)=ke
  xr=x(ieve)*p
  yr=y(ieve)*p
  zr=z(ieve)
  d=dini   !!!! d is the initial radius !!!!!
  
!!!! looks for the events in a cell of radius d around the imax event
  dd=0.
  call optrad(x,y,z,xr,yr,zr,iflag,npt,ntol,neve,ntry,d,delta,dd) 

  if(dd==0.)then
    call magminmax(q,npt,iflag,qmin,qmax,imax)
    cycle
  end if
  

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! assigns all the events to the cell !!!!!

  n(ke)=1
  qdummy(n(ke))=q(ieve)
  xdummy(n(ke))=x(ieve)
  ydummy(n(ke))=y(ieve)
  zdummy(n(ke))=z(ieve)
  tdummy(n(ke))=to(ieve)    
  qminc=10.
  qminc=10.
  xs=xdummy(1)
  ys=ydummy(1)
  zs=zdummy(1)
  do j=1,npt
    if(iflag(j).ne.0)cycle
    xx=x(j)*p
    yy=y(j)*p
    zz=z(j)
    dr=dist(xx,yy,zz,xr,yr,zr)
    if(dr<=dd)then
      iflag(j)=ke
      n(ke)=n(ke)+1
      qdummy(n(ke))=q(j)
      xdummy(n(ke))=x(j)
      ydummy(n(ke))=y(j)
      zdummy(n(ke))=z(j)
      tdummy(n(ke))=to(j) 
      xs=xs+xdummy(n(ke)) 
      ys=ys+ydummy(n(ke)) 
      zs=zs+zdummy(n(ke)) 
      if(q(j)<qminc)qminc=q(j)      
    end if
  end do
  xm=xs/dfloat(n(ke))
  ym=ys/dfloat(n(ke))
  zm=zs/dfloat(n(ke))  
  print*,ke,'cell found'
!  if(q(imax)-qminc<qrange)then
!    print*,'troppo piccolo',ke
!    go to 10
!  end if
  
  print*,n(ke)

  call subb(qdummy,ke,n(ke),b,sb)
  
  write(70,*)ke,b,sb,n(ke),xm,ym,zm
  
!writes on the output file all the events

  do i=1,n(ke)
    write(90,200)tdummy(i),xdummy(i),ydummy(i),zdummy(i),ke,qdummy(i)+qmintot,qc,b,sb
  end do
  print*,ke,'cell done'

!ending the program

10  ntot=ntot+n(ke)
  print*,ntot,npt
  if(npt-ntot<neve+ntol)exit !!!! if the number of events not assigned to a cell is smaller than neve-ntol the program ends


!!!!  looks for the next largest event

  qmax=-10.
  do j=1,npt
    if(iflag(j).ne.0)cycle
    if(q(j)>qmax)then
      qmax=q(j)
      imax=j
    end if
  end do
!  iflag(imax)=ke
end do

200 format(e12.6,1x,3(e12.6,1x),i4,1x,e10.4,2(1x,f5.2),1x,e10.4)
300 format(i3,1x,f4.2,1x,f5.3,1x,i3,1x,3(e12.6,1x))
end

!**********************************************************************
!calculates the distance between the point (xx,yy,zz) and the point (xr,yr,zr)

function dist(xx,yy,zz,xr,yr,zr)
double precision dr,xr,yr,zr,df,xx,yy,zz,dist,prad
   df=abs(yy-yr)
   dr=dacos(dsin(xx)*dsin(xr)+dcos(xx)*dcos(xr)*dcos(df))
   dr=dr*6370.d+0
   dist=dsqrt(dr**2.+(zr-zz)**2.)
   
return
end

!**********************************************************************
!optimize the radius of the cell for including neve+-ntol earthquaekes in the cell
!The inputs are:
!x,y,z = coordinates of the earthquakes in the catalogue
!xr,yr,zr = coordinates of the largest event in the cell assumed as the center of the cell
!iflag = the vector identifying the cell of the corresponding earthquake
!neve, ntol = the number of events per cell and the admitted tolerance
!ntry = the number of times for enlarging or reducing the cell radius
!d = the starting radius in km
!delta = the factor multiplyng the average distance between the centre of the cell and the earthquakes already included in the cell
!!!!
!The output is dd the radius allowing the inclusion of neve+-ntol events in the cell
!!!!

subroutine optrad(x,y,z,xr,yr,zr,iflag,npt,ntol,neve,ntry,d,delta,dd) 
parameter(k=3000000)
double precision x(k),y(k),z(k),dr,d,xx,yy,zz,dist,xr,yr,zr,dd,p
dimension iflag(k)
p=3.1415926535879/180. !this is the pi value
p=p*1.d0

  do iv=1,ntry
    ne=0
    ds=0.
    do j=1,npt
      if(iflag(j).ne.0)cycle
      xx=x(j)*p
      yy=y(j)*p
      zz=z(j)
      dr=dist(xx,yy,zz,xr,yr,zr)
      if(dr<d)then
        ds=ds+dr 
        ne=ne+1
      end if
    end do
    if(ne==0)then
      d=2.*d  !if the initial radius is too much short the program doubles it
      cycle
    end if
    dm=ds/float(ne)  !dm is the average distance between the centre of the cell and the earthquakes already included
    s=dm*delta
    if(ne==1)s=5.  !if the number of events is 1 increases d by a large amount here fixed at 5.
    if(neve-ntol<=ne.and.ne<=neve+ntol)then   !!!!  exits from the loop when ne=neve+-ntol
      dd=d
!      d=.2d+0
      delta=.1
      exit
    end if
    if(ne<neve)d=d+s 
    if(ne>neve)d=d-s  
  end do
return
end

!**********************************************************************
!find the largest and smallest event in the catalogue
subroutine magminmax(q,nq,iflag,qmin,qmax,imax)
dimension q(nq),iflag(nq)

qmin=10.
qmax=-10.
do i=1,nq
  if(iflag(i).ne.0)cycle
  if(q(i)<qmin)qmin=q(i)
  if(q(i)>qmax)then
    qmax=q(i)
    imax=i
  end if
end do
return
end

!**********************************************************************
!estimates the completeness magnitude using the cv method 
!Godano, C. and Petrillo, G. and Lippiello, E. 'Evaluating the incompleteness magnitude using an unbiased estimate of the $b$ value.'
!accepted on Geophysical Journal International, 2023

subroutine subcv(q,nq,qmax,qmin,cvth,qce)
dimension q(nq),q1(nq),iflag(nq)


!call magminmax(q,nq,iflag,qmin,qmax,imax)

qc=qmin
do iqc=1,55
  qc=qc+.1
  if(qc>=qmax-.2)exit
  qs=0.
  k=0
  do i=1,nq
    if(q(i)<qc)cycle
    k=k+1
    q1(k)=q(i)-qc
    qs=qs+q(i)-qc
  end do
  qm=qs/float(k)
  qqs=0.
  do i=1,k
    qqs=qqs+(q1(i)-qm)**2.
  end do
  sd=sqrt(qqs/float(k-1))
  cv=sd/qm
  if(cv>=cvth)then  !cv=0.93 is a threshold value fixed according to synthetic tests
    qce=qc
    exit
  end if
end do

return
end

!***********************************************************

subroutine subb(q,ke,nn,bs,sbs)
dimension q(nn),ngr(150),dm(nn)

do i=1,150
  ngr(i)=0
end do

dms=0.
k=0.
do i=1,nn-100
  do j=1,100
    dmm=q(i+j)-q(i)
    if(dmm>0.)then
      k=k+1
      dm(k)=dmm
      exit
    end if
  end do
  l=int(dmm*10.)
  ngr(l)=ngr(l)+1
  dms=dms+dmm
end do

dqm=dms/float(k)
bs=1./(dqm*log(10.))

qqs=0.
do i=1,k
  qqs=qqs+(dm(i)-dqm)**2.
end do
de=float(k)*float(k-1)
sbs=bs*bs*sqrt(qqs/de)

if(ke==8)print*,'88888888888',k,de,sbs

lu=ke*100
do i=1,150
  if(ngr(i)>0)write(lu,*)i*.1,float(ngr(i))/float(k)
end do

return    
end

