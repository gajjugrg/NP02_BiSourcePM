      program bi207PM
      real :: x
c
      drlong=205.                ! PM drift length in mm
      drshort=45.
      rainn=15.                 ! inner anode radius
      raout=30.                 ! outer anode radius
      resINNL=0.10             ! MeV 
      resOUTL=0.08              ! MeV
      resINNS=0.08              ! MeV
      resOUTS=0.08              ! MeV
      dtave=2*1000./(18.+1.4)     ! average time between events us
      srcfrac=1.4/(18.+1.4)
      drfiv=1./1.5              ! us/mm
      attdist=30000./drfiv      ! attenuation distance in mm
      thr=0.26
c     
      Aprob1=0.0015             ! ele 565
      Aprob2=0.0044+Aprob1      ! ele 556
      Aprob3=0.0152+Aprob2      ! ele 483
C                               ! gamma 570       
      EleEne1=0.482
      EleEne2=0.556
      EleEne3=0.566      
      GamEne1=0.570             ! MeV
      GamInt1=80.               ! mm
c     
      Bprob1=0.0054             ! ele 1060
      Bprob2=0.0184+Bprob1      ! ele 1049
      Bprob3=0.0703+Bprob2      ! ele 976
      Bprob4=0.7458+Bprob3      ! gamma 1064
      Bprob5=0.0687+Bprob4      ! gamma 1770
      Bprob6=0.0002+Bprob5      ! ele 1682
      Bprob7=0.0013+Bprob6      ! gamma 1440           
      
      EleEne4=1.060
      EleEne5=1.049
      EleEne6=0.976
      GamEne2=1.064             ! MeV
      GamInt2=120.              ! mm
      GamEne3=1.77              ! MeV
      GamInt3=160.              ! mm
      EleEne7=1.682
      GamEne4=1.440
      GamInt4=140.              ! mm
c      
      print *,Aprob3,Bprob7
c
      pi2=2.*3.1415927
      open(90,file="bi207stream.txt")
c
      ilosh=0
      timtot=0.
      timprv=0.
      do i=1,1000000
         call random_number(x)
      if (x .gt. srcfrac) then
            drmax=drlong
            ilosh=1
         else
            drmax=drshort
            ilosh=0
         endif
         
         dtim=-dtave*log(rndm(x))
         timtot=timtot+dtim
         
         Aprob=rndm(x)
         if(Aprob.le.Aprob1) then
            DInter=0.
            Energy=EleEne1
            iele=1
         elseif(Aprob.le.Aprob2) then
            DInter=0.
            Energy=EleEne2
            iele=1
         elseif(Aprob.le.Aprob3) then
            DInter=0.
            Energy=EleEne3
            iele=1
         else
            DInter=GamInt1
            Energy=GamEne1
            iele=0
         endif

         DinterA=DInter
         EnergyA=Energy
         ieleA=iele

         iok=0
         do while(iok.eq.0)
            x0=5.0*rndm(x)-2.5
            y0=5.0*rndm(x)-2.5
            if(x0*x0+y0*y0<6.25) iok=1
         enddo
         z0=0.
c
         d=-DinterA*log(rndm(x))
         ctheta=rndm(x)
         phi=rndm(x)*pi2         
         stheta=sqrt(1.-ctheta*ctheta)
         x1=x0+d*stheta*sin(phi)
         y1=y0+d*stheta*cos(phi)
         z1A=z0+d*ctheta
         r1A=sqrt(x1*x1+y1*y1)
c
         IokA=0
         if(EnergyA.gt.0.) then
            if(z1A.ge.0..and.z1A<drmax.and.r1A.lt.raout) then
               if(ieleA.eq.0) then
                  iok=0
                  gg=EnergyA/0.511
                  do while(iok.eq.0)
                     ct=2.*rndm(x)-1.
                     eps=1./(1.+gg*(1.-ct))
                     sctprob=0.5*eps**2*(eps+1/eps-(1.-ct**2))
                     call random_number(x)
      if (x .le. sctprob) then
                        iok=1
                        GammaE=EnergyA*eps
                        ElectEA=EnergyA-GammaE
                     endif
                  enddo
               else
                  ElectEA=EnergyA
               endif
c
               ElectEA=ElectEA*exp((z1A-drmax)/attdist)
               if(r1A.le.rainn) then
                  IokA=2
                  if(ilosh.eq.1) then
                     resol=resINNL
                  else
                     resol=resINNS
                  endif
               else
                  IokA=1
                  if(ilosh.eq.1) then
                     resol=resOUTL
                  else
                     resol=resOUTS
                  endif
               endif
               res=0.                 
               do k=1,12
                  res=res+rndm(x)
               enddo
               res=resol*(res-6.)
               rrA=ElectEA+res
c            print *,"A",i,ieleA,r1A,z1A,energyA,electeA,rrA
            endif
         else
            ElectEA=0.
            rrA=0.
         endif
            
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
c         
         Bprob=rndm(x)
         if(Bprob.le.Bprob1) then
            DInter=0.
            Energy=EleEne4
            iele=1
        elseif(Bprob.le.Bprob2) then
            DInter=0.
            Energy=EleEne5
            iele=1
         elseif(Bprob.le.Bprob3) then
            DInter=0.
            Energy=EleEne6
            iele=1
        elseif(Bprob.le.Bprob4) then
            DInter=GamInt2
            Energy=GamEne2
            iele=0
        elseif(Bprob.le.Bprob5) then
            DInter=GamInt3
            Energy=GamEne3
            iele=0
         elseif(Bprob.le.Bprob6) then
            DInter=0.
            Energy=EleEne7
            iele=1
         elseif(Bprob.le.Bprob7) then
            DInter=GamInt4
            Energy=GamEne4
            iele=0
         else
            DInter=0.
            Energy=0.
            iele=0
         endif
c
         DinterB=DInter
         EnergyB=Energy
         ieleB=iele
c    
         iok=0
         do while(iok.eq.0)
            x0=5.0*rndm(x)-2.5
            y0=5.0*rndm(x)-2.5
            if(x0*x0+y0*y0<6.25) iok=1
         enddo
         z0=0.         
         d=-DinterB*log(rndm(x))
         ctheta=rndm(x)
         phi=rndm(x)*pi2         
         stheta=sqrt(1.-ctheta*ctheta)
         x1=x0+d*stheta*sin(phi)
         y1=y0+d*stheta*cos(phi)
         z1B=z0+d*ctheta
         r1B=sqrt(x1*x1+y1*y1)
c
         IokB=0
         if(EnergyB.gt.0.) then
            if(z1B.ge.0..and.z1B<drmax.and.r1B.lt.raout) then
               if(ieleB.eq.0) then
                  iok=0
                  gg=EnergyB/0.511
                  do while(iok.eq.0)
                     ct=2.*rndm(x)-1.
                     eps=1./(1.+gg*(1.-ct))
                     sctprob=0.5*eps**2*(eps+1/eps-(1.-ct**2))
                     call random_number(x)
      if (x .le. sctprob) then
                        iok=1
                        GammaE=EnergyB*eps
                        ElectEB=EnergyB-GammaE
                     endif
                  enddo
               else
                  ElecteB=EnergyB
               endif
c            
               ElectEB=ElectEB*exp((z1B-drmax)/attdist)
               if(r1B.le.rainn) then
                  IokB=2
                  if(ilosh.eq.1) then
                     resol=resINNL
                  else
                     resol=resINNS
                  endif
               else
                  IokB=1
                   if(ilosh.eq.1) then
                     resol=resOUTL
                  else
                     resol=resOUTS
                  endif
               endif
               res=0.                 
               do k=1,12
                  res=res+rndm(x)
               enddo
               res=resol*(res-6.)
               rrB=ElecteB+res
c     print *,"B",i,ieleB,r1B,z1B,energyB,electeB,rrB
            endif
         else
            ElectEB=0.
            rrB=0.
         endif
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
c
         if(rrA.lt.thr) IokA=0
         if(rrB.lt.thr) IokB=0
         dtimA=timtot+(drmax-z1A)*drfiv
         dtimB=timtot+(drmax-z1B)*drfiv
         if(IokA.ne.0.and.IokB.ne.0) then
         if(z1B.gt.z1A) then
         write(90,*)i,ilosh,iokB,dtimB,ieleB,r1B,z1B,energyB,electeB,rrB
         write(90,*)i,ilosh,iokA,dtimA,ieleA,r1A,z1A,energyA,electeA,rrA
         else 
         write(90,*)i,ilosh,iokA,dtimA,ieleA,r1A,z1A,energyA,electeA,rrA
         write(90,*)i,ilosh,iokB,dtimB,ieleB,r1B,z1B,energyB,electeB,rrB
         endif
         elseif(IokA.ne.0) then
         write(90,*)i,ilosh,iokA,dtimA,ieleA,r1A,z1A,energyA,electeA,rrA
         elseif(IokB.ne.0) then
         write(90,*)i,ilosh,iokB,dtimB,ieleB,r1B,z1B,energyB,electeB,rrB
         endif
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
c         
      enddo
      print *,i,timtot
      end
 
