      program bi207PM
      dimension trace1(4000),trace2(4000),trace3(4000),trace4(4000)
      dimension spcOUTL(400),spcEINL(400),spcGINL(400)
      dimension spcOUTS(400),spcEINS(400),spcGINS(400)
      dimension signal(500)
c
      open(90,file="bi207stream.txt")
      open(91,file="outer_anode_long.txt")
      open(92,file="inner_anode_long.txt")
      open(93,file="outer_anode_short.txt")
      open(94,file="inner_anode_short.txt")
      open(95,file="bi207spectra.txt")
c
      do k=1,500
         x=float(k)*0.1
         signal(k)=exp(-(x-10.)/4.)/(0.56988*(1+exp(-(x-10.)/1.)))
      enddo
      
      do k=1,4000
         trace1(k)=0.
         trace2(k)=0.
         trace3(k)=0.
         trace4(k)=0.
      enddo
      
      do k=1,400
         spcOUTL(k)=0
         spcEINL(k)=0
         spcGINL(k)=0
         spcOUTS(k)=0
         spcEINS(k)=0
         spcGINS(k)=0
      enddo
c
      thr=0.26
      ioffset=1500
      timtrig=0.
      phmaxOUTL=0.
      phmaxINNL=0.
      phmaxOUTS=0.
      phmaxINNS=0.
      i=0
      ieleOUTL=-1
      ieleINNL=-1
      ieleOUTS=-1
      ieleINNS=-1
      do while(i.lt.999999)
         read(90,*,END=9999)i,ilosh,iok,dtim,iele,r1,z1,energy,electe,rr
         if(i.gt.999999) exit
         kstart=10*(dtim-timtrig)+ioffset
         if(dtim-timtrig.le.200.) then
c             print *,"R",i,ilosh,iok,dtim,rr
            if(ilosh.eq.1) then
               if(iok.eq.1) then
                  if(rr.gt.phmaxOUTL) then
                     phmaxOUTL=rr
                     ieleOUTL=iele
                  endif
                  kk=kstart
                  do k=1,500
                     kk=kk+1
                     trace1(kk)=trace1(kk)+signal(k)*rr
                  enddo
c                  print*,"B",i,"OUT-LNG",rr,kstart,ioffset
               elseif(iok.eq.2) then
                  if(rr.gt.phmaxINNL) then
                     phmaxINNL=rr
                     ieleINNL=iele
                  endif
                  kk=kstart
                  do k=1,500
                     kk=kk+1
                     trace2(kk)=trace2(kk)+signal(k)*rr
                  enddo
c                  print*,"B",i,"INN-LNG",rr,kstart,ioffset
               endif
            elseif(ilosh.eq.0) then
               if(iok.eq.1) then
                  if(rr.gt.phmaxOUTS) then
                     phmaxOUTS=rr
                     ieleOUTS=iele
                  endif
                  kk=kstart
                  do k=1,500
                     kk=kk+1
                     trace3(kk)=trace3(kk)+signal(k)*rr
                  enddo
c                  print*,"B",i,"OUT-SRT",rr,kstart,ioffset
               elseif(iok.eq.2) then
                  if(rr.gt.phmaxINNS) then
                     phmaxINNS=rr
                     ieleINNS=iele
                  endif
                  kk=kstart
                  do k=1,500
                     kk=kk+1
                     trace4(kk)=trace4(kk)+signal(k)*rr
                  enddo
c                  print*,"B",i,"INN-SRT",rr,kstart,ioffset
               endif
            endif
         endif
c
         if(dtim-timtrig.le.0.) then
            ioffset=kstart
            timtrig=dtim
         endif
c 
         if(dtim-timtrig.gt.200.) then
            phOUTL=0.
            phINNL=0.
            phOUTS=0.
            phINNS=0.
            if(ioffset.gt.1800.or.ioffset.le.0)
     $           PRint *,">>>>>>>",ioffset,"<<<<<<<<<"
            do k=ioffset,ioffset+2200
               if(trace1(k).gt.phOUTL) phOUTL=trace1(k)
               if(trace2(k).gt.phINNL) phINNL=trace2(k)
               if(trace3(k).gt.phOUTS) phOUTS=trace3(k)
               if(trace4(k).gt.phINNS) phINNS=trace4(k)
            enddo
            if(phmaxOUTL.gt.thr) then
               write(91,*) i,ieleOUTL,timtrig,phmaxOUTL,phOUTL
               k=phOUTL*200.
               if(k.lt.400) spcOUTL(k)=spcOUTL(k)+1
            endif
            if(phmaxINNL.gt.thr) then
               write(92,*) i,ieleINNL,timtrig,phmaxINNL,phINNL
               k=phINNL*200.
               if(k.lt.400) then
                  if(ieleINNL.eq.1) then
                     spcEINL(k)=spcEINL(k)+1
                  else
                     spcGINL(k)=spcGINL(k)+1
                  endif
               endif
            endif
            if(phmaxOUTS.gt.thr) then
               write(93,*) i,ieleOUTS,timtrig,phmaxOUTS,phOUTS
               k=phOUTS*200.
               if(k.lt.400) spcOUTS(k)=spcOUTS(k)+1
            endif
            if(phmaxINNS.gt.thr) then
               write(94,*) i,ieleINNS,timtrig,phmaxINNS,phINNS
               k=phINNS*200.
               if(k.lt.400) then
                  if(ieleINNS.eq.1) then
                     spcEINS(k)=spcEINS(k)+1
                  else
                     spcGINS(k)=spcGINS(k)+1
                  endif
               endif
            endif
c            print *,"W",phmaxOUTL,phmaxINNL,phmaxOUTS,phmaxINNS
c            print *,"W",phOUTL,phINNL,phOUTS,phINNS
c     print *,"R",i,ilosh,iok,dtim,rr

            timtrig=dtim
            ioffset=1500
            kstart=ioffset
            phmaxOUTL=0.
            phmaxINNL=0.
            phmaxOUTS=0.
            phmaxINNS=0.
            do k=1,4000
               trace1(k)=0.
               trace2(k)=0.
               trace3(k)=0.
               trace4(k)=0.
            enddo
            ieleOUTL=-1
            ieleINNL=-1
            ieleOUTS=-1
            ieleINNS=-1
c     
            kk=kstart
            if(ilosh.eq.1) then
               if(iok.eq.1) then
                  phmaxOUTL=rr
                  ieleOUTL=iele
                  do k=1,500
                     kk=kk+1
                     trace1(kk)=signal(k)*rr
                  enddo
c                  print*,"A",i,"OUT-LNG",rr,kstart,ioffset
               elseif(iok.eq.2) then
                  phmaxINNL=rr
                  ieleINNL=iele
                  do k=1,500
                     kk=kk+1
                     trace2(kk)=signal(k)*rr
                  enddo
c                  print*,"A",i,"INN-LNG",rr,kstart,ioffset
               endif
            elseif(ilosh.eq.0) then
               if(iok.eq.1) then
                  phmaxOUTS=rr
                  ieleOUTS=iele
                  do k=1,500
                     kk=kk+1
                     trace3(kk)=signal(k)*rr
                  enddo
c                  print*,"A",i,"OUT-SRT",rr,kstart,ioffset
               elseif(iok.eq.2) then
                  phmaxINNS=rr
                  ieleINNS=iele
                  do k=1,500
                     kk=kk+1
                     trace4(kk)=signal(k)*rr
                  enddo
c                  print*,"A",i,"INN-SRT",rr,kstart,ioffset
               endif
            endif
         endif
      enddo
 9999 do k=1,400
         write(95,*) k*0.005-0.0025,
     $        spcOUTL(k),spcEINL(k),spcGINL(k),spcEINL(k)+spcGINL(k),
     $        spcOUTS(k),spcEINS(k),spcGINS(k),spcEINS(k)+spcGINS(k)
      enddo
      close(90)
      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      end
