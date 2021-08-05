      include 'nlib.f'
      include 'olib.f'
      include 'clib.f'
      include 'resub.f'
      include 'rlib.f'
      include 'tlib.f'
      include 'flib.f'  
      include 'getxz1.f' 
      subroutine meemum_fun(pname_o,props_o,psys_o,iPs,ido,v_i,cblk_i)       

      implicit none
     
      include 'perplex_parameters.h'

      integer i, ier, idead

      logical bulk, bad

      character amount*6, yes*1

      integer itri(4), jtri(4), ijpt, ierr, nodata, l, j

      double precision wt(3), num ,icount,isnan
c----------------------------------------------------------------------
c                                 these common blocks are necessary to 
c                                 communicate between the MAIN program
c                                 and perplex:

c                                 iwt => 1 mass bulk comp, 0 molar bulk comp
      integer iwt
      common/ cst209 /iwt

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c                                 atwt -> g-formula wts of the chemical components
      double precision atwt
      common/ cst45 /atwt(k0) 
c                                 phase and system properties loaded by 
c                                 subroutine getloc, the call to getloc
c                                 is unnecessary if only the molar phase
c                                 proportions and compositions are of interest.
c                                 phase (prop(i,1:ntot)) and system (psys(i))
c                                 are for index i (this list may not be exhaustive):
c                                  1  - molar volume
c                                  2  - molar enthalpy
c                                  3  - gruneisen thermal parm
c                                  4  - K_S
c                                  5  - Mu_S
c                                  6  - v_phi
c                                  7  - v_p
c                                  8  - v_s
c                                  9  - v_p/v_s
c                                  10 - rho
c                                  11 - G
c                                  12 - cp
c                                  13 - alpha
c                                  14 - beta
c                                  15 - S
c                                  16 - molar amount
c                                  17 - molar weight
c                                  18 - KS_T
c                                  19 - MuS_T
c                                  20 - KS_P
c                                  21 - MuS_P
c                                  22 - vphi_T
c                                  23 - vp_T
c                                  24 - vs_T
c                                  25 - vphi_P
c                                  26 - vs_P
c                                  27 - vp_P
c                                  28 - heat capacity ratio (cp/cv)
      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)
c                                 v -> the standard perplex potential variables
      double precision v
      common/ cst5  /v(l2)
      double precision tr,pr,r,ps
      common/ cst5  /tr,pr,r,ps
c                                 ipot -> the number of potential variables in use
c                                 jv -> the indices of the potential variables
      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
c                                 vname -> the name of the potential variables
      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)
c                                 cname -> the names of the components
      character*5 cname
      common/ csta4 /cname(k5)
c                                 phase names
      character pname*14
      common/ cxt21a /pname(k5)
c                                 phase compositions
      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

c                                 jbulk -> the number of chemical components

      integer jbulk
c                                 cblk -> the molar bulk composition
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
c                                 ntot - number of phases stable
c                                 np - number of solution phases
c                                 ncpd - number of compounds
c                                 kkp(ntot) - pointer to cpd or solution model
c                                 cp3(1:jbulk,1:ntot) - molar phase compositions
      integer kkp, np, ncpd, ntot
      double precision cp3, amt
c      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                 iam -> a variable indicating the Perple_X program
      integer iam 
      common/ cst4 /iam

c----------------local variable declaration for NEW STUFF------------------
      double precision cblk_i(k5)!! k5 = 12
      double precision v_i(l2)!!l2 = 5
      double precision props_o(k5,9+k5) !! array of numerical output properties per phase
      double precision psys_o(1,9+k5) !! array of numerical output properties of bulk system     
      character*14 pname_o(k5)
      double precision iPs, ido

      integer, save :: iprim = 0

      props_o=0d0
      psys_o=0d0
      ips=0d0
      v=0d0
      cblk=0d0
      pname = ' '
      pname_o = ' '
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 2
      rxn = .false.
c                                 initialization, read data files etc. 

      
      if(iprim == 0) then
          call iniprp
c          write(9990,*) iprim
          iprim = 1
c          else
c          write(9991,*) iprim
      endif
  
      bulk = .true. !!            bulk is true, composition and p-t conditions
      amount = 'weight'
      icount = -10e0
c                                 computational loop
c      print *, ier
      v = v_i
      cblk = cblk_i
c      write(*,'(a)') 'start of assigned v & cblk'
c      print *, v
c      print *, cblk  
c      write(*,'(a)') 'end of assigned v & cblk'
         if (bulk) then 
            if (iwt.eq.1) then 
c                                 convert mass to molar 
               do i = 1, jbulk
                  cblk(i) = cblk(i)/atwt(i)
               end do 
            end if
c            ctotal = 0d0
c            do i = 1, icp
c               ctotal = ctotal + cblk(i)
c            end do 

c            do i = 1, icp
c               b(i) = cblk(i)/ctotal
c            end do

         end if 

        call meemum (bad,idead)

c ------------next lines are equivalent to meemum (bad)-----
c         call incdp0 !! set dependent variables
c         call lpopt0 (idead) !! minimization and outputs the results (to the print file???)

c         if (idead.gt.0) then
c             write(*,667) v(jv(1)), v(jv(2)) !! can't print to command line in MATLAB
c         else
c            call getloc (itri,jtri,ijpt,wt,nodata) !! compute derivative properties
c         end if 
c ------------------up till here -----------------------
         if (goodc(1)+badc(1).gt.0d0) then

            num = badc(1)/(badc(1)+goodc(1))*1d2
            if (num.gt.1d-1) call warn (53,num,i,'MEEMUM')
            goodc(1) = 0d0
            badc(1) = 0d0 

         end if   
c----------------NEW STUFF------------------
c-----------Computing array of numerical output properties per phase
!! why this post-process is done inside the do loop???
c         write(*,*) 'Number of stable phases', ntot
         do i = 1, ntot                               
            props_o(i,1) =  props(17,i)*props(16,i)/psys(17)*1d2 !! weight %
            props_o(i,2) =  props(1,i)*props(16,i)/psys(1)*1d2 !! vol %   
            props_o(i,3) =  props(16,i)/psys(16)*1d2 !! mol %
            props_o(i,4) =  props(17,i) !! molar weight
	    props_o(i,5) =  props(2,i)  !! molar enthalpy
	    props_o(i,6) =  props(15,i) !! S
	    props_o(i,7) =  props(1,i)  !! molar volume
	    props_o(i,8) =  props(12,i) !! cp
            props_o(i,9) =  props(13,i) !! alpha
	    props_o(i,10)=  props(10,i) !! rho
	    props_o(i,11)=  props(7,i)  !! v_p
	    props_o(i,12)=  props(8,i)  !! v_s            
	    do l = 1,icomp
              props_o(i,12+l) =  pcomp(l,i)*atwt(l)/props(17,i)*1d2 !! wt% composition
            end do       
         end do

c----------Computing array of numerical output properties of bulk system
	psys_o(1,1) = v(jv(2))
	psys_o(1,2) = v(jv(1))
	psys_o(1,3) = psys(11) !! G
	psys_o(1,4) = psys(17) !! molar weight
	psys_o(1,5) = psys(2)  !! molar enthalpy
	psys_o(1,6) = psys(15) !! S
	psys_o(1,7) = psys(1)  !! molar volume
	psys_o(1,8) = psys(12) !! cp
	psys_o(1,9) = psys(13) !! alpha
	psys_o(1,10) = psys(10) !! rho
	psys_o(1,11) = psys(7)  !! v_p
	psys_o(1,12)= psys(8)  !! v_s
        do j = 1, icomp
	  psys_o(1,12+j) = fbulk(j)*atwt(j)/psys(17)*1d2 !! wt% composition 
        end do  
        
c        write(*,'(a)') 'Props out'
c        print *, props_o  
c        write(*,'(a)') 'Psys out'
c        print *, psys_o     

c      print *,   psys_o

c 667   format(/,'minimization failed at this point   ',2(f8.4,4x),/)
       pname_o = pname
       iPs = ntot
       ido = idead
c      write(*,'(a)') 'Thermodynamic phases'
c      print *, pname_o

c      print *,   pname_o
      return
      end
