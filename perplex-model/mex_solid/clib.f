c routines common to all programs? could be in tlib.f?

      subroutine setau1 (output)
c----------------------------------------------------------------------
c setau1 sets autorefine dependent parameters. called by vertex, werami,
c pssect, convex, and meemum.

c output is set to false if autorefine mode is not auto (i.e., iopt(6) = 2) 
c or it is auto and in the second cycle.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical output
 
      character y*1, badnam(h9)*10

      integer ibad2, ibad1, igood, i, j, ier

      character*100 n10nam,n11nam,n8nam

      character*100 prject,tfname
      common/ cst228 /prject,tfname
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      character tname*10
      logical refine, resub
      common/ cxt26 /refine,resub,tname

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      integer iam
      common/ cst4 /iam
c-----------------------------------------------------------------------
      refine = .false.
c                                 only use autorefine if solutions
c                                 are present and it is requested.
      if (isoct.ne.0) then 

         call mertxt (n10nam,prject,'.arf',0)
         open (n10, file = n10nam, iostat = ier, status = 'old')

         call mertxt (n8nam,prject,'.tof',0)

         if (iam.eq.1.or.iam.eq.2.or.iam.eq.15) then
c                                 VERTEX, MEEMUM, or CONVEX:
            if (iam.eq.1.or.iam.eq.15) then 

               open (n8, file = n8nam, status = 'unknown')
c                                 user friendly text version 
               if (lopt(11)) then 
                  call mertxt (n11nam,prject,'_auto_refine.txt',0)
                  open (n11, file = n11nam, status = 'unknown')
               end if 

            end if 

            ibad1 = 0 

            if (ier.ne.0.and.(iam.eq.1.or.iam.eq.15)) then 
c                                 no auto_refine data
               open (n10, file = n10nam, status = 'unknown')

            else if (ier.eq.0.and.(iam.eq.1.or.iam.eq.15)) then 

               read (n10,*,iostat=ier) ibad1, ibad2, igood
               if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)

               if (iopt(6).ne.2.or.output) write (*,1030) n10nam

               if (iopt(6).eq.1) then 
c                                 manual mode, allow reinitialization
c                                 or suppression.
                  write (*,1060) 
                  read (*,'(a)') y

                  if (y.eq.'y'.or.y.eq.'Y') then

                     iopt(6) = 0

                  else 

                     refine = .true.  

                  end if

                  output = .true.
 
               else if (output) then  
c                                 second cycle of automated mode
                  refine = .true.

               end if  

               write (n8,*) refine

            else if (ier.eq.0.and.iam.eq.2) then 
c                                 MEEMUM, ask the user if he wants
c                                 to use the data 
               write (*,'(/,a,a,/,a)') 'Auto-refine data exists from a',
     *                  ' previous calculation with VERTEX.',
     *                   'Do you want MEEMUM to use this data (y/n)?'
               read (*,'(a)') y

               if (y.ne.'y'.and.y.ne.'Y') then

                  iopt(6) = 0

               else 

                  refine = .true.  
                  read (n10,*,iostat=ier) ibad1, ibad2, igood
                  if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)
                  iopt(6) = 1

                  write (*,1030) n10nam

               end if

            end if 
c                                 set cycle dependent parameters
            if (refine.or.iam.eq.2) then 

               i = 2 

            else 

               i = 1

            end if
c                                 solvus tolerance 
            if (lopt(9)) nopt(8) = 1.5d0*rid(3,i)
c                                 number of iterations
            iopt(10) = grid(6,i)
c                                 bound relaxation rate
c                                 the initial resolution of the exploratory stage
            nopt(10) = rid(3,i) 
c                                 speciation tolerance
            nopt(5) = rid(5,i)

         else if (iam.eq.13) then
c                                 the global level of unsplt, which should generate 
c                                 neither arf or tof files; at this point if ier is
c                                 zero, the arf file has been successfully opened:
c                                 kill and close it:
            if (ier.ne.0) close (n10, status = 'delete')
c                                 open and kill the tof file
            open (n8, file = n8nam, status = 'unknown')
            close (n8, status = 'delete')
c                                 open and kill the irf file
            call mertxt (n8nam,prject,'.irf',0)
            open (n8, file = n11nam, iostat=ier, status = 'unknown')
            close (n8,status = 'delete')

         else
c                                 werami/pssect if refine, get the 
c                                 solution models to be rejected
            open (n8, file = n8nam, iostat=ier, status = 'old')
        
            if (ier.eq.0) then 
c                                 write a flag to indicate if auto-refine
c                                 has been used, this is necessary so that other
c                                 perplex programs know whether to reject the
c                                 badnam phases:
               read (n8,*,iostat=ier) refine
c                                 read phases to be rejected if in auto-refine
               if (refine) then 
                  read (n10,*,iostat=ier) ibad1, ibad2, igood
                  if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)
               end if 

            end if 

         end if
c                                 only want *_auto_refine.txt for the exploratory
c                                 stage. VERTEX or CONVEX:
         if (refine) then

            lopt(11) = .false.

         else if (.not.refine.and.(iam.eq.1.or.iam.eq.15)) then
c                                 user friendly text version of the exploratory stage
c                                 auto_refine file:
            if (lopt(11)) then 
               call mertxt (n11nam,prject,'_auto_refine.txt',0)
               open (n11, file = n11nam, status = 'unknown')
            end if
c                                 write blurb
            write (n11,1000)

         end if 

      end if 

      close (n8)
c                                 just to be sure
      if (iopt(6).eq.0) refine = .false.

      if (refine) then 
c                                 reject solution models that were 
c                                 not found to be stable and set parameters 
c                                 that depend on refinement
         ibad2 = 0 

         do 50 i = 1, isoct

            do j = 1, ibad1
               if (fname(i).eq.badnam(j)) then
                  if (iam.eq.1.or.iam.eq.15) write (*,1070) fname(i)
                  goto 50
               end if 
            end do 

            ibad2 = ibad2 + 1
            fname(ibad2) = fname(i)

50       continue 

         isoct = ibad2 

         write (*,'(/)')

      end if

      if (iopt(6).eq.2.and..not.refine) then
c                                 this means it must be in the exploratory
c                                 stage
         output = .false.

      else

         output = .true.

      end if

      if (.not.(iopt(6).eq.2.and.refine).and.iopt(34).ne.0.and.
     *    iam.eq.1) then
c                                  initialize (i.e., delete prior) intermediate 
c                                  results file 
         call mertxt (n11nam,prject,'.irf',0)
         open (1000, file = n11nam, iostat=ier, status = 'unknown')
         close (1000,status = 'delete')

      end if

1000  format (//,'NOTE: this file echoes the auto-refine data after ',
     *       'the exploratory stage. If',/,'the composition of a phase',
     *       ' has been relaxed (**warning ver993**) during this stage,'
     *    /,'then the appropriate subdivision scheme* should be modifi'
     *      ,'ed and the exploratory',/,'stage calculation repeated un'
     *      ,'til the warning has been eliminated. This process can be',
     *     /,'expedited by setting the auto_refine option = man or off',
     *    //,'For the less critical compositional ranges at the end of',
     *       ' the auto-refine stage refer',/,'to the console output ',
     *       'written to the console at the end of the calculation.',//,
     *      '*see the header of the solution model file for a brief ex',
     *       'planation of subdivision schemes',//)
1030  format (/,'Reading data for auto-refinement from file: ',a,/)
1060  format ('Suppress or reinitialize auto-refinement (y/n)?')
1070  format ('Eliminating solution model: ',a,' in auto-refinement.')

      end 

      subroutine setau2 (output)
c----------------------------------------------------------------------
c setau2 sets/resets autorefine parameters after the solution models have
c been read. setau1 must be called first.

c output is set to true if autorefine mode is auto (i.e., iopt(6) = 2) 
c but no solutions are present (isoct = 0). 
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical output

      integer i,index
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      logical oned
      common/ cst82 /oned

      character tname*10
      logical refine, resub
      common/ cxt26 /refine,resub,tname
c-----------------------------------------------------------------------
      if (isoct.eq.0) then 
     
         index = 2
         output = .true.

      else if (.not.output) then

         index = 1

      else 

          if (refine) then

             index = 2

          else 

             index = 1

          end if 

      end if 
c                                 set auto-refine dependent parameters
      if (icopt.eq.5) then 
c                                 gridded minimization
         if (oned) then 
            jlow = grid(4,index)
            loopx = 1
         else 
            jlow = grid(2,index)
            loopx = grid(1,index) 
         end if

         jlev = grid(3,index) 
          
      else if (icopt.gt.5) then 
c                                 1d/2d phase fractionation
         jlow = grid(4,index)

      else if (icopt.eq.1) then 
c                                 schreinemakers diagrams

c                                 max variance of curves to be traced
          isudo = grid(5,index)
c                                 default variable tracing increment
          do i = 1, 2
             dv(iv(i)) = (vmax(iv(i)) - vmin(iv(i)))*rid(1,index)
          end do 

      else if (icopt.eq.3) then 
c                                 mixed variable diagrams 

c                                 no variance restriction
          isudo = 99
c                                 default search increment
          dv(iv(1)) = (vmax(iv(1)) - vmin(iv(1)))*rid(1,index)

      end if 

      end 

      subroutine input1 (first,output,err)
c-----------------------------------------------------------------------
c input1 reads data from a file on unit n1, this data controls the
c computational options and is modified frequently.

c iam - indicates calling program 1 - vertex
c                                 2 - meemum
c                                 3 - werami
c                                13 - unsplt, global call
c                                14 - unsplt, local call
c                                 any other values no output
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      logical output, eof, first, err

      character*100 blank*1,string(3)*8,rname*5,name*8,strg*80,n2name,
     *              n9name,y*1,sname*10,prt*3,plt*3

      integer idum, nstrg, i, j, k, ierr, icmpn, jcont, kct

      double precision dip

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      character*100 cfname
      common/ cst227 /cfname

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      character*162 title
      common/ csta8 /title(4)

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character*8 xname,vname
      common/ csta2 /xname(k5),vname(l2)

      character*5 cname
      common/ csta4 /cname(k5) 

      integer icp2
      common/ cst81 /icp2

      character*5 zname
      common/ cst209a /zname

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)

      double precision buf
      common/ cst112 /buf(5)

      integer iwt
      common/ cst209 /iwt

      integer ivfl
      common/ cst102 /ivfl

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer ibuf,hu,hv,hw,hx 
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer isoct
      common/ cst79 /isoct

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      logical oned
      common/ cst82 /oned

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      character*8 eoscmp
      common/ cst98 /eoscmp(2)

      integer iam
      common/ cst4 /iam

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      save blank
      data blank/' '/
c-----------------------------------------------------------------------
c                             output = .false. then in 1st cycle of
c                             autorefine.
      if (.not.output) then 
c                                 read computational option file
c         call fopen1 


c        write(*,*) ' the name of the input file is in clib'

c           open (n1,file='test_bob_aux.dat',iostat=ierr,status='old')
c           open (n1,file='test_kelly.dat',iostat=ierr,status='old')
c           open (n1,file='input.dat',iostat=ierr,status='old')
           open (n1,file='test_solid.dat',iostat=ierr,status='old')
c           open (n1,file='test_marthe.dat',iostat=ierr,status='old')
c           open (n1,file='test_bob_FeO.dat',iostat=ierr,status='old')
c           open (n1,file='test_trinity_H_TiMnNi.dat',iostat=ierr,status='old')
c           open (n1,file='test_trinity_JH_TiMnNi.dat',iostat=ierr,status='old')
c           open (n1,file='test_trinity_H_Ti.dat',iostat=ierr,status='old')
c           open (n1,file='test_trinity_JH_Ti.dat',iostat=ierr,status='old')
c           open (n1,file='test_trinity_H.dat',iostat=ierr,status='old')
c           open (n1,file='test_trinity_JH.dat',iostat=ierr,status='old')
c           open (n1,file='test_bob_O2.dat',iostat=ierr,status='old')
           
c           open (n1,file='test_bob.dat',iostat=ierr,status='old')
c           open (n1,file='test_BOB_hp622.dat',iostat=ierr,status='old')
c	 dummy name for 'internal' files
	 prject = 'my_project'
      
      else 
c                                 create the file name
         call mertxt (tfname,prject,'.dat',0)
         open (n1, file = tfname, iostat = ierr, status = 'old')
         if (ierr.ne.0) call error (120,r,n1,tfname)

      end if 
c                                 begin reading input:

c                                 read name of thermodynamic data file
      read (n1,'(a)') n2name
c      write (*,1100) n2name
      
c1100  format (' thermodynamic data file ', a)
      call enblnk (n2name)
c                                 read print and graphic file names
      read (n1,'(a)') prt

      read (n1,'(a)') plt

      read (n1,'(a)') n9name
      
      
      call enblnk (n9name)
c
      do i = 1, 4
         title(i) = ' '
      end do 
c                                 read title for the calculation:
      read (n1,'(a)') title(1)
c                                 read computational option or option file name
c                                 use error condition to determine which:
      read (n1,'(a)') tfname
c                                 get first non-blank string 
      call getstg (tfname)

      read (tfname,'(i2)',iostat=ierr) icopt 

      if (ierr.eq.0) then 
c                                 if no error, old version
         tfname = 'perplex_option.dat'

      else
c                                 new version, read icopt
         read (n1,*,err=998) icopt

      end if 
c                                 if fractionation path from data 
c                                 file, get name:
      fileio = .false.

      if (icopt.eq.10.or.icopt.eq.11) then 

         fileio = .true.

         read (n1,'(a)') cfname
         call enblnk (cfname)

         if (icopt.eq.10) then 
            icopt = 7
         else 
            icopt = 9
         end if

      end if 
c                                 if meemum, override whatever computational option
c                                 is set in the input file. 
      if (iam.eq.2) icopt = 5
c                                 dummy variable place holders
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum

      read (n1,*,err=998) itrans
      read (n1,*,err=998) icmpn
c                                 read new component definitions:
      do i = 1, itrans
         read (n1,'(a,1x,i2)') tcname(i), ictr(i)
         read (n1,*) (ctrans(j,i), j = 1, icmpn)
      
      end do

      read (n1,*,err=998) iwt
c                                 dummy variable place holders
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum  
      read (n1,*,err=998) idum
c                                 read code for choice of fluid equation
c                                 of state from terminal. 
      read (n1,*,err=998) ifug
      
      if (ifug.eq.8 .or.ifug.eq.10.or.ifug.eq.12.or.ifug.eq.16.or.
     *    ifug.eq.17.or.ifug.eq.19.or.ifug.eq.20.or.ifug.eq.24.or.
     *    ifug.eq.25) then

         read (n1,*,err=998) ibuf,hu,dlnfo2,elag

      else if (ifug.eq.6 .or.ifug.eq.7 .or.ifug.eq.11.or.ifug.eq.18.or.
     *         ifug.eq.21.or.ifug.eq.22.or.ifug.eq.23) then

        call error (77,0d0,0,' the input file specifies a disabled '//
     *                       'or ivalid internal fluid EoS')

      end if 

      if (ibuf.eq.5) read (n1,*,err=998) buf

      if (hu.eq.1) then 
c                                 hardwired fluid EoS endmember names
         eoscmp(1) = 'H2      '
         eoscmp(2) = 'O2      '

      else 

         eoscmp(1) = 'H2O     '
         eoscmp(2) = 'CO2     '

      end if 
c                                 no dependent variable
      iind = 0 
c                                 dummy variable
      read (n1,*,err=998) idum
c                                 idum is just a 1d/2d flag for 
c                                 gridded minimization, for backwards 
c                                 compatibility set the to 2d if > 2 or < 1.
      if (idum.eq.1) then 
         oned = .true.
      else
         oned = .false.
      end if 

      read (n1,*,err=998) idep
      read (n1,*,err=998) c0,c1,c2,c3,c4

      if (idep.eq.1) then 
         iind = 2
      else if (idep.eq.2) then 
         iind = 1
      end if 
c                                 decode thermodynamic components
c                                 read to the beginning of the component list
      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 
c                                 count (icp) and save names (cname)
      icp = 0
      jbulk = 0

      do 

         read (n1,'(a,a)') rname,strg

         if (rname.eq.'end t') then 
c                                 finished, check for no components
            if (icp.eq.0) then
               write (*,*) 'No thermodynamic components'
               goto 998
            else if (icopt.eq.5.and.jbulk.lt.icp) then 
               write (*,*) 'All thermodynamic components must be ',
     *                     'constrained.'
               goto 998
            end if 
        
            exit 

         else if (rname.eq.blank) then 
 
            cycle 

         else if (rname.eq.'Volume'.or.rname.eq.'Entropy') then

            usv = .true.

         else

            icp = icp + 1
            cname(icp) = rname
c                                 encode a graphics names for the
c                                 compositional variables, this is kind of
c                                 pointless, but it looks good.
            write (xname(icp),'(a,a,a)') 'x(',rname,')'
c                                 unblank the name
            call unblnk (xname(icp))
            if (icp.gt.k5) call error (197,r,icp,'INPUT1')

         end if 
c                                 check for compositional constraints
         read (strg,*,err=998) icont

         if (icopt.eq.12) then 
            k = 2
         else 
            k = icont
         end if 

         if (k.ne.0) then 
            jbulk = jbulk + 1
            read (strg,*,err=998) j, (dblk(i,jbulk), i = 1, k)
         end if 

      end do           

      icp1 = icp + 1
      icp2 = icp + 2

      if (usv) then

         hcp = icp2
         tindex = icp1
         pindex = icp2
         cname(tindex) = 'T(K) '
         cname(pindex) = '-P(b)'

      else

         hcp = icp

      end if 
c                                 decode saturated components
c                                 isat is the saturated component counter
      isat = 0
      io2  = 0 

      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 


      do 

         read (n1,'(a,a)') rname,strg
         if (rname.eq.blank) cycle 

         if (rname.eq.'end s') then 

            icomp = icp + isat
            exit 

         else if (rname.eq.blank) then 

            cycle 

         end if 
c                                 check for compositional constraints
         read (strg,*,err=998) jcont

         if (jcont.ne.0) then 
            jbulk = jbulk + 1
            read (strg,*,err=998) j, (dblk(i,jbulk), i = 1, jcont)
         end if

         isat = isat + 1
         if (isat.gt.h5) call error (15,r,i,'BUILD')
         cname(icp+isat) = rname
         if (rname.eq.'O2') io2 = isat

      end do 
c                                 decode saturated phase components
      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 
c                                 ifct is the saturated phase component counter
      ifct = 0

      do 

         read (n1,'(a)') rname

         if (rname.eq.'end s') then 
            icomp = icomp + ifct
            exit 
         else if (rname.eq.blank) then 
            cycle 
         end if 
      
         ifct = ifct + 1
         if (ifct.gt.2) call error (44,r,i,'BUILD')
c                                 save the component if only one
c                                 for use in input2.
         if (ifct.eq.1) zname = rname
         cname(icomp+ifct) = rname
      end do 
      
!!            write (*,1600) rname

1600  format (' thermodynamic data file  5', i2)
c                                  decode mobile components
c                                  jmct - mobile component counter
      jmct = 0 
      ifact = 0
      jmuct = 0 

      do 

         call rdstrg (n1,nstrg,string,eof)

         if (eof) then 

            goto 998

         else if (string(1).eq.'begin') then

            cycle 

         else if (string(1).eq.'end') then

            icomp = icomp + jmct
            exit 

         else 

            read (string(1),'(a5)') rname
            jmct = jmct + 1
            if (jmct.gt.2) call error (45,r,i,'BUILD')
            cname(icomp+jmct) = rname

            if (nstrg.eq.1) then 
c                                 old format, create variable name
               write (vname(3+jmct),'(a,a)') 'mu_',rname
               imaf(jmct) = 1
               jmuct = jmuct + 1

            else 
c                                 new format
               read (string(2),'(a1)') y
               vname(3+jmct) = string(2)
               afname(jmct) = string(3)

               if (y.eq.'m') then 
c                                 chemical potential
                  imaf(jmct) = 1
                  jmuct = jmuct + 1

               else if (y.eq.'f') then 

                  imaf(jmct) = 2

               else if (y.eq.'a') then 

                  imaf(jmct) = 3

               end if 

               if (imaf(jmct).gt.1) ifact = ifact + 1 

            end if 
               
         end if 

      end do 
c                             the ifct flag can probably be set later if fluid
c                             is in the thermodynamic composition space.   
      jfct = icp + isat 
c                             jprct+1..icomp -> (jmct.ne.0) mobile components 
      jprct = icomp - jmct 
c                             excluded phases
      ixct = 0
c                             decode excluded phases
      do 
         read (n1,'(a)',end=998) name
         if (name.eq.'begin ex') exit
      end do

      do 

        read (n1,'(a)') name

         if (name.eq.'end excl') then
            exit
         else if (name.eq.blank) then
            cycle 
         end if 

         ixct = ixct + 1
         if (ixct.gt.h8) call error (13,r,i,'BUILD')
         exname(ixct) = name

      end do  
c                             solution phases:
      do 
         read (n1,'(a)',end=998) sname
         if (sname.eq.'begin solu') exit
      end do
c                             isoct - solution phase counter,
c                             io9 is a flag = 0 no solution file
      isoct = 0

      do 

         read (n1,'(a)') sname
 
         if (sname.eq.'end soluti') then 
            if (io9.eq.1) isoct = 0 
            exit 
         else if (sname.eq.blank) then 
            cycle  
         end if 

         isoct = isoct + 1
         if (isoct.gt.h9) call error (25,r,i,'BUILD')
         fname(isoct) = sname

      end do  
c                             read the maximum pressure, temper-
c                             ature, xco2, u1, and u2; the minimum
c                             pressure temperature, xco2, u1, and u2;
c                             and the default pressure, temperature,
c                             xco2, and chemical
c                             potential increments use kelvins, bars and
c                             joules as units (if no mobile components
c                             enter two zeroes for each read).
      read (n1,*,err=998) vmax
      read (n1,*,err=998) vmin
      read (n1,*,err=998) dv
c                             read the default indices of the
c                             dependent, independent, and secondary
c                             independent intensive variables, p = 1,
c                             t = 2, and xco2 = 3, respectively.
      read (n1,*,err=998) (iv(i), i = 1, 5)
c                             check variable ranges are consistent,
c                             variable iv(1):
      if (icopt.ne.0.and.icopt.ne.4.and.iam.ne.2) then

         if (iv(1).eq.3.and.ifct.eq.0) call error (110,r,i,'I')

         if (iv(1).eq.3.and.ifct.eq.1) then 

            if (icopt.ne.7.and.iv(2).ne.3) call error (111,r,i,'I')

         end if 

         if (vmin(iv(1)).ge.vmax(iv(1)).and.icopt.lt.5) then 

            call error (112,r,i,'less than or equal')

         else if (vmin(iv(1)).eq.vmax(iv(1)).and.
     *            icopt.eq.5.and.icont.lt.3) then

            call error (112,r,i,'equal')

         end if 

         if (vname(iv(1)).eq.blank) call error (116,dip,i,'I')

      end if
c                             variable iv(2):
      if (iam.ne.2.and.(icopt.eq.1.or.
     *                  (icopt.eq.5.and.icont.eq.1.and..not.oned))) then

         if (iv(2).eq.3.and.ifct.eq.0) call error (110,r,i,'INPUT1')

         if (iv(2).eq.3.and.ifct.eq.1) call error (111,r,i,'INPUT1')

         if (icopt.eq.1) then 

            if (vmin(iv(2)).ge.vmax(iv(2))) call error (112,r,i,
     *                                            'less than or equal')

         else 

            if (vmin(iv(2)).eq.vmax(iv(2))) call error (112,r,i,'equal')

         end if 

         if (vname(iv(2)).eq.blank) call error (116,r,i,'INPUT1')

      end if
c                             if a chemical potential is specified as an
c                             independent variable (iv(1-3)), check if
c                             the variable is defined:
      kct = 0

      do i = 1, 3
         if (iv(i).gt.3) kct = kct + 1
      end do 
c                             identify the variable used to determine
c                             which phases lie on the left hand side
c                             of a reaction equation.
      if (icopt.eq.3) then
         ivfl = iv(1)
      else if (iv(1).eq.2.or.iv(2).eq.2) then
c                             choose T
         ivfl = 2
      else if (iv(1).eq.1.or.iv(2).eq.1) then
c                             no T, so choose P
         ivfl = 1
      else
c                             no P or T, choose independent V
         if (iv(2).ne.3) then
            ivfl = iv(2)
         else
            ivfl = iv(1)
         end if
      end if
c                             ok, now find out which variables are
c                             dummies and store the indexes of the
c                             non-dummy variables in jv.
      ipot = 0

      do i = 1, 5
c                             variables v(1) (p) and v(2) (t) are
c                             only dummies if idep is set.
         if ((iv(i).ne.idep.or.icopt.eq.7.or.icopt.eq.9).and.
     *       (iv(i).eq.1.or.iv(i).eq.2)) then
            ipot = ipot+1
            jv(ipot) = iv(i)
c                             variable v(3) is a dummy if ifct = 0:
         else if ((iv(i).eq.3).and.ifct.gt.0) then
            ipot = ipot+1
            jv(ipot) = iv(i)
c                             variables v(4) and v(4) are dummies if
c                             imyn = 1:
         else if (jmct.ne.0) then
            if (iv(i).eq.4) then
               ipot = ipot+1
               jv(ipot) = iv(i)
            else if (iv(i).eq.5.and.jmct.eq.2) then
               ipot = ipot+1
               jv(ipot) = iv(i)
            end if
         end if

      end do 
c                                 if dependent variable add to jv list, could
c                                 increment ipot, but maybe it's better not to.
      if (idep.ne.0) jv(ipot+1) = idep
c                                 set convergence criteria for routine univeq
      if (icopt.le.3) then 

         call concrt

      else if (icopt.eq.12) then 
c                                 0-d infiltration
         read (n1,*,err=998) iopt(36), nopt(36)

      end if 

      if (icopt.ne.0) close (n1)
c                                 open files requested in input
      call fopen (n2name,prt,n9name,err)
c                                 err only set for unsplt (iam.eq.14)
      if (err) return
c                                 read auxilliary input for 2d fractionation
      if (icopt.eq.9) call rdain
c                                 get runtime parameters
      if (first.or.(.not.first).and.(.not.output).or.iam.eq.13) 
     *   call redop1 (first,tfname)

      goto 999
c                                 archaic error trap
998   call mertxt (n2name,prject,'.dat',0)
      call error (27,r,i,n2name)

999   end

      subroutine input2 (first)
c----------------------------------------------------------------------
c input2 reads the thermodynamic data file for most perplex programs, 
c a (the?) notable exception being frendly that calls the parallel 
c routine jnput2.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*5 mnames(k16*k17)*8

      double precision twt(k5),tsel(k5),tcox(k5),cst
 
      integer i,j,im, ict, k, ifer,inames, jphct, imak(k16), iox
 
      logical eof, good, first

      integer iff,idss,ifug
      common / cst10 /iff(2),idss(h5),ifug

      double precision ctot
      common/ cst3  /ctot(k1)

      integer iwt
      common/ cst209 /iwt

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
     
      character*5 zname
      common/ cst209a /zname

      character*5 cname
      common/ csta4 /cname(k5)

      character*8 names
      common/ cst8 /names(k1)

      character*8 name
      common/ csta6 /name

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ic
      common/ cst42 /ic(k0)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision atwt
      common/ cst45 /atwt(k0) 

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer make
      common / cst335 /make(k10)

      integer eos
      common/ cst303 /eos(k10)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer iam
      common/ cst4 /iam

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec
c-----------------------------------------------------------------------
c                               initialization for each data set
c                               for k10 endmembers
      do i = 1, k10
         make(i) = 0 
         names(i) = ' '
      end do
c                               for k1 phases:
      do i = 1, k1
         ikp(i) = 0
      end do 
c                               other counters and flags:
      do i = 1, h5
         isct(i) = 0
      end do 
c                               counters for bounds
      iphct = 0
      lamin = 0 
      idsin = 0 
      idfl = 0
c                               read data base header, do component
c                               transformations, read make definitions.
      call topn2 (0)
c                               general input data for main program

c                               reorder thermodynamic components
c                               if the saturated phase components are 
c                               present
      if (lopt(7)) then

         do k = 1, ispec 
                             
            do i = 1, icp

               if (cname(i).eq.cmpnt(idspe(k))) then 

                  if (i.eq.k) exit 

                  cname(i) = cname(k)

                  do j = 1, 3
                     cst = dblk(j,i)
                     dblk(j,i) = dblk(j,k) 
                     dblk(j,k) = cst
                  end do 

                  cname(k) = cmpnt(idspe(k))

                  exit            

               end if 

            end do 

         end do 

      end if  
c                              load the old cbulk array
      if (ifct.gt.0) iphct = 2
c                               identify nonzero components.
c                               initialize icout(i) = 0
      do i = 1, icmpn
         icout(i) = 0
      end do

      do i = 1, icomp

         im = 0

         do j = 1, icmpn

            if (cname(i).eq.cmpnt(j)) then 

               twt(i) = atwt(j)
               tsel(i) = sel(j)
               tcox(i) = cox(j)

               ic(i) = j
               icout(j) = 1

               do k = 1, ispec
                  if (j.eq.idspe(k)) then 
                     iff(k) = i
                     idfl = idfl + 1
                  end if 
               end do 
 
               im = 1

            end if 

         end do 
c                               write error message if a component
c                               was not found:

         if (im.eq.0) then 
            write (*,1230) cname(i), (cmpnt(k), k = 1, icmpn)
            write (*,1240)
            stop
         end if 
 
      end do 
c                                 this segment is to check if
c                                 a possible saturated phase component
c                                 has been made a mobile component,
c                                 if there is also a saturated phase
c                                 component idfl is the identity of the
c                                 mobile component otherwise idfl = 0.
      if (ifct.eq.1.and.idfl.eq.2) then

         do i = 1, ispec
            if (zname.ne.cmpnt(idspe(i))) cycle 
            idfl = i
            exit 
         end do 

      else 
         idfl = 0
      end if
c                                 load atwts, sel in updated order
      do i = 1, icomp
         atwt(i) = twt(i)
         sel(i)  = tsel(i)
         cox(i)  = tcox(i)
         if (cox(i).lt.0d0) iox = i 
      end do 
c                                 convert weight to molar amounts
      if (jbulk.ne.0) then 

         if (iwt.eq.1) then 
            do i = 1, jbulk
               do j = 1, 3
                  dblk(j,i) = dblk(j,i)/atwt(i)
               end do 
            end do 
         end if 

         do i = 1, jbulk
            cblk(i) = dblk(1,i)
         end do   

      end if 
      
c                                 get composition vectors for entities
c                                 defined by a make definition:
      call makecp (inames,mnames,first)
c                                 loop to read reference phase data for
c                                 activity/fugacity variables
      ict = 0 

      if (ifact.gt.0) then
c                                 rewind and read 'til end of header
         call eohead (n2)

         good = .false.

         do

            call getphi (name,.false.,eof)

            if (eof) then 

               write (*,1000) (afname(i),i=1,jmct)
               write (*,1010)
               call errpau

            end if 
c                                 now look for a match with the 
c                                 reference phase names
            do i = 1, jmct

               if (name.eq.afname(i)) then 
c                                 got a match, count
                  iphct = iphct + 1

                  ict = ict + 1

                  idaf(i) = iphct
c                                 store thermodynamic parameters:
                  call loadit (iphct,.false.,.true.)
c                                 zero the component
c                 vnumu(i,iphct) = 0d0

                  if (imaf(i).eq.2) then 
c                                 if some cretin chooses fugacity, prevent
c                                 gphase from calling the EoS.   
                     eos(iphct) = ieos 

                  else if (lopt(7)) then 
c                                 check for special component names
c                                 this is necessary because loadit 
c                                 will not set isfp if ifct > 0.
                     do k = 1, ispec
                        if (name.ne.cmpnt(idspe(k))) cycle
                        eos(iphct) = 100 + k 
                        exit 
                     end do 
 
                  end if 
c                                 blank the name, this has two purposes,
c                                 it prevents problems if an entry is 
c                                 replicated in the data file, and flags
c                                 tagged entries 
                  afname(i) = ' '

                  if (ict.eq.ifact) good = .true.

                  exit 

               end if 

            end do 

            if (good) exit 

         end do 

      end if 
c                                 begin first read loop for data on
c                                 saturated components.
      if (isat.eq.0.and.ifct.eq.0) goto 40
c                                 read 'til end of header
      call eohead (n2)
c                                 loop to read real saturated
c                                 entities:
      ifer = 0

      do 

         call getphi (name,.false.,eof)

         if (eof) exit
 
         call chkphi (0,name,good)

         if (good) call sattst (ifer,good)

      end do 
c                                 loop to load made saturated entities
      do i = 1, nmak

         if (.not.mksat(i)) cycle
c                                 load make data 

         do j = 1, icmpn
            comp(j) = mcomp(i,j)
         end do 

         name = mknam(i,mknum(i)+1)
c                                 redundant check:
         call chkphi (2,name,good)
c                               
         if (.not.good) call error (57,comp(1),iphct,name)
c                                 set eos flag
         ieos = meos(i)

         call sattst (ifer,good)

         if (good) then 
            make(iphct) = i
c                                 pointer used for iemod.
            imak(i) = iphct
         end if 

      end do 
c                                 check that there is data for
c                                 every fluid component.
      if (ifct.gt.0.and.ifer.ne.ifct) call error (36,r,i,'INPUT2')
c                                 check that there is one phase
c                                 for each saturation constraint
40    do i = 1, isat
         if (isct(i).lt.1) call error (15,r,i,cname(icp+i))
      end do

      if (isat.gt.1.and.first.and.(iam.lt.4.or.iam.eq.15)) then

         write (*,'(a)') 'Summary of saturated-component entities:'

         do i = 1, isat

            write (*,1040) (cname(icp+j),j=1, i)
            write (*,1050) (names(ids(i,j)), j = 1, isct(i))

         end do

         if (iam.eq.15) write (*,'(a)') 
     *         '* solutions may also have composition'
     *       //'s consisting entirely of saturated components'

         write (*,'(/)')

      end if 
c                                 save endmembers that consist entirely 
c                                 of saturated phase or mobile components:
      kphct = iphct 

      if (ifct+jmct.gt.0) then 

         call eohead (n2)

         do 

            call getphi (name,.false.,eof)

            if (eof) exit

            call chkphi (4,name,good)

            if (.not.good) cycle 
c                                 reject phases already in the list
            do i = 1, kphct
               if (names(i).eq.name) then
                  good = .false.
                  exit
               end if 
            end do 

            if (.not.good) cycle             
c                                 matched a name
            iphct = iphct + 1
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)

         end do

      end if 
c                                 -------------------------------------
c                                 real entities in the thermodynamic 
c                                 composition space:
      istct = iphct + 1
c                                 read till end of header
      call eohead (n2)
c                                 loop to load normal thermodynamic data:
      do  
    
         call getphi (name,.false.,eof)

         if (eof) exit 
c                                 check if valid phase:
         call chkphi (1,name,good)

         if (good) then 
c                                 acceptable data, count the phase:
            iphct = iphct + 1
c                                 for normalized composition:
            ctot(iphct) = tot
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)
         end if 
      end do 
c                                 -------------------------------------
c                                 made entities (as opposed to the required
c                                 data read later):
      do i = 1, nmak

         if (mksat(i)) cycle
c                                 load make data 
         do j = 1, icmpn
            comp(j) = mcomp(i,j)
         end do 

         name = mknam(i,mknum(i)+1)
c                                 redundant check, but makes ctot.
         call chkphi (3,name,good)
c                               
         if (.not.good) call error (57,comp(1),iphct,name)

         iphct = iphct + 1
         ctot(iphct) = tot
c                                 set ieos flag to that of the first
c                                 real entity in the make definition
         ieos = meos(i)

         call loadit (iphct,.true.,.true.)

         make(iphct) = i
c                                 pointer used for iemod.
         imak(i) = iphct

      end do 
c                                 load thermodynamic data for make definitions and
c                                 solute species, at this point iphct points to the 
c                                 last real entity, save this value and restore it later.
      jphct = iphct
c                                 -------------------------------------
c                                 make definition data: this
c                                 data is saved in the arrays thermo
c                                 and cp by loadit, but are not counted,
c                                 i.e., the counters ipoint and iphct
c                                 are reset. soload will then load the
c                                 cp array over the values loaded here,
c                                 but thermo should not be affected. gmake
c                                 then gets the data using the array 
c                                 mkind. the names array will also be 
c                                 overwritten.
      call eohead (n2)

      do 

         call getphi (name,.true.,eof)

         if (eof) exit

         do i = 1, inames

            if (name.ne.mnames(i)) cycle
c                                 matched a name
            iphct = iphct + 1
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.false.)

         end do

      end do 

      do i = 1, nmak
c                                remake pointer array for makes 
         do j = 1, mknum(i)
            do k = jphct + 1, iphct
               if (names(k).ne.mknam(i,j)) cycle
               mkind(i,j) = k
            end do
         end do 
      end do  
c                                 -------------------------------------
c                                 aqueous species, thermo data, as is the
c                                 case for make data is loaded in thermo;
c                                 names and composition loaded into 
c                                 aqnam and aqcp.
      aqst = iphct 
c
      call eohead (n2)
c                                 loop to load solute data:
      do  
    
         call getphi (name,.true.,eof)

         if (eof) exit
c                                 skip non-solute standard state data
         if (ieos.ne.15.and.ieos.ne.16) cycle
c                                 check if valid species:
         call chkphi (1,name,good)
c                                 check for oxidation state of aqueous
c                                 data if aq_oxides is set:
c         if (good.and.lopt(36).and.oxchg) then

c            qchg = thermo(6,k10)

c            if (qchg.eq.0d0.and.comp(ic(iox)).ne.0d0.or.
c     *          qchg-cox(iox)*comp(ic(iox)).ne.0d0) then 

c               call warn (100,r,102,
c     *              name//' has been rejected; to retain '//name//
c     *              ' set aq_oxide_components to false.')

c               good = .false.

c            end if

c         end if 

         if (good) then 
c                                 acceptable data, count the phase:
            iphct = iphct + 1
c                                 for normalized composition, probably
c                                 con't need this, but could be used to
c                                 save molar wt or something like that:
            ctot(iphct) = tot
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)

         end if 

      end do
c                                write summary and checks
      if (aqct.gt.0) then 

         if (lopt(25).and.(ihy.eq.0.or.ioh.eq.0)) then 
            call warn (99,0d0,0,'missing H+ or OH- species, '//
     *                          'aq_output set = F (INPUT2)')
         end if 

         ichg = 0
         
         do i = 1, aqct 

            q(i) = thermo(6, aqst + i)
            q2(i) = q(i)**2

            if (q(i).ne.0d0) then 
               ichg = ichg + 1
               jchg(ichg) = i 
            end if

         end do 

         if (first.and.iam.lt.3) then 
            write (*,1020)
            do i = 1, aqct, 6
               k = i + 5
               if (k.gt.aqct) k = aqct
               write (*,1030) (aqnam(j),int(thermo(6,j+aqst)),j=i,k)
            end do 
            write (*,'(//)')
         end if

      else if (lopt(32).or.lopt(25)) then 

         if (first.and.iam.lt.4) 
     *       call warn (99,0d0,0,' no data for aqueous species, '
     *                 //'aq_output and aq_lagged_speciation disabled.')

         lopt(32) = .false.
         lopt(25) = .false.
  
      end if 
c                                reset ipoint counter, but do not 
c                                reset iphct, because the compositions
c                                of the make phases are necessary for
c                                chemical potential variables.
c                                really? then why was it reset here?
      iphct = jphct
      ipoint = jphct

      do i = 1, nmak
c                                make an iemod flag for made
c                                endmembers:   
         do j = 1, mknum(i)
            if (iemod(mkind(i,j)).eq.0) exit 
         end do 

         if (j.le.mknum(i)) cycle

         iemod(imak(i)) = iemod(mkind(i,1))

      end do 

1000  format ('**error ver007** at least one of the reference ',
     *        'endmembers:',/,5(a,1x))
1010  format ('needed to define an independent fugacity/activity ',
     *    'variable is missing',/,'most likely the endmember has ',
     *    'been rejected, if so then set',/,'the auto_exclude ',
     *    'option to FALSE.',/)
1020  format (/,'Summary of aqueous solute species:',//,
     *        6('name     chg   ')) 
1030  format (6(a,2x,i2,3x))
1040  format ('for ',15(a,1x))
1050  format (6(a,2x))
1230  format ('**error ver013** ',a,' is an incorrect component'
     *       ,' name, valid names are:',/,12(1x,a))
1240  format ('check for upper/lower case matches or extra blanks',/)

      close (n2)

      end

      subroutine setvr0 (i,j)
c--------------------------------------------------------------------
c setvr1 computes nodal variables for node ij, three cases:

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables

c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont
c----------------------------------------------------------------------

      if (icont.eq.1) then 

         v(iv1) = vmin(iv1) + (i-1)*dv(iv1)
         v(iv2) = vmin(iv2) + (j-1)*dv(iv2)
         call incdp0

      else if (icont.eq.2) then 

         v(iv1) = vmin(iv1) + (j-1)*dv(iv1)
         call incdep (iv1)

         cx(1) =  (i-1)*dvr(1)
         call setblk 

      else 

         cx(1) = (i-1) * dvr(1)
         cx(2) = (j-1) * dvr(2)
         call setblk

      end if 

      end

      subroutine setblk
c-----------------------------------------------------------------------
c for gridded minimization setblk computes the bulk composition
c and initializes the arrays for lpopt.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x0

      integer i,j

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
      x0 = 1d0

      if (lopt(1)) then 
c                                 closed composition
         do j = 1, icont-1
            x0 = x0 - cx(j)
         end do 

      end if 

      do j = 1, jbulk
         cblk(j) = x0*dblk(1,j)
      end do 
         
      do j = 1, jbulk
         do i = 2, icont 
            cblk(j) = cblk(j) + cx(i-1)*dblk(i,j)
         end do 
      end do
c                                 modify cblk here to change the 
c                                 composition before minimization.
      ctotal = 0d0 
c                                 get total moles to compute mole fractions             
      do i = 1, hcp
         ctotal = ctotal + cblk(i)
      end do

      do i = 1, hcp 
         b(i) = cblk(i)/ctotal
      end do

      end 

      subroutine setvar 
c--------------------------------------------------------------------
c setvar initializes the variables for gridded minimization, three
c cases:

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j

      double precision rloopy,rloopx

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc 

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      integer iam
      common/ cst4 /iam

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
      if (iam.eq.3) then 
c                                 WERAMI (3), PSSECT (7):
c                                 jinc will only be ~= 1 only for 
c                                 2d intermediate grid results
         rloopy = dfloat((loopy-1)/jinc)
         rloopx = dfloat((loopx-1)/jinc)

      else

         rloopy = dfloat(loopy-1)
         rloopx = dfloat(loopx-1)

      end if 
c                                 for 1d calculations
      if (loopx.eq.1.or.loopx.eq.0) rloopx = rloopy

      do i = 1, ipot
         v(jv(i)) = vmin(jv(i))
      end do

      call incdp0

      if (icopt.eq.7.and.fileio) then 
c                                using nodal coordinate system
         dvr(1) = 1d0

      else if (icopt.eq.9.or.icopt.eq.11) then 
c                                using non-thermodynamic coordinate frame
         dvr(1) = (vmx(1) - vmn(1))/rloopx
         dvr(2) = (vmx(2) - vmn(2))/rloopy

      else if (icopt.eq.12) then 

         dvr(1) = nopt(36)
         dvr(2) = 1
         loopx = iopt(36)
         rloopx = dfloat(loopx)

      else if (icont.eq.1) then 
c                                v(iv1) on x, v(iv2) on y
         dv(iv1) = (vmax(iv1) - vmin(iv1))/rloopx
         dvr(1) = dv(iv1)

         dv(iv2) = (vmax(iv2) - vmin(iv2))/rloopy
         dvr(2) = dv(iv2)

      else if (icont.eq.2) then 
c                               composition is on x, v(iv1) on y
         dvr(1) = 1d0/rloopx

         dv(iv1) = (vmax(iv1) - vmin(iv1))/rloopy
         dvr(2) = dv(iv1)

      else 
c                                compositions on both axes
         dvr(1) = 1d0/rloopx
         dvr(2) = 1d0/rloopy 
         cx(1) = 0d0
         cx(2) = 0d0

      end if 
c                                set the bulk composition:
      do j = 1, jbulk
         if (icont.ne.0) then 
            cblk(j) = dblk(1,j)
         else 
            cblk(j) = 1d0
         end if 
      end do 

      end 

      subroutine inipot 
c--------------------------------------------------------------------
c setvar initializes the independent potential variables to their 
c minimum values
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------
c                                 initialize potentials
      do i = 1, ipot
         v(jv(i)) = vmin(jv(i))
      end do 
c                                 set dependent potential, if it exists
      call incdp0

      end

      subroutine getcmp (jd,id,ids)
c-----------------------------------------------------------------------
c getcmp gets the composition of pseudocompund id, where:
c  if ids < 0, -ids points to the composition of a true compound in array cp
c  if ids > 0, id points to the composition of a solution defined in terms
c              on endmember fractions defined and saved by routine resub
c              in array zcoor.
c the composition is saved in arrays cp3 and x3, entry jd

c getcmp is called by both WERAMI and MEEMUM/VERTEX

c this is an attempt to use the compositions stored in cp2....
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, id, jd, ids

      logical bad

      double precision xx
c                                 -------------------------------------
c                                 global variables:
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 single site solution coordinates:
      integer jend
      common/ cxt23 /jend(h9,m4)
c                                 refined compositions and solution 
c                                 pointer
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision cp0
      common/ cst71 /cp0(k0,k5)

      integer iam
      common/ cst4 /iam

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd
c----------------------------------------------------------------------
      kkp(jd) = ids
      cptot(jd) = 0d0

      if (ids.lt.0) then
c                                 simple compounds
         if (iam.ne.5) then
c                                 all programs except frendly 
            do i = 1, icomp
               cp3(i,jd) = cp(i,-ids)
            end do 
c                                 check if it's a solution endmember
            if (ikp(-ids).ne.0) call endcp (jd,-ids,ikp(-ids))
   
         else 
c                                 frendly 
            do i = 1, k0
               cp3(i,jd) = cp0(i,-ids)
            end do 

         end if

      else 
c                                 solutions, initialize
         do i = 1, icomp
            cp3(i,jd) = 0d0
         end do
c                                 if getcmp is being called by WERAMI:
c                                 GETXZ (dlib.f) gets the x(i,j) coordinates for the
c                                 composition from the x3(jd,i,j) array and
c                                 the id argument is irrelevamt. 
c                                 if getcmp is being called by MEEMUM/VERTEX:
c                                 GETXZ (getxz1.f) gets both the x(i,j) and 
c                                 x3(jd,i,j) compositional coordinates from the
c                                 zcoor array.
         call getxz (jd,id,ids)
c                                 convert the x(i,j) coordinates to the
c                                 geometric y coordinates
         call xtoy (ids,jd,.true.,bad)

         if (lopt(32).and.ksmod(ids).eq.39) then 

            if (iam.ne.3) then
c                                 MEEMUM:
c                                 cp2 works for meemum/vertex, but not werami
c                                 the id index on cp2 is intentional.
               do j = 1, icomp 
                  cp3(j,jd) = cp2(j,id)*c2tot(id)
               end do

               if (cp2(k5,id).gt.1d1) then
                  write (*,*) 'BAZORK'
                  write (*,*) 'BAZORK'
                  write (*,*) 'hydroxyl solution stable ',cp2(k5,id)
                  write (*,*) 'BAZORK'
                  write (*,*) 'BAZORK'
               end if 

            else
c                                  WERAMI:
               if (caq(jd,na1).eq.0d0) then
c                                  pure solvent, use the y array to be safe
                  do i = 1, ns
                     do j = 1, icomp 
                        cp3(j,jd) = cp3(j,jd) + y(i) * cp(j,jnd(i))
                     end do 
                  end do

               else 
c                                  impure solvent
                  do i = 1, ns
                     do j = 1, icomp 
                        cp3(j,jd) = cp3(j,jd) + caq(jd,i) * cp(j,jnd(i))
                     end do 
                  end do

                  do i = sn1, nsa

                     k = i - ns
c                                 convert molality to mole fraction (xx)
                     xx = caq(jd,i)/caq(jd,na2)

                     do j = 1, icomp
                        cp3(j,jd) = cp3(j,jd) + xx * aqcp(j,k)
                     end do  

                  end do

               end if

            end if

         else if (lrecip(ids)) then
c                                 get the p' coordinates (amounts of 
c                                 the independent endmembers)     
            call getpp (ids) 

            do i = 1, lstot(ids)
               do j = 1, icomp 
                  cp3(j,jd) = cp3(j,jd) + p0a(i) * cp(j,jend(ids,2+i))
               end do 
            end do          

         else if (ksmod(ids).eq.20) then 
c                                 electrolyte:
c                                 solute species  
            do i = sn1, nqs
               do j = 1, icomp
                  cp3(j,jd) = cp3(j,jd) + y(i) * aqcp(j,jnd(i) - aqst)
               end do
            end do 
c                                 solvent species 
            do i = 1, ns 
               do j = 1, icomp
                  cp3(j,jd) = cp3(j,jd) + y(i) * cp(j,jnd(i))
               end do
            end do

         else 
c                                 solutions with no dependent endmembers:
c                                 y coordinates used to compute the composition
            do i = 1, mstot(ids)
               do j = 1, icomp
                  cp3(j,jd) = cp3(j,jd) + y(i) * cp(j,jend(ids,2+i))
               end do
            end do

         end if 

      end if 

      do i = 1, icp
         cptot(jd) = cptot(jd) + cp3(i,jd)
      end do 

      end 

      subroutine inblnk (text,char)
c----------------------------------------------------------------------
c inblnk - scan text to last '/' or '\' and insert char after.
 
c     text - character string 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar
 
      character text*(*), bitsy(lchar)*1, char*1 
c----------------------------------------------------------------------
      nchar = len(text) 
      read (text,1000) (bitsy(i), i = 1, nchar)
c                                 scan for blanks:

      do i = nchar,1,-1
c                                 this line may cause problems
c                                 on some operating systems that 
c                                 recognize the backslash as an escape
c                                 character.
         if (bitsy(i).eq.'/') goto 10
         bitsy(i+1) = bitsy(i)
      end do 

      i = 0

10    bitsy(i+1) = char

      write (text,1000) (bitsy(i), i = 1, nchar)
 
1000  format (400a)

      end

      subroutine matchj (unnown,itis)
c----------------------------------------------------------------------
 
c matchj - subroutine to determine if the string unnown is a valid
c          solution or compound name.
 
c   itis = -id if compound
c   itis = ikp if solution 
c   itis = 0 if invalid
c----------------------------------------------------------------------
      implicit none

      integer i, itis
 
      character*10 unnown
 
      include 'perplex_parameters.h'
 
      integer isoct
      common/ cst79 /isoct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character names*8
      common/ cst8 /names(k1)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c---------------------------------------------------------------------- 
 
      itis = 0

      do i = 1, isoct
         if (unnown.eq.fname(i)) then
             itis = i
             goto 99
         end if
      end do

      do i = 1, iphct
         if (unnown.eq.names(i)) then
            itis = -i
            goto 99
         end if
      end do 

99    end

      subroutine maktit 
c-----------------------------------------------------------------------
c create a title for graphics output, the title consists of the 
c calculation title + saturation hierarchy (provided one is 
c specified) and is the first two elements of title (csta8).
c if icopt = 1 or 3, also adds a blurb about reaction convention.

c title is max 3 lines, but four lines are written to be consistent
c with old plot file formats written by frendly, pt2curv etc.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      character*162 title
      common/ csta8 /title(4)

      character*8 vname,xname     
      common/ csta2  /xname(k5),vname(l2)

      integer ivfl
      common/ cst102 /ivfl

      character*5 cname
      common/ csta4 /cname(k5)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
      do i = 2, 4
         title(i) = ' '
      end do                              
c                               saturated and buffered component names:
      if (isat.gt.0) then 
         write (title(2),1070) (cname(i+icp), i= 1, isat)
      else 
         write (title(2),1000) ' '
      end if 
c                                 reaction convention
      if (icopt.eq.1.or.icopt.eq.3) write (title(3),1080) vname(ivfl)

      do i = 1, 3
         call deblnk (title(i))
      end do 

1000  format (a)
1070  format ('Component saturation hierarchy: ',7(a,1x))
1080  format ('Reaction equations are written with the high ',
     *         a,'assemblage to the right of the = sign')

      end

      subroutine rdain
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for 2d fractionation 
c calculations, called by VERTEX and WERAMI
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer maxbox,lay, mpol, mord
      parameter (maxbox=1760,lay=6,mpol=3,mord=4) 

      logical dynam, titrat, qfile 

      integer i,j,k,ier

      double precision zlayer

      character*100 name

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc 

      logical pzfunc
      integer ilay,irep,npoly,nord
      double precision abc0,vz,iblk
      common/ cst66 /abc0(0:mord,mpol),vz(6),iblk(lay,k5),ilay,
     *               irep(lay),npoly,nord,pzfunc

      double precision a,b
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23  /a(k8,k8),b(k8),ipvt(k8),idv(k8),iophi,idphi,
     *                iiphi,iflg1

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      character*100 cfname
      common/ cst227 /cfname

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c-----------------------------------------------------------------------
c                                 look for input data from a file 
c                                 of type aux
      call mertxt (name,prject,'.aux',0)
      open (n8,file=name,status='old',iostat=ier)

      if (ier.ne.0) call error (51,vz(1),ilay,name)

      call mertxt (name,prject,'.fld',0)
      open (n12,file=name)
c                                 set the number of independent variables
c                                 to 1 (the independent path variable must
c                                 be variable jv(1), and the dependent path
c                                 variable must be jv(2), the path variables
c                                 can only be pressure and temperature
      ipot = 1

      jbulk = icp
c                                 true => flush model, ~true => subducting column
      read (n8,*) dynam
c                                 this is just a trick to avoid changing the variable 
c                                 names, the i/o needs to be rewired for dynam + titrat
      flsh = .not.dynam
c                                 true => basal mass flux
      read (n8,*) titrat
c                                 true => anneal the column
      read (n8,*) anneal
c                                 true => don't output console info on the fractionated phase
      read (n8,*) short
c                                 true => p-t field from file/internal function
      read (n8,*) pzfunc
c                                 Perple_X assumes upward directed depth, but to 
c                                 make the input intuitive, the input is specified
c                                 in downward coordinates, hence the sign changes 
c                                 below:
c                                 thickness of a box in column
      read (n8,*) vz(1)
c                                 gradient in variable jv(1) with z, jv(1)
c                                 is the independent variable, for subduction
c                                 this is logically pressure, i.e., dp(bar)/dz(m)
      read (n8,*) vz(2)
c                                 z if flush or or dzmax if not flush
      if (flsh) read (n8,*) vz(3)
c                                 value of the x-coordinate at the origin
      vz(4) = 0d0
      if (.not.flsh) read (n8,*) vz(4)
c                                 max value of the x-coordinate
      read (n8,*) vz(5)
c                                 Perple_X assumes an upward directed column coordinate
c                                 but in frac2d 
c                                 make the input intuitive, the input is specified
c                                 in downward coordinates, hence the sign changes 
c                                 below:

      if (flsh) then
c                                 specification n t-z points to fit n-1^th order
c                                 polynomial 
         read (n8,*) npoly

         if (npoly.gt.mpol) call error (77,b(1),i,'too many t-z '/
     *                     /'coordinates increase mpol in common cst66')

         do i = 1, npoly

            read (n8,*) b(i), a(i,1)

            do j = 2, npoly - 1
               a(i,j) = a(i,1)**j
            end do

            a(i,j) = 1d0

         end do

         call factor (a,npoly,ipvt,i)

         if (i.eq.0) call subst (a,ipvt,npoly,b,i)

         if (i.ne.0) call error (77,b(1),i,'degenerate t-z'//
     *                                     ' coordinates, FRAC2D')
         do i = 1, npoly
            abc0(1,i) = b(i)
         end do

      else
c                                 now we need a path function for the dependent
c                                 variable, here we take a function defined in
c                                 terms of the absolute depth of the top of the
c                                 column (z0) and the relative depth (dz) within
c                                 the column
         if (.not.pzfunc) then
c                                 slab dip (degree)
            read (n8,*) vz(6)
c                                 number of geothermal polynomials
            read (n8,*) npoly
            if (npoly.gt.mpol) call error (77,b(1),i,'too many '/
     *                       /'geotherms increase mpol in common cst66')
c                                 order of geothermal polynomials
            read (n8,*) nord
            if (nord.gt.mord) call error (77,b(1),i,'geothermal '/
     *      /'polynomial order too high, increase mord in common cst66')

            do i = 1, npoly
c                                 depth in column for the i'th geotherm
              read (n8,*) abc0(nord+1,i) 
c                                 convert orthogonal depth to vertical depth
              abc0(nord+1,i) = abc0(nord+1,i) / 
     *                         dcos(vz(6)*.1745329252d-1)
c                                 polynomial coefficients for the geotherm
              read (n8,*) (abc0(j,i), j = 0, nord)

            end do

         end if

      end if 
c                                 get the initial global composition array
c                                 consisting of ibox compositions defined 
c                                 in terms of icp components. this read
c                                 statement assumes that H2O an CO2 (if 
c                                 thermodynamic components) are the 1st and
c                                 2nd components (if present). 
      ilay = 0
      ncol = 0
c                                 number of nodes with appended composition
c                                 end of data indicated by zero 
      do 

         read (n8,*) zlayer

         if (zlayer.eq.0) exit 

         ilay = ilay + 1

         if (ilay.eq.lay) call error (77,b(1),i, 
     *                               'increase lay in common cst66')

         read (n8,*) (iblk(ilay,i),i=1,icp)

         irep(ilay) = idint(zlayer/vz(1)) 

         ncol = ncol + irep(ilay)

         if (ncol.gt.maxbox) call error (77,b(1),i, 
     *                            'increase maxbox in common cst66')

      end do
c                                 read aliquot composition
      if (titrat) then 

         if (ilay+1.eq.lay) call error (77,b(1),i, 
     *                               'increase lay in common cst66')

         read (n8,*) qfile

         if (qfile) then
            write (*,*) 'oink'
            call errpau
         else 
            read (n8,*) (iblk(ilay+1,i),i=1,icp)
         end if 
      end if 

      close (n8)
c                                 two cases, file input or analytical
      if (fileio) then 
c                                 file input of nodal p-t coordinates
         open (n8,file=cfname,status='old',iostat=ier)
c                                 read header info
         read (n8,*) i, nrow

         if (ncol*nrow.gt.k2) call error (77,b(1),i,'too many'/
     *      /' coordinates, increase k2 to ncol*nrow in routine FRAC2D')

         if (i.ne.ncol) call error (77,b(1),i,'the number of'//
     *     'nodes in a column specified in: '//cfname//'must equal the'/
     *    /' number of nodes specified in the aux file.')

         do i = 1, nrow

            k = (i-1) * ncol

            do j = 1, ncol
               read (n8,*) vn(k+j,1),vn(k+j,2)
            end do 

         end do

         close (n8)

      end if

      end

      subroutine fr2dpt (p0,dz)
c----------------------------------------------------------------------
c subroutine to set p-t variables from i-j coordinates in 2d-fractionation
c calculations, called by VERTEX and WERAMI
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lay,i,j,mpol,mord

      parameter (lay=6,mpol=3,mord=4) 

      double precision p0, z0, dz, z2, z3, z4, z5, z6, t0, t1, t2,aa,bb

      logical pzfunc
      integer ilay,irep,npoly,nord
      double precision abc0,vz,iblk
      common/ cst66 /abc0(0:mord,mpol),vz(6),iblk(lay,k5),ilay,
     *               irep(lay),npoly,nord,pzfunc

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      double precision a,b
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23  /a(k8,k8),b(k8),ipvt(k8),idv(k8),iophi,idphi,
     *                iiphi,iflg1


      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
      if (fileio) then 
c                                convert p0-dz coordinate to nodal 
c                                values
         i = idint((p0 - vmin(1))/dv(1)) + 1
         j = ncol - idint(dz/vz(1))

         v(1) = vn((i-1)*ncol + j, 1)
         v(2) = vn((i-1)*ncol + j, 2)

      else if (pzfunc) then
c                                 this could be made a lot more efficient by
c                                 making the quadratic coeffs once for each column
         z0 = p0/1d3
         z2 = z0*z0
         z3 = z2*z0
         z4 = z3*z0
         z5 = z4*z0
         z6 = z5*z0

         t2 = -0.1099312D-6*z4 +0.5065153D-4*z3 -0.3902580D-2*z2 
     *        +0.3024415D0 *z0 +0.8107985D3

          if (z0.lt.75d0) then
c                                t0, t1 shallow
             t0 = 0.1255734D-5*z5 -0.2000554D-3*z4 +0.1180485D-1*z3 
     *           -0.3163565D0 *z2 +0.6026698D1 *z0 + 0.276185544D3
             t1 = 0.1409099D-4*z4 -0.1603057D-2*z3 + 0.5553760D-1*z2 
     *           +0.2762566D0 *z0 +0.4401928241D3
          else if (z0.lt.78.99d0) then
c                                t0 deep
             t0 = -0.2059655D-9*z6 +0.2323113D-6*z5 - 0.1076535D-3*z4 
     *            +0.2625959D-1*z3 -0.3566382D1 *z2 + 0.2582593D3 *z0 
     *            -0.6916326D4
c                                t1 shallow
             t1 = 0.1409099D-4*z4 -0.1603057D-2*z3 + 0.5553760D-1*z2 
     *           +0.2762566D0 *z0 +0.4401928241D3
          else
c                                t0, t1 deep
             t0 = -0.2059655D-9*z6 +0.2323113D-6*z5 - 0.1076535D-3*z4 
     *            +0.2625959D-1*z3 -0.3566382D1 *z2 + 0.2582593D3 *z0 
     *            -0.6916326D4

             t1 = -0.3998088D-6*z4 +0.3672092D-3*z3 - 0.1290587D0*z2 
     *            +0.2181334D2 *z0 -0.5161647D3
          end if

         aa = -t1 / 272d0 + t2 / 850d0 + t0 / 400d0
         bb = -dsqrt(2d0) * (64d0*t2 - 625d0*t1 + 561d0*t0)/6800d0

         v(1) = (p0 - dz) * vz(2)

         v(2) = aa*dz**2/1d6 - bb*dz/1d3 + t0

      else if (flsh) then 

         z0 = (vz(3) -dz)
         v(1) = z0 * vz(2)
         v(2) = abc0(1,npoly)

         do i = 1, npoly-1
            v(2) = v(2) + abc0(1,i) * z0 ** i
         end do

      else
c                                 compute the npoly t-corrdinates
         do i = 1, npoly 
c                                 b - geotherm t
            b(i) = abc0(0,i)
c                                 depth for geotherm
            z0 = p0 + abc0(nord+1,i)

            do j = 1, nord
               b(i) = b(i) + abc0(j,i) * z0**j
            end do

            do j = 1, npoly-1
               a(i,j) = z0**j
            end do 

            a(i,j) = 1d0

         end do

         call factor (a,npoly,ipvt,i)

         if (i.eq.0) call subst (a,ipvt,npoly,b,i)

         if (i.ne.0) call error (77,b(1),i,'degenerate t-z'//
     *                                     ' coordinates, FRAC2D')
c                                  true depth is 
         z0 = p0 - dz
c                                  pressure is
         v(1) = z0 * vz(2)
c                                  temperature is
         v(2) = b(npoly)

         do i = 1, npoly-1
            v(2) = v(2) + b(i) * z0** i
         end do

      end if

      end

      subroutine getpp (id)
c-----------------------------------------------------------------------
c getpp computes the amounts of the indepdendent edmembers of a reciprocal
c solution in terms of the disordered endmembers (i.e., the p coordinates
c corrected for the amounts of the ordered species if present [ksmod=8]).
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
c                                  first convert the istot disordered
c                                  endmember coordinates to the 
c                                  kstot + nord p0 coordinates
      call y2p0 (id) 
c                                  decompose ordered species
      if (nord(id).gt.0) call p0dord (id)

      end

      subroutine fopen (n2name,prt,n9name,err)
c-----------------------------------------------------------------------
c open files for subroutine input1.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, err, tic 

      integer ier

      character n2name*100, prt*3, name*100, n9name*100

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer iam
      common/ cst4 /iam

      save first

      data first/.true./
c----------------------------------------------------------------------
c                                 open thermodynamic data file
      call fopen2 (0,n2name)

      tic = .false.
      err = .false.

      if (iam.eq.3.or.iam.eq.7.or.iam.eq.14) then
c                                 use existing plt/blk files
c                                 iam - 14 - unsplt (local)

c                                 plt/blk files for werami/pssect opened 
c                                 later by redplt to allow interim results
         if (iam.eq.14) then 
c                                 open the plot file
            call mertxt (name,prject,'.plt',0)

            open (n4, file = name, iostat = ier, status = 'old')

            if (ier.ne.0) err = .true.
c                                 open assemblage file
            call mertxt (name,prject,'.blk',0)

            open (n5, file = name, iostat = ier, status = 'old')

            if (ier.ne.0) err = .true.

         end if

      else if (iam.eq.1.or.iam.eq.2.or.iam.eq.13.or.iam.eq.15) then 
c                                 iam -  1 - vertex
c                                 iam -  2 - meemum
c                                 iam - 13 - unsplt (global)
c                                 iam - 15 - convex

         if (first) then
            tic = .true.
            call mertxt (name,prject,'.dat',0)
!!            write (*,1160) name
!!            write (*,1170) n2name
         end if



!!        write(*,*) ' the name of the input file is in clib'

c                                 open print/plot files if requested
         if (prt.ne.' '.and.prt.ne.'no_'.and.iam.ne.13) then 

            io3 = 0 
            call mertxt (name,prject,'.prn',0)
            open (n3, file = name)

         else

            io3 = 1
            name = 'none requested'

         end if

         if (first.and.iam.ne.2) then
c                                 plt output file
            io4 = 0
            call mertxt (name,prject,'.plt',0)
            if (iam.ne.13) write (*,1180) name

            open (n4, file = name, iostat = ier, status = 'new')
            if (ier.ne.0) then 
               open (n4, file = name)
               close (n4, status = 'delete')
               open (n4, file = name)
            end if

            write (*,1190) name

            if (iam.ne.15) then 
c                                 blk output file
               call mertxt (name,prject,'.blk',0)
               open (n5, file = name, iostat = ier, status = 'new')
               if (ier.ne.0) then 
                  open (n5, file = name)
                  close (n5, status = 'delete')
                  open (n5, file = name)
               end if

               write (*,1220) name

            end if

         else if (iam.ne.15) then 

            rewind (n5)

         end if

      else

         call error (999,0d0,n9,'oops fopen')

      end if 

      if (n9name.ne.' ') then

         io9 = 0 
c                                 open solution model file
         open (n9,file = n9name,iostat = ier,status = 'old')
         if (ier.ne.0) call error (120,0d0,n9,n9name)

!!         if (tic) write (*,1210) n9name

      else

         io9 = 1
         if (tic) write (*,1210) 'not requested'

      end if

      first = .false.

1160  format (/,'Reading problem definition from file: ',a)
1170  format ('Reading thermodynamic data from file: ',a)
1180  format ('Writing print output to file: ',a)
1190  format ('Writing plot output to file: ',a)
1210  format ('Reading solution models from file: ',a)
1220  format ('Writing phase assemblage data to file: ',a)

      end

      subroutine outgrd (loopx,loopy,jinc,lun,ind2)
c----------------------------------------------------------------------
c output grid data to the plot file, called by VERTEX and UNSPLT
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer loopx,loopy,jinc,i,j,jst,kst,kd,ltic,iend,lun,jgrd(l7),
     *        ind1, ind2, ier

      logical ext 

      character string*(lchar), name*170, text*3

      integer igrd
      common/ cst311 /igrd(l7,l7)

      integer iam
      common/ cst4 /iam

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      character tname*10
      logical refine, resub
      common/ cxt26 /refine,resub,tname
c----------------------------------------------------------------------
      if (lun.ne.n4) then 
c                                 write interim result file list
         call mertxt (name,prject,'.irf',0)
         inquire (lun, opened = ext)

         open (lun, file = name, iostat = ier, position = 'append')

         if (ier.eq.0) then 

            if (refine) then 
               ind1 = 1
            else 
               ind1 = 0
            end if 

            write (lun,*) ind1, ind2

            close (lun)

c                                 writing interim blk file
            write (text,'(a,i1,i1)') '_',ind1, ind2
            call mertxt (name,prject,text,0)
            call mertxt (name,name,'.blk',0)

            open (lun, file = name)

            rewind (n5)
c                                 the length of text should be able to 
c                                 handle format 1010 in outbl1
            do

               read (n5,'(a)',end=99) name
               write (lun,'(a)') name

            end do

99          close (lun)
            backspace (n5)
c                                 and the interim plt file
            call mertxt (name,prject,text,0)
            call mertxt (name,name,'.plt',0)
            open (lun, file = name)

         else 

            write (*,1000)
            return

         end if

      end if

      if (iam.eq.1) then
         write (lun,*) loopx, loopy, jinc
      else 
c                                 unsplt, set jinc to flag for werami
         write (lun,*) loopx, loopy, -1
      end if
c                                 fill in grid
      do i = 1, loopx, jinc

         if (i.ne.1.and.igrd(i,1).eq.0) then 
            jgrd(1) = igrd(i-jinc,1)
         else 
            jgrd(1) = igrd(i,1)
         end if 

         kst = 1

20       jst = kst

         if (i.ne.1.and.igrd(i,jst).eq.0) then
            jgrd(jst) = igrd(i-jinc,jst)
         else 
            jgrd(jst) = igrd(i,jst)
         end if 

         kd = jgrd(jst)

         ltic = -1

         do j = jst, loopy

            if (i.ne.1.and.igrd(i,j).eq.0) then 
               jgrd(j) = igrd(i-jinc,j)
            else 
               jgrd(j) = igrd(i,j)
            end if 

            if (jgrd(j).eq.0.or.jgrd(j).eq.kd) then
               ltic = ltic + 1
               if (j.eq.loopy) write (lun,*) ltic,kd
            else 
               write (lun,*) ltic,kd
               kst = j
               goto 20
            end if 
         end do 
      end do
c                                 write assemblage list
      write (lun,*) iasct

      do i = 1, iasct
         write (lun,*) iavar(1,i),iavar(2,i),iavar(3,i)
         write (lun,*) (idasls(j,i), j = 1, iavar(3,i))
      end do 

      if (lun.ne.n4) then
c                                 close interim plt file
         close (lun)

      else if (io3.eq.0) then 
c                                 write assemblages to print file
         write (n3,'(/,1x,a,a,/)') 'Stable assemblages identified ',
     *                         'by assemblage index:'
         do i = 1, iasct
            call psbtxt (i,string,iend)
            write (n3,'(i4,a,400a)') i,' - ',(chars(j), j = 1, length)
         end do

      end if

1000  format (/,'**warning ver999** OS fileio error, the irf file is'
     *        ,' corrupt and interim results',/,'for this calculation'
     *        ,' will be unreadable unless the irf file is edited',/)

      end

      subroutine psbtxt (id,string,iend)
c----------------------------------------------------------------------
c subprogram to write a text labels for bulk composition output 
c id identifies the assemblage

      implicit none

      include 'perplex_parameters.h'

      character string*(*), pname*14

      integer i, j, ist, iend, id, np, ntot, ids

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)
c----------------------------------------------------------------------
      iend = 0

      string = ' '

      ist = 1
      np = iavar(1,id)
      ntot = iavar(3,id)

      do i = 1, lchar
         chars(i) = ' '
      end do
c                                 first solution names:
      do i = 1, ntot
             
         ids = idasls(i,id)

         call getnam (pname,ids) 

         ist = iend + 1
         iend = ist + 14
         read (pname,'(400a)') (chars(j),j=ist,iend)

         call ftext (ist,iend)

      end do 

      write (string,'(400a)') (chars(j),j=1,iend) 

      length = iend

      end 

      subroutine redplt (name,err)
c-----------------------------------------------------------------------
c open/read plt/blk files for PSSECT and WERAMI.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character name*100

      integer ier

      logical err

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
      err = .false.
c                                 open the plot file
      call mertxt (tfname,name,'.plt',0)
      open (n4, file = tfname, iostat = ier, status = 'old')
      if (ier.ne.0) then 
         err = .true.
         return
      end if 
c                                 open assemblage file
      call mertxt (tfname,name,'.blk',0)
      open (n5, file = tfname, iostat = ier, status = 'old')
      if (ier.ne.0) then 
         err = .true.
         return
      end if 
c                                 read grid data:
      call plinp (err)
      if (err) return
c                                 read bulk composition data:
      call bplinp (err)

      end

      subroutine interm (finish,err)
c-----------------------------------------------------------------------
c if finish (only vertex) close plt/blk and delete interim results else 
c if ~finish open/read plt/blk files for PSSECT, UNSPLT, and WERAMI.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character yes*1, text*3, name*100

      integer ier, jnd(12,2), i, j, ind1, ind2

      logical err, finish, inter

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      character tname*10
      logical refine, resub
      common/ cxt26 /refine,resub,tname

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
      if (finish) then

         close (n4)
         close (n5)

         if (iopt(34).eq.1) then 
c                                 delete interim results
            call mertxt (tfname,prject,'.irf',0)
            open (1000, file = tfname, status = 'old', iostat = ier)
            if (ier.ne.0) return

            do

               read (1000,*,iostat=ier) i,j
c                                 file is in use or end of irf file
               if (ier.ne.0) exit 
c                                 make the root
               write (text,'(a,i1,i1)') '_',i,j
               call mertxt (name,prject,text,0)

               call mertxt (tfname,name,'.plt',0)
               open (1001, file = tfname, status = 'old', iostat = ier)
               if (ier.ne.0) exit 
               close (1001, status = 'delete')

               call mertxt (tfname,name,'.blk',0)
               open (1001, file = tfname, status = 'old', iostat = ier)
               if (ier.ne.0) exit 
               close (1001, status = 'delete')

            end do

            close (1000, status = 'delete')

         end if 

         return

      end if 

      if (iopt(34).ne.2.or.icopt.ne.5.or.iam.eq.14) then 
c                                 for all calculations other than 2d gridded 
c                                 min OR if interim_results (iopt(34)) < 2
c                                 try to open final plt and blk files
         name = prject

         call redplt (name,err)

         if (err) then

            if (iam.eq.14) then

               return

            else if (icopt.ne.5.or.iopt(34).eq.0) then 

               call error (72,nopt(1),i,'missing/corrupt plt/blk files '
     *                     //'VERTEX may still be running or the files'
     *                     //' are locked by another program')

            end if

         else 

            return

         end if

      end if
c                                 the only paths here are
c                                 1) iopt(34) = 2 => man.
c                                 2) iopt(34) = 1 => auto and no final results.
      inter = .false.
c                                 only hope is interim results:
      call mertxt (tfname,prject,'.irf',0)
      open (1000, file = tfname, status = 'old', iostat = ier)

      if (ier.ne.0) then 

         if (iopt(2).eq.1) then 
c                                  end of the line
            call error (72,nopt(1),i,'no IRF file: interim '//
     *                               'results are not available')
         else 
c                                  maybe the user deleted the irf file
            call warn (99,nopt(1),i,'no IRF file: interim '//
     *                              'results are not available')

            i = 0 

         end if

      else 
c                                 make a list of the result files
         i = 1
c                                 make a list of the result files
         do

            read (1000,*,iostat=ier) jnd(i,1),jnd(i,2)

            if (ier.ne.0) then

               if (i.eq.1) then

                  call error (72,nopt(1),i,'empty IRF file: interim '//
     *                                     'results are not available')

               else

                  i = i - 1
                  exit

               end if

            end if

            i = i + 1

         end do

      end if 

      if (iopt(34).eq.1) then 

         if (i.eq.0) then 
            write (*,'(a)') 'VERTEX has not completed the calculation '
     *                    //'and no interim results are available.'

            stop

         end if 
c                                 interim_results is auto, and the final results
c                                 are not available, find/use last interim result:
         write (*,'(a)') 'VERTEX has not completed the calculation, '//
     *                  'continue with the latest interim result (Y/N)?'

         if (refine.and.jnd(i,1).eq.0) write (*,'(2(/,a))')
     *      'WARNING: VERTEX is currently in, or was interrupted '//
     *      'during, the auto-refine stage, but the','latest interim '//
     *      'result is from the exploratory stage, the result may be '//
     *      'inconsistent or unreadable.'

         read (*,'(a)') yes

         if (yes.ne.'y'.and.yes.ne.'Y') stop

         write (text,'(a,i1,i1)') '_',jnd(i,1),jnd(i,2)
         call mertxt (name,prject,text,0)

      else
c                                 if here must be auto and an irf file exists
         if (i.gt.0) then 

            write (*,'(a)') 'Do you want to plot/analyze interim '//
     *                        'results (Y/N)?'
            read (*,'(a)') yes

            if (yes.eq.'y'.or.yes.eq.'Y') then
c                                 use intermediate results
               write (*,'(/,a,/)') 'Choose from the following interim'//
     *                             ' results [default is the last]:'

               do j = 1, i 

                  if (jnd(j,1).eq.0) then

                     write (*,'(4x,i1,a,i1)') j,
     *                      ' - exploratory stage, grid level ',jnd(j,2)
                  else

                     write (*,'(4x,i1,a,i1)') j,
     *                      ' - auto-refine stage, grid level ',jnd(j,2)

                  end if

               end do

               call rdnumb (nopt(1),0d0,i,i,.false.)
               write (*,'(/)')

               ind1 = jnd(i,1)
               ind2 = jnd(i,2)

               if (refine.and.ind1.eq.0) write (*,'(3(a,/))')
     *            'WARNING: VERTEX is in, or has completed, the '//
     *            'auto-refine stage, interim results ',
     *            'from the exploratory stage may be '//
     *            'inconsistent or unreadable.','if VERTEX has been '//
     *            'terminated and the next message is **error ver072'//
     *            '**, then edit T to F in the TOF file'

               write (text,'(a,i1,i1)') '_',ind1, ind2
               call mertxt (name,prject,text,0)
               inter = .true.

            else 

               name = prject

            end if

         else
c                                 use final results
            name = prject

         end if

      end if

      call redplt (name,err)

      if (err) then
         if (inter) then 
            call error (72,nopt(1),i,'corrupt interim results, '//
     *                             'use auto-refine stage results.')
         else
            call error (72,nopt(1),i,'missing/corrupt plt/blk files '
     *                     //'VERTEX may still be running or the files'
     *                     //' are locked by another program')
         end if
      end if

      end

      subroutine bplinp (err)
c-----------------------------------------------------------------------
c read the b-plot file that contains the information on the assemblages
c stable at each grid node
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical err

      integer jxco, kxco, i, j, ids, ier
c                                 -------------------------------------
c                                 global variables
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 global assemblage data
      integer icog,jcog
      common/ cxt17 /icog(k2),jcog(k2)

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)

      double precision bg
      common/ cxt19 /bg(k5,k2)

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jtest,jpot
      common/ debug /jtest,jpot

      double precision amu
      common/ cst48 /amu(k8,k2)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer iam
      common/ cst4 /iam

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
c                                 assemblage counter
      ibulk = 0
c                                 pointer to solution compositional coordinates
      jxco = 0 
      kxco = 0

      err = .false. 

      do 

         ibulk = ibulk + 1

         if (ibulk.gt.k2) call error (183,0d0,k2,'BLINP')

         read (n5,*,end=99) icog(ibulk),jcog(ibulk),iap(ibulk)

         ias = iap(ibulk)
c                                if ias = 0, probably reading 
c                                an inconsistent blk file in unsplt
         if (ias.le.0) then 
            ier = 1
            exit 
         end if 
c                                phase molar amounts
         read (n5,*,iostat=ier) (bg(i,ibulk),i=1,iavar(3,ias))
         if (ier.ne.0) goto 99

         ico(ibulk) = jxco

         do i = 1, iavar(1,ias)

            ids = idasls(i,ias)     

            kxco = jxco + ncoor(ids) 
            jxco = jxco + 1

            if (kxco.gt.k18) call error (61,0d0,k18,'BPLINP')

            read (n5,*,iostat=ier) (xco(j), j = jxco, kxco)
            if (ier.ne.0) goto 99

            if (lopt(32).and.ksmod(ids).eq.39) then 
c                                lagged speciation

               jxco = kxco + 1
               kxco = kxco + nat

               if (kxco.gt.k18) call error (61,0d0,k18,'BPLINP')

               read (n5,*,iostat=ier) (xco(j), j = jxco, kxco)
               if (ier.ne.0) goto 99

            end if  
         
            jxco = kxco

         end do 

         jxco = kxco  
c                                 read mu's if available
         if (jpot.ne.1) then
 
            read (n5,*,iostat=ier) (amu(i,ibulk), i = 1, jbulk)

            if (ier.ne.0) then 
c                                 if error on read most probably its
c                                 because of NaN's for the chemical 
c                                 potentials
               do i = 1, jbulk
                  amu(i,ibulk) = nopt(7)
               end do 
 
               ier = 0 

            end if 
         end if 

      end do

99    ibulk = ibulk - 1

      if (ier.ne.0) err = .true.

      end


      subroutine plinp (err)
c---------------------------------------------------------------------- 
c plinp - subroutine to read assemblage info for gridded min calculations.
c if icopt = 7 and fileio also reads nodal coordinates.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, jst, irep, kd, jend, ier

      logical count, err 

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer igrd
      common/ cst311/igrd(l7,l7)

      integer iam
      common/ cst4 /iam

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer idasls,iavar,iasct,ias
      common/ cst75 /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vip
      common/ cst28 /vip(l2,k2)

      character*100 cfname
      common/ cst227 /cfname

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      integer idstab,nstab,istab
      common/ cst34 /idstab(i11),nstab(i11),istab

      integer idsol,nrep,nph
      common/ cst38/idsol(k5,k3),nrep(k5,k3),nph(k3)
c----------------------------------------------------------------------
      err = .false.
c                                 top of plot file
      read (n4,*,iostat=ier) loopx, loopy, jinc
c                                 check if the file was generated by unsplt
c                                 if so set unsplt flag for sample_on_grid
      if (jinc.eq.-1) then 
         jinc = 1
         lopt(47) = .true.
      else 
         lopt(47) = .false.
      end if 

      if (ier.ne.0) goto 99
c                                 prior to 6.8.5 vertex did not write 
c                                 the final value of jinc to the plot 
c                                 file, reset it here for back-compatibility
      if (loopx.eq.1.or.loopy.eq.1) jinc = 1
c                                 decompress the grid data
      do i = 1, loopx, jinc
         jst = 1
         do while (jst.le.loopy)
            read (n4,*,iostat=ier) irep, kd
            if (ier.ne.0) goto 99
            if (kd.eq.0) write (*,*) 'bad un at i, j',i,j
            jend = jst + irep 
            do j = jst, jend
               if (j.gt.l7) call error (2,nopt(1),j,
     *                      'coordinates (routine PLINP), increase L7')
               igrd(i,j) = kd
            end do 
            jst = jend + 1
         end do 
      end do 
c                                 read assemblages
      read (n4,*,iostat=ier) iasct
      if (ier.ne.0) goto 99

      istab = 0 

      do i = 1, iasct
         read (n4,*,iostat=ier) iavar(1,i),iavar(2,i),iavar(3,i)
         if (ier.ne.0) goto 99
         read (n4,*,iostat=ier) (idasls(j,i), j = 1, iavar(3,i))
         if (ier.ne.0) goto 99
c                                 make a cumulative list of stable phases
c                                 first get the number of occurrences of 
c                                 each phase in the assemblage
         nph(i) = 0
         do j = 1, k5
            idsol(j,i) = 0
            nrep(j,i) = 0
         end do 

         do j = 1, iavar(3,i) 

            count = .true.

            if (j.le.iavar(1,i)) then 

               do k = 1, nph(i)
                  if (idsol(k,i).eq.idasls(j,i)) then 
                     count = .false.
                     nrep(k,i) = nrep(k,i) + 1
                     exit 
                  end if 
               end do
 
            end if 

            if (count) then
               nph(i) = nph(i) + 1
               idsol(nph(i),i) = idasls(j,i)
               nrep(nph(i),i) = 1
            end if

         end do  
c                                 make an array in which each id 
c                                 occurs only once
         
c                                 next compare to the existing list
         do k = 1, nph(i)  

            count = .true.

            do j = 1, istab

               if (idsol(k,i).eq.idstab(j)) then 
                  if (nrep(k,i).gt.nstab(j)) nstab(j) = nrep(k,i)
                  count = .false.
                  exit
               end if 

            end do 

            if (count) then 
               istab = istab + 1
               if (istab.gt.k10) call error (999,0d0,istab,'ISTAB ')
               nstab(istab) = nrep(k,i)
               idstab(istab) = idsol(k,i)
            end if 

         end do 

      end do 
c                                 make the "null" assemblage
      iap(k2) = k3
      iavar(1,k3) = 0
      iavar(2,k3) = 0 
      iavar(3,k3) = 0 

      if (icopt.eq.7.and.fileio) then 
c                                 if coodinates from a file, read
c                                 coordinate file.
         open (n8,file=cfname,status='old',iostat=ier)
         if (ier.ne.0) call error (6,vip(1,1),i,cfname)
         if (loopy.gt.k2) call error (1,vip(1,1),loopy,'k2')
         do j = 1, loopy
            read (n8,*,iostat=ier) (vip(i,j), i = 1, ipot)
            if (ier.ne.0) then 
               write (*,1000) cfname
               stop
            end if 
         end do 
         close (n8)

      end if 

99    if (ier.ne.0) err = .true.

1000  format (/,'**error ver635** Coordinate file ',a,/,
     *       'is inconsistent with plot file, re-run VERTEX.',/)

      end


      subroutine getvar  
c--------------------------------------------------------------------
c getvar makes a list of variables to be used for i/o:

c if icopt = 10 -> using nodal coordinates else, 

c if icopt =  9/11 -> using 2d frac coordinates else:

c one-dimensional diagram (oned = .true.) then 

c the vertical (real) axis is variable iv(2), the horizontal axis
c is dummy.

c two-dimensional diagram (oned = .false.) and not icopt = 9, then 

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables

c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,lay,mpol,mord

      parameter (lay=6,mpol=3,mord=4) 

      integer iam
      common/ cst4 /iam

      logical pzfunc
      integer ilay,irep,npoly,nord
      double precision abc0,vz,iblk
      common/ cst66 /abc0(0:mord,mpol),vz(6),iblk(lay,k5),ilay,
     *               irep(lay),npoly,nord,pzfunc

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)   

      character vname*8, xname*8
      common / csta2 /xname(k5),vname(l2)  

      logical oned
      common/ cst82 /oned

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------

      if (icopt.eq.7.and.fileio) then 
c                                 1d-fractionation with file input,
c                                 use nodal coordinates:
         vnm(1) = 'node #'
         vmn(1) = 1
         vmx(2) = 1d0
         vmn(2) = 0d0 
         vmx(1) = loopy 
         oned = .true.
         
         jvar = ipot + 1

         do i = 2, jvar
            vnm(i) = vname(jv(i-1))
         end do  

      else if (icopt.lt.9) then 

         jvar = ipot

         if (idep.ne.0) jvar = ipot + 1

         if (icont.eq.1) then 

            do i = 1, jvar
               vnm(i) = vname(jv(i))
               vmx(i) = vmax(jv(i))
               vmn(i) = vmin(jv(i))
               var(i) = vmin(jv(i))
            end do   

         else 

            if (icont.eq.2) then 

               jvar = jvar + 1

               vnm(1) = ' X(C1)  '
               vmx(1) = 1d0
               vmn(1) = 0d0

               do i = 2, jvar
                  vnm(i) = vname(jv(i-1))
                  vmx(i) = vmax(jv(i-1))
                  vmn(i) = vmin(jv(i-1))
                  var(i) = vmin(jv(i-1))
               end do   
 
            else 

               jvar = jvar + 2

               vnm(1) = ' X(C1)  '
               vmx(1) = 1d0
               vmn(1) = 0d0

               vnm(2) = ' X(C2)  '
               vmx(2) = 1d0
               vmn(2) = 0d0

               do i = 3, jvar
                  vnm(i) = vname(jv(i-2))
                  vmx(i) = vmax(jv(i-2))
                  vmn(i) = vmin(jv(i-2))
                  var(i) = vmin(jv(i-2))
               end do   

            end if 

         end if 

         if (oned) then 
c                                 make a fake y-axis for 1-d plots
            vmx(2) = 1d0
            vmn(2) = 0d0

         end if

      else if (icopt.eq.9) then 
c                                using non-thermodynamic coordinate frame
         vmn(1) = vz(4)
         vmx(1) = vz(5)

         if (iam.ne.1) then
c                                 vertex sets the number of y-nodes in frac2d
            ncol = loopy
         else 
c                                 pssect/werami get the number of nodes from the plt file
            loopy = ncol

         end if 

         if (flsh) then
c                                  flush calculations: 
            vnm(1) = 'Q,kg/m^2'
            vnm(2) = 'dz,m   '
c                                  set the base to 
            vmn(2) = vz(1)/2d0
            vmx(2) = vmn(2) + dfloat(ncol-1)*vz(1)

         else
c                                  frac2d calculations.
            vnm(1) = 'z0,m'
            vnm(2) = 'dz,m'
c                                  set y = 0 ti be the top
            vmx(2) = -vz(1)/2d0
            vmn(2) = vmx(2) - dfloat(ncol-1)*vz(1)

         end if

         jvar = 4

         do i = 3, 4
            vnm(i) = vname(jv(i-2))
         end do

      else if (icopt.eq.12) then 

         vnm(1) = 'n,alqt. '
         vnm(2) = 'node#      '

         vmn(2) = 1d0
         vmx(2) = dfloat(iopt(36)) + 1d0
         var(2) = 1d0
 
         vmn(1) = 0d0
         vmx(1) = nopt(36)*dfloat(iopt(36))
         var(1) = 0d0

         v(1) = vmin(1)
         v(2) = vmin(2)

         jvar = ipot + 2

         do i = 3, jvar
            vnm(i) = vname(jv(i-2))
            vmx(i) = vmax(jv(i-2))
            vmn(i) = vmin(jv(i-2))
            var(i) = vmin(jv(i-2))
         end do

      end if 

      end
      
      
            subroutine grxn (gval)
c--------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer j

      double precision fo2, fs2, gval, gphase

      external gphase

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c---------------------------------------------------------------------
c                                 compute free energy change of the rxn
      gval = 0d0
      fh2o = 0d0
      fco2 = 0d0
 
      if (ifct.gt.0) call cfluid (fo2,fs2)
 
      do j = 1, iphct
         gval = gval + vnu(j) * (gphase(j) + r * t * dlog(act(j)))
      end do 
 
      gval = gval + vuf(1) * fh2o*r*t + vuf(2)*fco2*r*t

      if (idf(3).ne.0.and.ifct.gt.0) gval = gval + vnu(idf(3))*r*t*fo2
 
      end
      
      
      
