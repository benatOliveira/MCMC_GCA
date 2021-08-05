      program meemmum_call_fun       
c----------------------------------------------------------------------
      character*14 pname_o(1,12)!! k5 = 12
      double precision cblk_i(1,12)!! k5 = 12
      double precision v_i(1,2)!!should be of size 2 (not l2) because otherwise it detects numerical 0's  (~1E-300) as input data
      double precision psys_o(1,10+12) !! array of numerical output properties of bulk system  
      double precision props_o(12,12+12) !! array of numerical output properties per phase  
      double precision iPs

      v_i(1,1) = 19967.015
      v_i(1,2) = 1599.0901
      cblk_i(1,1) = 44.760000
      cblk_i(1,2) = 3.9800000
      cblk_i(1,3) = 8.1900000
      cblk_i(1,4) = 0.1300000
      cblk_i(1,5) = 38.770000
      cblk_i(1,6) = 3.1700000
      cblk_i(1,7) = 0.2800000
      cblk_i(1,8) = 0.5700000
      cblk_i(1,9) = 0.1300000

      write(*,'(a)') 'Thermodynamic potentials (P;T)'
      print *, v_i
      write(*,'(a)') 'Chemical concentrations'
      print *, cblk_i

      call meemum_fun(pname_o,props_o,psys_o,iPs,v_i,cblk_i)
      
      write(*,*) 'Number of stable phases', iPs
      write(*,'(a)') 'Psys out'
      print *, psys_o
      write(*,'(a)') 'Thermodynamic phases'
      print *, pname_o
      write(*,'(a)') 'Props out'
      print *, props_o

      end 
