!==============================================================================
	  module parameters
             implicit none
             integer, parameter:: Np = 14   ! length of a partition
             integer, parameter:: max_domin_parti=72000 ! number of dominated partitions
             integer, parameter:: nb_particles = 4! how many particles (variables)
             integer, parameter:: nb_orbit =14 ! how many orbitals on the sphere ! = Np?
      end module parameters
!=========================
      module formats
      		character(LEN=6), parameter :: party = '(40I1)'   !LEN!!!
      		character(LEN=6), parameter :: party2 = '(11I3)'   !LEN!!!
      		character(LEN=20), parameter :: maci = '(8F6.1)'
      		character(LEN=20), parameter :: wspp = '(F15.10)'
      end module formats

!============================================================================== END Zmienne + Formaty


!============================================================================== FUNKCJE
      module funkcje
      contains
!============================================ rho

      double precision Function rho(lambda)
      use parameters
      use formats
      implicit none
      integer lambda(0:Np)
      integer i,j

      integer how_many_variabl_l, sum_sq, sum_i, index
	  double precision alpha

      common /alp/ alpha


      sum_i = 0
      sum_sq = 0
      how_many_variabl_l =0
      index = 1  ! For sums in LLL

!============================== Particle number counting
      do i=0,Np
         how_many_variabl_l = how_many_variabl_l + lambda(i)
      enddo

      if (how_many_variabl_l .ne. nb_particles) then
          write(*,*) 'Problem with rho - number of particles different than variables'
          write(*,*) 'particles:', nb_particles
          write(*,*) 'variables:', how_many_variabl_l

          write(*,party) lambda
          read(*,*)
          call EXIT
      endif
!============================== END: Particle number counting

      do i=0,Np
         sum_sq = sum_sq + i*i*lambda(i)
      enddo

      do i=0,Np
         j = 0
         do while (lambda(Np-i).ne. j)
                 sum_i = sum_i + index*(Np-i)
                 index = index + 1
                 j = j+1
         enddo
      enddo


      rho = dble(sum_sq) + dble(sum_i)*2d0*(1d0 - 1d0 / alpha) ! PLUS - MINUS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      return
      end function



!============================================ END rho


!============================================ C_{\lambda, \mu} Coefficients

      integer Function C_lm (theta,mu)

      use parameters
      implicit none

      integer theta(0:Np), mu(0:Np)
      integer down(1:2), up(1:2)     ! numbers that are changed while acting with raising operator
      integer tmp
      integer ell ! difference of changed numbers
      integer isdown, isup ! how many numbers would need to be changed to obtain mu from lambda by lowering )
      integer i
      integer swap

      swap = 0

      isdown = 0
      isup = 0
      down(1:2)= 0
      up(1:2)= 0

      do tmp = 0,Np

         if(theta(tmp) .ne. mu(tmp)) then

                        !How many parts of lambda are bigger them mu
            select case (theta(tmp) - mu(tmp))
!                   case (0)
						! not going to happen
                   case (1)
                        isup = isup+1

			if (isup > 2) then
                        	C_lm = 0
                            goto 111
                        endif

                        up(isup) = tmp
                   case (2)
                        isup = isup+2

			if (isup > 2) then
                        	C_lm = 0
                            goto 111
                        endif

                        down(1) = tmp
                        down(2) = tmp
                   case (-1)
                        isdown=isdown+1

			if (isdown > 2) then
                        	C_lm = 0
                            goto 111
                        endif

                        down(isdown)=tmp
                   case (-2)
                        isdown = isdown+2

			if (isdown > 2) then
                        	C_lm = 0
                            goto 111
                        endif

                        down(1) = tmp
                        down(2) = tmp
                   case default
                        C_lm = 0
                        goto 111
            end select

            if (isdown > 2 .or. isup > 2) then !if more than one can not be obtain by one raising operation
               C_lm = 0
               goto 111
            endif

         endif

      enddo

      if (up(1)-down(1) .eq. down(2)-up(2)) then

          C_lm = (down(2)-down(1))

          do i = up(1)+1, down(1)-1
             swap = swap + mu(i)
          enddo

          do i = down(2)+1, up(2)-1
             swap = swap + mu(i)
          enddo

          if (mod(swap,2) .eq. 1) then
             C_lm = (-1)*C_lm
          endif

      else
          C_lm = 0
      endif

111   return
      end Function


!============================================ END  C_{\lambda, \mu} Coefficients


!============================================ Is partition dominated in natural order by another partition?


      Logical Function is_dominated(lambda, mu)
      use parameters
      implicit none


      integer lambda(0:Np)
      integer mu(0:Np)
      integer i,j
      integer sumL, sumM
      integer partL, partM
      integer howmanyL, howmanyM
      integer how_many_variablL, how_many_variablM

      how_many_variablL = 0
      how_many_variablM = 0

      sumL = 0
      sumM = 0

      do i=0,Np
         how_many_variablL = how_many_variablL + lambda(i)
      enddo

      do i=0,Np
         how_many_variablM = how_many_variablM + mu(i)
      enddo
!	   Wont work for notsqueezingable partition cause mu(0)=-1
      if (mu(0) .ne. -1 .and. how_many_variablL .ne. how_many_variablM) then
         write(*,*) 'partitions have different length'
         write(*,*) lambda
         write(*,*)
         write(*,*) mu
         read(*,*)
      endif



      do i = 0, Np
         if (lambda(Np-i) .ne. 0) then
            partL = Np-i
            howmanyL = lambda(Np-i)
            EXIT
         endif
      enddo


      do i = 0, Np
         if (mu(Np-i) .ne. 0) then
            partM = Np-i
            howmanyM = mu(Np-i)
            EXIT
         endif
      enddo

      do j = 1, how_many_variablL



             if (howmanyL .eq. 0) then
                do i = (Np - partL +1) ,Np
                   if (lambda(Np-i).ne. 0) then
                      partL = Np-i
                      howmanyL = lambda(Np-i)
                      EXIT
                   endif
                enddo
             endif


             if (howmanyM .eq. 0) then
                do i = (Np - partM +1) ,Np
                   if (mu(Np-i).ne. 0) then
                      partM = Np-i
                      howmanyM = mu(Np-i)
                      EXIT
                   endif
                enddo
             endif

             sumL = sumL + partL
             sumM = sumM + partM


             if (sumL < sumM) then
                is_dominated = .false.


                goto 222
             endif



             howmanyL = howmanyL -1
             howmanyM = howmanyM -1



      enddo




      is_dominated = .true.

222   return
      end function
!============================================ END Is partition dominated in natural order by another partition?



!============================================ Next partition in natural order


      Function dominated(lambda) result(mu)
      use parameters
      implicit none

      integer lambda(0:Np)
      integer mu(0:Np)
      integer i,j
      integer how_many_variabl,how_many_to_fill  ! How many variables are available , parts of partition (to fill)
      integer flag_1
      flag_1 =-1

      mu(0:Np) = lambda(0:Np)

      how_many_variabl = lambda(0) + lambda(1)
      how_many_to_fill = lambda(1)


      do i = 2,Np
         if(lambda(i) .ne. 0) then

            if ( (i-1)*how_many_variabl  .ge. (how_many_to_fill+1)) then !<=> (how_many_variabl+1)(i-1) .ge. how_many_to_fill+i
                flag_1 =1 ! partition is changing
                mu(i)=mu(i)-1
                mu(i-1) = (how_many_to_fill+1)/(i-1) + 1 ! +1 cause one is from mu(i) ! <=> (how_many_to_fill+i)/(i-1)

                if ((how_many_to_fill + 1 - ((mu(i-1) -1 ) * (i-1))) .eq. 0) then ! how many variables one has to use ! <=> how_many_to_fill + i - mu(i-1) *(i-1)  = 0
                   do j = 0,(i-2)
                      mu(j)=0
                   enddo
                    mu(0)= (how_many_variabl+1) - mu(i-1)
                   goto 333
                else
                    do j = 0,(i-2)
                       mu(j)=0
                    enddo
                    mu(how_many_to_fill+1 - ((mu(i-1)-1) * (i-1)) ) = 1 ! zero everything exept for number that has left
                    mu(0)= how_many_variabl - mu(i-1) ! ile zmiennych niewykorzystanych idzie do zera
                    goto 333
                endif
            else
                how_many_to_fill = how_many_to_fill+ lambda(i)*i
                how_many_variabl =  how_many_variabl + lambda(i)
            endif

         endif

      enddo



      if (flag_1 .eq. -1) then
         mu(0) = -1
      endif


333   return
      end function

!============================================ END Next partition in natural order


!============================================ List of partitions dominated by fixed partition


      Subroutine  list_zdom_B (list_domin_B, lambda)
      use parameters
      implicit none

      integer list_domin_B(0:Np,0:max_domin_parti)
      integer lambda(0:Np)
      integer Msize
      integer i
      integer temppart (0:Np)
      integer temppart2 (0:Np)



      Msize = -1
      list_domin_B(0:Np,0:max_domin_parti) = 0
      list_domin_B(0:Np,1) = lambda   ! We are starting with lambda

      do i = 2,max_domin_parti

 	 temppart2 =  list_domin_B(0:Np, i-1)
         temppart(0:Np) = dominated(temppart2) !lambda!dominated(lista(0:Np, i-1))

555      if ( .not. (is_dominated(lambda, temppart)) ) then
             temppart(0:Np) = dominated(temppart)
             goto 555
         endif

         if  (temppart(0) .eq. -1 ) then   ! can partition be squezzed? if not = -1
            Msize = i-1
            goto 444
         endif

         list_domin_B(0:Np, i)  = temppart
      enddo


      if (Msize .eq. -1) then
         write(*,*) 'Increase max partitions'
         read(*,*)
      endif

444   list_domin_B(0,0) = Msize
	  write(*,*) 'Size of partition matrix:  ', Msize


      return
      end Subroutine



!============================================ END List of partitions dominated by fixed partition

!============================================ List of partitions dominated by fixed partition


      Subroutine  list_zdom_F (list_domin_B, list_domin_F, Msize)
      use parameters
      implicit none

      integer list_domin_F(0:Np,0:max_domin_parti)
      integer list_domin_B(0:Np,0:max_domin_parti)
      integer temppart (0:Np)
      integer index
      integer i,j,k
      integer Msize

      list_domin_F(0:Np,0:max_domin_parti) = 0

      do i = 1, Msize
         index = 0
         do j = 0, Np
            k = 0
            do while (list_domin_B(j,i) .ne. k)

                 list_domin_F(j+index,i) =  list_domin_F(j+index,i) + 1

                 index = index + 1
                 k = k + 1
            enddo
         enddo
      enddo

      return
      end Subroutine



!============================================ END List of partitions dominated by fixed partition



!============================================ Coefficients (on the plane)

      Subroutine wsp(Msize, list_domin_F, c_i)
      use parameters
      implicit none

      integer Msize

      double precision prod, suma
      integer i,j

      double precision rho_lam
      double precision c_i(0:max_domin_parti)
      integer list_domin_F(0:Np,0:max_domin_parti)
      integer temppart (0:Np)
      integer temppart2 (0:Np)
      double precision alpha

      common /alp/ alpha


      c_i(0:max_domin_parti) = 0d0
      c_i(1) = 1d0

      prod = 1d0
      temppart(0:Np) =   list_domin_F(0:Np,1)
      rho_lam = rho(temppart)


      do i = 2, Msize

         temppart = list_domin_F(0:Np,i)

         do j = 1, i-1
            temppart2 =   list_domin_F(0:Np,j)
            c_i(i)  = c_i(i)+ c_i(j) * C_lm(temppart2, temppart)
         enddo

         if (abs(rho_lam - rho(temppart)) < 0.0000001d0) then
            write(*,*) 'Eigenvalues equal: Jack has a pole'
            read(*,*)
         endif

         c_i(i)=  c_i(i)*2d0*(1d0/alpha-1d0) /(rho_lam - rho(temppart))
      enddo

      c_i(0)=real(Msize)

      write(*,*) 'Max and min in c_i: ', maxval(real(c_i)), minval(real(c_i))
      write(*,*)



      return
      end Subroutine

!============================================ END Coefficients (on the plane)

!============================================ Monopole Harmonics - Coefficients

      Subroutine  MH_S (MH, orbitals)
      use parameters
      implicit none

      double precision MH(0:Np) ! Monopole harmonica coeffiecients (inverse stereographic projection)
      integer orbitals ! Np?

      integer i
      MH(0:Np) =0d0
      MH(0) = 1d0

      MH(orbitals) = MH(0)

      do i = 1, (orbitals)/2
         MH(i) = (MH(i-1))*sqrt(real(i) / real(orbitals+1-i))
         MH(orbitals-i) = MH(i)
      enddo

      write(*,*) 'Max and min in MH_S: ', maxval(real(MH)), minval(real(MH))
      write(*,*)

      return
      end Subroutine

!============================================ END Monopole Harmonics - Coefficients

!============================================ Jack: from plane to sphere

      Subroutine  coeff_Sp (MH, c_i, b_i, list_domin_F,Msize )
      use parameters
      implicit none

      double precision MH(0:Np)
      double precision c_i(0:max_domin_parti)
      double precision b_i(0:max_domin_parti)
      integer list_domin_F(0:Np,0:max_domin_parti)

      integer Msize
      integer i,j
      double precision sphere
      integer temppart(0:Np)


      do i = 1, Msize

      	 temppart =  list_domin_F(0:Np, i)
	 sphere = 1d0

	 do j = 0, Np
	    if(temppart(j) .ne. 0) then
            	sphere = sphere * (MH(j))**(temppart(j))
            endif
         enddo

	 b_i(i) = c_i(i) * (sphere)

      enddo

      write(*,*) 'Max and min in b_i: ', maxval(real(b_i)), minval(real(b_i))
      write(*,*)

      return
      end Subroutine

!============================================ END Jack: from plane to sphere
!============================================ Normalization

      Subroutine  nor (b_i, Msize)

      use parameters
      implicit none

      double precision b_i(0:max_domin_parti)
      integer Msize
      double precision  suma

      integer i

      suma = 0d0

      do i = 1, Msize
	suma = suma + b_i(i)*b_i(i)
      enddo

      write(*,*) 'SUM: ', suma
      suma = sqrt(suma)

      do i = 1, Msize
	 b_i(i) = b_i(i)/suma
      enddo

      return
      end Subroutine

!============================================ END Normalization
!============================================ Ordering

      Subroutine  order (ordering, b_i, list_domin_F, Msize)

      use parameters
      use formats
      implicit none

      integer ordering(0:max_domin_parti)
      double precision b_i(0:max_domin_parti)
      integer list_domin_F(0:Np,0:max_domin_parti)
      integer Msize
      integer temppart (0:Np)

      integer i,j, suma,max,position
      double precision alpha

      common /alp/ alpha

      open(1,file='FJ.out',status='unknown')
      open(2,file='alpha.in',status='unknown')
      open(3,file='lambda.in',status='unknown')

      ordering(0:max_domin_parti) = 0



      suma = 0
      max = 0
      position = 0

      do i = 1, Msize

      	 temppart =  list_domin_F(0:Np, i)
         suma =0

         do j = 0, Np
            suma = suma + 2**(Np-j) * abs( temppart(j)-1 )
         enddo

         ordering (i) = suma
      enddo

      write(1,*) '***Fermionic Jack polynomial (sphere)***'
            write(1,*) ' '
      write(1,*) 'Partition:'
      write(1,party) list_domin_F(0:Np,1)
      write(1,*) ' '
      write(1,*) 'Jack parameter (alpha):'
      write(1,*) alpha
      write(1,*)
      write(1,*) 'Dim. of the Hilbert subspace (Slater det.):', Msize
      write(1,*)
      write(1,*) 'Coefficients of Jack fermionic polynomial (partition, coeff):'
      write(1,*)

      do j = 1,Msize
         max = 0
         do i = 1,Msize

            if (max < ordering(i) ) then
               max = ordering(i)
               position = i
            endif

         enddo

         ordering(position) = -1

         write(1,party, advance="no" )  list_domin_F(0:Np, position)
         write(1,*) b_i(position)


      enddo



      return
      end Subroutine

!============================================ END Ordering



!============================================ zero

      Subroutine  zero (lambda,c_i, b_i, MH, list_domin_B,list_domin_F)

      use parameters
      use formats
      implicit none

      integer ordering(0:max_domin_parti)
      integer list_domin_F(0:Np,0:max_domin_parti)
      integer list_domin_B(0:Np,0:max_domin_parti)
      integer lambda(0:Np)
      integer Msize
      double precision c_i(0:max_domin_parti) ! coefficients of Jack (plane)
      double precision MH(0:Np)
      double precision b_i(0:max_domin_parti) !
      integer temppart (0:Np)

      integer i,j, suma,max,position,io
      double precision alpha

      common /alp/ alpha

      open(1,file='FJ.out',status='unknown')
      open(2,file='alpha.in',status='unknown')
      open(3,file='lambda.in',status='unknown')

      lambda(0:Np) =0


         READ(3, *) lambda(0:Np)




      c_i(0:max_domin_parti) = 0.0
      b_i(0:max_domin_parti) = 0.0
      MH(0:Np)  = 0.0

!==================================

      read(2,*) alpha


      return
      end Subroutine

!============================================ END zero






      end module funkcje



!============================================================================== Program JACK ON SPHERE


      Program JackSphere
      use parameters
      use funkcje
      use formats
      implicit none

      integer lambda(0:Np)
      integer mu (0:Np)

      integer Msize

      integer list_domin_B(0:Np,0:max_domin_parti)
      integer list_domin_F(0:Np,0:max_domin_parti)
      double precision c_i(0:max_domin_parti) ! coefficients of Jack (plane)
      double precision MH(0:Np)
      double precision b_i(0:max_domin_parti) ! coefficients of Jack (sphere)
      integer i ,io
      integer ordering(0:max_domin_parti)
      double precision alpha

      common /alp/ alpha


      open(1,file='FJ.out',status='unknown')
      open(2,file='alpha.in',status='unknown')
      open(3,file='lambda.in',status='unknown')








      !do i=0,Np
      !   READ(3, *, iostat=io) lambda(0:i)
      !   if(io==0) exit
      !   write(*,*) '--------------------------'
      !   write(*,*) i
      !   write(*,*) lambda
      !   read(*,*)
      !enddo



      list_domin_B(0:Np,0:max_domin_parti) = 0
      list_domin_F(0:Np,0:max_domin_parti) = 0




      call zero (lambda,c_i, b_i, MH, list_domin_B,list_domin_F)
      call list_zdom_B(list_domin_B, lambda)
      Msize = list_domin_B(0,0)
      call list_zdom_F(list_domin_B, list_domin_F, Msize)
      call wsp(Msize, list_domin_F, c_i)
      call MH_S (MH, nb_orbit)
      call coeff_Sp (MH, c_i, b_i, list_domin_F,Msize)
      call nor (b_i, Msize)
      call order (ordering, b_i, list_domin_F, Msize)





      End Program


!============================================================================== END Program JACK ON SPHERE
