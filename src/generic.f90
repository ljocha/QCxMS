module qcxms_use_generic 
   use common1
   use readcommon
   use qcxms_utility, only: getspin
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only: toSymbol 
   implicit none

   contains

     subroutine genericout(nat,xyz,iat,chrg,spin,etemp,grad,ECP)

        integer ::  nat
        integer ::  chrg
        integer ::  iat(nat)
        integer ::  i,spin
        integer ::  io_in

        real(wp) :: xyz (3,nat)
        real(wp) :: etemp

        logical :: grad,wrb
        logical :: ECP

        ECP = .False.
        if(spin == 0) call getspin(nat,iat,chrg,spin)

        open(file='generic.in',newunit=io_in)

        write(io_in,'(''energy'')')
        if(grad) write(io_in,'(''gradient'')')
   
        write(io_in,'(''* xyz ''2i3)')chrg,spin
        do i=1,nat
           write(io_in,300)toSymbol(iat(i)),xyz(1,i)*autoaa,xyz(2,i)*autoaa,xyz(3,i)*autoaa
        enddo
     300 format(a2,3(f22.12))
        write(io_in,'(''*'')')
        close(io_in)
    
     end subroutine genericout
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     subroutine rdgenericgrad(fname,nat,g,achrg,aspin,edum)

        integer :: nat
        integer :: i,j,nn
        integer :: io_out

        real(wp) :: g(3,nat),edum
        real(wp) :: achrg(nat),aspin(nat)
        real(wp) :: xx(20)
     
        character(len=*)  :: fname

        logical :: ex
     
        open(file=fname, newunit=io_out, status='old')
        do
           read ( unit=io_out, FMT='(a)', iostat=iocheck ) line
           if ( iocheck > 0 ) stop 'error in rdgenericgrad!'
           if ( iocheck < 0 ) exit !EOF

           if ( index(line,'energy') /= 0 ) then
              call readl(line,xx,nn)
              edum=xx(1)
           endif

           if ( index(line,'charges and spin populations') /= 0 ) then
              do i=1, nat
                 read(unit=io_out,FMT='(a)') line
                 call readl(line,xx,nn)
                 achrg(i) = xx(2)
                 if ( nn == 3 ) aspin(i) = xx(3)
              enddo
           endif

           if ( index(line,'gradient') /= 0 ) then
              do i=1, nat
                 read(unit=io_out,FMT='(a)') line
                 call readl(line,xx,nn)
                 g(1,i) = xx(2)
                 g(2,i) = xx(3)
                 g(3,i) = xx(3)
              enddo
           endif
        enddo 

        close(io_out)
     end subroutine rdgenericgrad
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     subroutine rdgenericen(fname,energy)

        integer :: nn
        integer :: io_out

        real(wp) :: xx(20)
        real(wp) :: energy
     
        character(len=*)  :: fname
     
        open(file=fname, newunit=io_out, status='old')

        do
          read (unit=io_out, FMT='(a)', iostat=iocheck) line
          if ( iocheck > 0 ) stop 'error in rdgenericen!'
          if ( iocheck < 0 ) exit ! EOF

          if (index(line,'energy') /= 0 ) then
             call readl(line,xx,nn)
             energy=xx(1)
             exit
          endif

        enddo

        close(io_out)
     
     end subroutine rdgenericen



end module qcxms_use_generic
