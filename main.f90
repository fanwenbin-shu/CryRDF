module para
use constant
implicit none

! file operation
character(icha) struFileName

! main parameter
integer :: iter ! iteration, the amount of computing RDF
character(icha), allocatable :: struFile(:)  ! structure file name
character(icha), allocatable :: ipair(:,:)   ! two atoms in each structure
character(icha) :: fileName ! current file name
integer :: px, py, pz ! extra periodical, default p* = 1
real(inum) :: bw ! bin width, default bw = 0.1d0

! crystal information
logical :: icrystal ! 1: periodical box, 0: no box
real(inum) :: v(3,3) ! xyz \times abc
real(inum) :: volume, rcutoff, meanDensity ! \rho
integer :: Natoms
real(inum), allocatable :: coord(:,:)
character(2), allocatable :: eleName(:)
real(inum), allocatable :: eleMass(:)
integer, allocatable :: eleNum(:)
real(inum) :: InertiaTensor(3,3)

! RDF
integer :: scatterCount
real(inum), allocatable :: scatterPoints(:)
integer :: Aatoms, Batoms
integer, parameter :: derCountMax = 8

end module



program main
!use para
implicit none

call ReadInputFile
call Compute

write(*,*)
write(*,*) 'All jobs done! '

end program


subroutine Compute
use para
implicit none
integer i

write(*,*)
do i = 1, iter
write(*,*)
write(*,'(A,I4,A)') ' Calculating ', i, ' structure...'
fileName = struFile(i)
call crystalCoord
call crystalInfo
call atomCount(trim(ipair(1,i)), trim(ipair(2,i)))
call computeRadialScatter(trim(ipair(1,i)), trim(ipair(2,i)))
call computeDistribution(trim(ipair(1,i)), trim(ipair(2,i)))

end do

end subroutine


subroutine computeDistribution(A, B)
use para
implicit none
character(len=*), intent(in) :: A, B
real(inum), allocatable :: RDFx(:), RDFy(:), RDFint(:)
integer :: RDFlength, i, loc
real(inum), external :: sphereShellVolume
character(icha) :: RDFFileName, RDFIntFileName
real(inum) :: der, derOld
integer :: derCount

    write(*,*)
    write(*,*) 'Computing RDF from scatter points... '
    RDFFileName = trim(fileName)//'_'//trim(A)//'_'//trim(B)//'_RDF.txt'
    open(723, file=trim(RDFFileName), status='replace', action='write') ! RDF -> 723

    ! create radius array
    RDFlength = ceiling(rcutoff / bw) + 1
    allocate( RDFx(RDFlength), RDFy(RDFlength), RDFint(RDFlength) )
    do i = 1, RDFlength
        RDFx(i) = (i-1) * bw
    end do

    ! histogram
    RDFy(1:RDFlength) = 0d0
    do i = 1, scatterCount
        call hunt(RDFx, RDFlength, scatterPoints(i), int(floor(scatterPoints(i)/bw)), loc)
        loc = loc - 1
        RDFy(loc) = RDFy(loc) + 1d0
    end do
    
    ! divide by volume
    do i = 1, RDFlength - 2 ! remove first and last points
        RDFy(i) = RDFy(i) / sphereShellVolume(RDFx(i), bw) / meanDensity
        write(723,'(F12.4,F18.10)') RDFx(i), RDFy(i)
    end do
    close(723)
    
    ! integration using trapezoidal rule
    write(*,*) 'Integrating RDF... '
    RDFIntFileName = trim(fileName)//'_'//trim(A)//'_'//trim(B)//'_RDF_int.txt'
    open(467, file=trim(RDFIntFileName), status='replace', action='write') ! int -> 467
    RDFint(1:RDFlength) = 0d0
    ! $\int_0^0 f(x)\mathrm dx = 0$, therefore the fisrt integration value must be zero. 
    write(467,'(F12.4,F18.10)') RDFx(1), RDFint(1)
    ! Let's integrate from the second. 
    ! formula : $\int r^2 g(r) \mathrm dr$
    do i = 2, RDFlength - 2
        RDFint(i) = RDFint(i-1) + & 
                    !(RDFy(i) * RDFx(i) * RDFx(i) + RDFy(i-1) * RDFx(i-1) * RDFx(i-1)) & 
                    !* bw * 0.5d0
                    RDFy(i) * RDFx(i) * RDFx(i) * bw
        write(467,'(F12.4,F18.10)') RDFx(i), RDFint(i)
    end do
    close(467)
    
    ! compute derivates
    derOld = 0d0
    derCount = 0
    write(*,*) 'Peak of RDF (No., Location)'
    do i = 2, RDFlength - 2
        der = ( RDFy(i) - RDFy(i-1) ) / bw
        if ( derOld .gt. 0d0 .and. der .lt. 0d0 ) then
            derCount = derCount + 1
            write(*,'(I3,F9.3,E12.3)') derCount, RDFx(i-1)!, der
        end if
        derOld = der
        if ( derCount .eq. derCountMax ) exit
    end do
    if ( derCount .eq. 0 ) write(*,*) 'No peak found! '

end subroutine

function sphereShellVolume(R, dr)
use constant
implicit none
real(inum), intent(in) :: R, dr
real(inum) sphereShellVolume

    ! formula : 
    ! $V(R+\mathrm dr) - V(R) = \frac43 \pi \mathrm dr (3R^2 + 3\mathrm drR + \mathrm dr^2)$
    sphereShellVolume = (3d0*R*R + 3d0*R*dr + dr*dr) * 4d0*pi*dr / 3d0

end function

subroutine computeRadialScatter(A, B)
use para
implicit none
character(len=*), intent(in) :: A, B
character(icha) :: scatterFileName
integer :: i,j,k,ii
integer :: lx, ly, lz ! in periodical conditions
real(inum), external :: dist
real(inum) :: shiftedCoord(3), shiftValue(3)
real(inum) rm ! r in mirror image

    !scatterFileName = trim(fileName)//'_'//trim(A)//'_'//trim(B)//'_scatter.txt'
    !open(728, file=trim(scatterFileName), status='replace', action='write')
    write(*,*)
    write(*,*) 'Start calculating RDF scatter... '
    
    allocate( scatterPoints(Aatoms*Batoms) )
    scatterCount = 0

    do i = 1, Natoms ! A atom, center atom
        if ( trim(adjustl(eleName(i))) == trim(adjustl(A)) ) then
            do lx = -px, px
            do ly = -py, py
            do lz = -pz, pz
            
                do ii = 1, 3
                    shiftValue(ii) = v(ii, 1) * lx + v(ii, 2) * ly + v(ii, 3) * lz
                end do
                
                do j = 1, Natoms
                    if ( trim(adjustl(eleName(j))) == trim(adjustl(B)) ) then
                        do k = 1, 3 ! xyz
                            shiftedCoord(k) = coord(k,j) + shiftValue(k)
                        end do
                        rm = dist(shiftedCoord(1:3), coord(1:3,i))
                        if (rm .lt. rcutoff .and. rm .gt. 1d-5) then
                            scatterCount = scatterCount + 1
                            scatterPoints(scatterCount) = rm 
                        end if
                    end if
                end do
                
            end do
            end do
            end do
        end if
    end do    
    !close(728)
    
    write(*,'(A,I10)') ' Total effective RDF scatter points : ', scatterCount

end subroutine


function dist(a, b)
use constant
implicit none
real(inum), intent(in) :: a(3), b(3)
real(inum) :: dist
real(inum) :: tmp(3)

tmp(1:3) = a(1:3) - b(1:3)
dist = dot_product(tmp(1:3), tmp(1:3))
dist = dsqrt(dist)

end function

subroutine atomCount(A, B)
use para
implicit none
character(len=*), intent(in) :: A, B
integer :: i
    
    Aatoms = 0; Batoms = 0;
    if ( trim(A) .eq. trim(B) ) then
        do i = 1, Natoms
            if ( trim(adjustl(eleName(i))) .eq. trim(A) ) Aatoms = Aatoms + 1
        end do
        Batoms = Aatoms
    else
        do i = 1, Natoms
            if ( trim(adjustl(eleName(i))) .eq. trim(A) ) then
                Aatoms = Aatoms + 1
            else if ( trim(adjustl(eleName(i))) .eq. trim(B) ) then
                Batoms = Batoms + 1
            end if
        end do
    end if
     
    write(*,*) 'The number of each atoms : '
    write(*,'(X,2A,I6,5X,2A,I6)') trim(A), ' : ', Aatoms, trim(B), ' : ', Batoms
    
    ! mean density
    ! formula : $\rho = N_A N_B / V$
    if (icrystal) then
        meanDensity = Aatoms * Batoms / volume
    else
        meanDensity = Aatoms! * Batoms
    end if
    write(*,'(A,F18.6)') ' Mean density : ', meanDensity

end subroutine


subroutine crystalCoord
use para
implicit none
integer fileStatus, atom
character(icha) tempLine, temp

    Natoms = 0
    open(732, file=trim(fileName), status='old', action='read') ! pdb -> 732
    
    ! count the total number of atoms
    do while (.true.)
        read(732,'(A)',iostat=fileStatus) tempLine
        if (tempLine(1:4) .eq. 'ATOM' .or. tempLine(1:6) .eq. 'HETATM') then
            Natoms = Natoms + 1
        end if
        if (fileStatus /= 0) exit
    end do
    if (Natoms .le. 1) then
        write(*,'(A,I2,3A)') ' There only ', Natoms, ' atom in `', trim(fileName), '`.'
        call logError('Too few atoms were detected! ')
    end if
    write(*,*)
    write(*,'(A,I8)') ' The number of total atoms : ', Natoms
    allocate(coord(3,Natoms), eleName(Natoms), eleMass(Natoms), eleNum(Natoms))
    rewind(732)

    ! read coordinates
    atom = 1
    do while (.true.)
        read(732,'(A)',iostat=fileStatus) tempLine
        if (tempLine(1:4) .eq. 'ATOM' .or. tempLine(1:6) .eq. 'HETATM') then
            read(tempLine,'(A30, 3F8.3, A22, A2)') &
                temp, coord(1,atom), coord(2,atom), coord(3,atom), temp, eleName(atom)
            atom = atom + 1
        end if
        if (fileStatus /= 0) exit
    end do
    close(732)
    
    ! get mass, index, and normalize all name
    call eleNameParse

end subroutine


subroutine crystalInfo
use para
implicit none
integer fileStatus
character(icha) tempLine, temp
! cystal info
real(inum) :: a, b, c, alpha, beta, gamma
real(inum) :: cosa, cos2a, cosb, cos2b, cosg, cos2g, sing, Delta
real(inum) :: r(3) ! distance from center to each plane

    open(732, file=trim(fileName), status='old', action='read') ! pdb -> 732
    do while (.true.)
        read(732,'(A)',iostat=fileStatus) tempLine
        if (tempLine(1:6) .eq. 'CRYST1') then
            read(tempLine, '(A6,3F9.3,3F7.2)') temp, a, b, c, alpha, beta, gamma
            icrystal = .true.
            exit
        end if
        if (fileStatus /= 0) then
            write(*,*) 'This .pdb file does not contain periodical conditions. '
            icrystal = .false.
            close(732)
            call computeInertiaTensor
            return
        end if
    end do
    close(732)
    
    write(*,*)
    write(*,*) 'Structure infomation : '
    write(*,*) '     |a|      |b|      |c|   alpha    beta   gamma'
    write(*,'(3F9.3,3F8.2)') a, b, c, alpha, beta, gamma
    
    ! judge the rationality of gamma
    if (gamma .gt. alpha + beta .or. gamma .lt. alpha - beta) then
        write(*,'(A,F7.2,A,F7.2,A)') &
            ' [NOTE] The \gamma must be from ', alpha - beta, ' to ', alpha + beta, ' degree. '
        call logError('The \gamma is impossible! Check the three angles! ')
    end if
    
    ! calculate trigonometric function
    alpha = alpha * pi / 180d0
    beta  = beta  * pi / 180d0
    gamma = gamma * pi / 180d0
    cosa = dcos(alpha); cos2a = dcos(2d0 * alpha)
    cosb = dcos(beta); cos2b = dcos(2d0 * beta)
    cosg = dcos(gamma); cos2g = dcos(2d0 * gamma); sing = dsin(gamma)
    Delta = - cos2a - cos2b - cos2g + 4d0 * cosa * cosb * cosg - 1d0
    if (Delta .lt. 0d0) then
        call logError('Delta is negetive! ')
    else
        Delta = dsqrt(Delta)
    end if
    
    ! calculate lattice vector
    v(1:3, 1:3) = 0d0
    v(1,1) = a
    v(1,2) = b * cosg; v(2,2) = b * sing;
    v(1,3) = dsqrt(2d0) * c * cosb * sing / Delta
    v(2,3) = dsqrt(2d0) * c * (cosa - cosb * cosg) / Delta
    !v(3,3) = c * dsqrt(-cos2a-cos2b+4d0*cosa*cosb*cosg-1d0) ! less efficient
    v(3,3) = dsqrt(c**2d0 - v(1,3) * v(1,3) - v(2,3) * v(2,3))
    write(*,*)
    write(*,*) 'Basic vector : '
    write(*,*) '\vec      x/Ang       y/Ang       z/Ang'
    write(*,'(A2, 2X, 3F12.6)') 'a', v(1:3,1)
    write(*,'(A2, 2X, 3F12.6)') 'b', v(1:3,2)
    write(*,'(A2, 2X, 3F12.6)') 'c', v(1:3,3)
    
    ! calculate the volume of the cell
    volume = a * b * sing * v(3,3)
    write(*,*) ! formula : $V = |\vec a| |\vec b| \sin\gamma \ c_z$
    write(*,'(A,F15.6)') ' Volume of the cell (Angstrom^3): ', volume
    
    ! compute distance to each plane
    call computeMinDistance(v(1:3,1:3), r(1:3))
    rcutoff = minval(r(1:3))
    
    write(*,*)
    write(*,*) 'Distances from center to each plane: '
    write(*,*) '         bc          ac          ab'
    write(*,'(3F12.6)') r(1:3)
    write(*,'(2(A,F12.6))') ' r_cutoff : ', rcutoff, ', r^2 = ', rcutoff * rcutoff

end subroutine


subroutine computeMinDistance(v, r)
use constant
implicit none
real(inum), intent(in) :: v(3,3)
real(inum), intent(out) :: r(3)
real(inum) :: O(3), rc(3) ! origin and center of cell

    O(1:3) = 0d0
    rc(1:3) = ( v(1:3,1) + v(1:3,2) + v(1:3,3) ) * 0.5d0
    call computeVerticalDistance(O(1:3), v(1:3,2), v(1:3,3), rc(1:3), r(1))
    call computeVerticalDistance(O(1:3), v(1:3,1), v(1:3,3), rc(1:3), r(2))
    call computeVerticalDistance(O(1:3), v(1:3,1), v(1:3,2), rc(1:3), r(3))

end subroutine


subroutine computeVerticalDistance(a, b, c, rc, r)
use constant
implicit none
real(inum), intent(in) :: a(3), b(3), c(3), rc(3)
real(inum), intent(out) :: r
real(inum) :: AA, BB, CC, DD

    call computePlaneCoeff(a,b,c,AA,BB,CC,DD)
    r = abs(AA * rc(1) + BB * rc(2) + CC * rc(3) + DD)
    r = r / dsqrt(AA * AA + BB * BB + CC * CC)

end subroutine

subroutine computePlaneCoeff(a,b,c,AA,BB,CC,DD)
use constant
implicit none
real(inum), intent(in) :: a(3), b(3), c(3)
real(inum) :: AA, BB, CC, DD, maxCoeff

    AA = (b(2) - a(2)) * (c(3) - a(3)) - (c(2) - a(2)) * (b(3) - a(3))
    BB = (b(3) - a(3)) * (c(1) - a(1)) - (c(3) - a(3)) * (b(1) - a(1))
    CC = (b(1) - a(1)) * (c(2) - a(2)) - (c(1) - a(1)) * (b(2) - a(2))
    DD = - (AA * a(1) + BB * a(2) + CC * a(3))
    
    if (abs(AA+BB+CC+DD) .lt. 1D-25) then
        write(*,*) '          A           B           C           D'
        write(*,'(4E12.2)') AA,BB,CC,DD
        call logError('The coefficients of plane could not be zero at same time! ')
    end if
    
    maxCoeff = maxval(abs( (/AA, BB, CC, DD/) ))
    AA = AA / maxCoeff
    BB = BB / maxCoeff
    CC = CC / maxCoeff
    DD = DD / maxCoeff

end subroutine


subroutine ReadInputFile
use para
implicit none
character(icha) inputFileName, tempLine
integer inputFileStatus, inputFileReadStatus, inputFileLineCount
integer i

    write(*,*) 'CryRDF'
    write(*,*) 'Author : Wenbin FAN (fanwenbin@shu.edu.cn)'
    write(*,*) 

    call getarg(1, inputFileName)
    inquire(file=trim(inputFileName), exist=inputFileStatus)
    if (inputFileStatus) then
        open(478, file=trim(inputFileName), status='old', action='read') ! input -> ipt -> 478
        write(*,*) 'Input file name : ', trim(inputFileName)
    else
        call logError("Your input file `" // trim(inputFileName) // "` doesn't exist! ")
    end if

    ! counting the amount of line
    inputFileLineCount = 0
    do while (.true.)
        read(478, *, iostat=inputFileReadStatus)
        if (inputFileReadStatus/=0) exit
        inputFileLineCount = inputFileLineCount + 1
    end do
    rewind(478) ! go back to the beginning of input file
    
    if (inputFileLineCount .gt. 0) then
        if (mod(inputFileLineCount, 2) .eq. 0) then
            iter = inputFileLineCount / 2
            
            write(*,'(A,I4)') ' Number of structure(s) : ', iter
            allocate(struFile(iter), ipair(2, iter))
            
            write(*,*)
            write(*,*) 'Structure             Atom 1  Atom 2'
            write(*,*) '--------------------  ------  ------'
            
            do i = 1, iter
                read(478, *) struFile(i)
                read(478, *) ipair(1:2,i)
                
                inquire(file=trim(struFile(i)), exist=inputFileStatus)
                if (not(inputFileStatus)) then
                    call logError('Your structure file `' // trim(struFile(i)) // "` doesn't exist! ")
                end if
                
                write(*,'(A20,2X,A2,6X,A2)') trim(struFile(i)), trim(ipair(1,i)), trim(ipair(2,i))
            end do
        else
            call logError("Check contents in your input file `" // trim(inputFileName) // "`! ")
        end if
    else
        call logError("Your input file `" // trim(inputFileName) // "` is blank! ")
    end if
    close(478)
    
    ! read calculating parameters (optional)
    write(*,*) 
    inquire(file='parameter', exist=inputFileStatus)
    if (inputFileStatus) then
        open(7272, file='parameter', status='old', action='read') ! para -> 7272
        write(*,*) '`parameter` file exists. Loading parameters...'
        read(7272,*) px, py, pz
        read(7272,*) bw
    else
        write(*,*) 'Default parameters loaded. '
        px = 1; py = 1; pz = 1
        bw = 0.1d0
    end if
    close(7272)
    
    write(*,'(A,3I3)') 'Extra periodical box in x,y,z : ', px, py, pz
    write(*,'(A,F12.6)') 'Bin width : ', bw
    
end subroutine

! calculate inertia tensor
subroutine computeInertiaTensor
use para
implicit none
integer i
real(inum) :: m,x,y,z
real(inum) :: ITEValue(3), ITEV(3,3) ! eigenvector of inertia tensor
real(inum) :: totalMass, COM(3), rg, maxBox(3)

    ! total mass
    totalMass = sum(eleMass(1:Natoms))
    write(*,'(A,F18.3)') ' Total mass (amu): ', totalMass
    ! center of mass
    do i = 1,3
        COM(i) = sum( coord(i,1:Natoms) * eleMass(1:Natoms) ) / totalMass
    end do
    write(*,'(A,3F12.6)') ' Center of mass (X,Y,Z): ', COM(1:3)
    ! shifted to the origin
    do i = 1,3
        coord(i,:) = coord(i,:) - COM(i)
    end do
    
    ! radius of gyration
    rg = sum( eleMass(:)*( coord(1,:)**2 + coord(2,:)**2 + coord(3,:)**2 ) )
    rg = dsqrt( rg / totalMass)
    write(*,'(A,F12.6)') ' Radius of gyration (Angstrom): ', rg

    InertiaTensor(1:3,1:3) = 0d0
    do i = 1, Natoms
        m = eleMass(i)
        x = coord(1,i); y = coord(2,i); z = coord(3,i)
        InertiaTensor(1,1) = InertiaTensor(1,1) + m * (y*y + z*z)
        InertiaTensor(2,1) = InertiaTensor(2,1) - m * x*y
        InertiaTensor(3,1) = InertiaTensor(3,1) - m * x*z
        InertiaTensor(1,2) = InertiaTensor(1,2) - m * y*x
        InertiaTensor(2,2) = InertiaTensor(2,2) + m * (x*x + z*z)
        InertiaTensor(3,2) = InertiaTensor(3,2) - m * y*z
        InertiaTensor(1,3) = InertiaTensor(1,3) - m * z*x
        InertiaTensor(2,3) = InertiaTensor(2,3) - m * z*y
        InertiaTensor(3,3) = InertiaTensor(3,3) + m * (x*x + y*y)
    end do
    
    !write(*,*) 'Inertia tensor : '
    !do i = 1,3
    !    write(*,'(3E12.3)') InertiaTensor(:,i)
    !end do
    
    call calcEigenVector(3, InertiaTensor, ITEValue, ITEV)
    !write(*,*) 'Eigenvector of inertia tensor : '
    !do i = 1,3
    !    write(*,'(3F12.6)') ITEV(:,i)
    !end do
    
    ! rotate molecule and export to .xyz
    open(999, file=trim(fileName)//'_rotated.xyz',status='replace',action='write')
    write(999, *) Natoms
    write(999, *) 'Symmetrised geometry generated by CryRDF - Wenbin, FAN'
    do i = 1, Natoms
        coord(1:3,i) = matmul( coord(1:3,i), ITEV(1:3,1:3) )
        write(999, '(A2, 3F12.6)') trim(eleName(i)), coord(1:3,i)
    end do
    close(999)
    
    ! measure the max box
    do i = 1,3
        maxBox(i) = maxval(coord(i,:)) - minval(coord(i,:)) + 4d0 !vdW radius
    end do
    write(*,'(A,3F12.6)') ' Max box size (plused 4 Ang): ', maxBox(1:3)
    rcutoff = dsqrt(sum(maxBox(1:3)*maxBox(1:3))) * 0.5d0 ! diameter -> radius
    write(*,'(A,F12.6)') ' Maximum of r : ', rcutoff
    px = 0; py = 0; pz = 0; volume = 1d0

end subroutine

! parse all element name to number, mass, and default name
subroutine eleNameParse
use para
implicit none
integer i
integer, external :: eleName2Num

    do i = 1, Natoms
        eleNum(i) = eleName2Num( trim(adjustl(eleName(i))) )
        eleMass(i) = allMass( eleNum(i) )
        eleName(i) = trim( allName( eleNum(i) ) ) ! default name
    end do

end subroutine

! get number of element
function eleName2Num(nameIn)
use constant
implicit none
character(len=*), intent(in) :: nameIn
integer :: eleName2Num
integer j

    do j = 1, supNum
        if ( trim(nameIn) == trim(allName(j)) ) then
            eleName2Num = j
            return
        else if ( trim(nameIn) == trim(allNameUp(j)) ) then
            eleName2Num = j
            return
        end if
    end do
    
    if ( j .gt. supNum ) then
        call logError( 'Element not supported: '//trim(nameIn) )
    end if

end function

! log error and exit whole program
subroutine logError(content)
implicit none
character(len=*), intent(in) :: content

    print *, '[ERROR] ', content
    stop

end subroutine

subroutine logWarn(content)
implicit none
character(len=*), intent(in) :: content

    print *, '[WARNING] ', content

end subroutine