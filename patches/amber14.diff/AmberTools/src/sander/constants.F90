! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++
!This module contains various parameters
!and constants used by the different 
!routines that make up sander.
!
!If you want to use one of the constants
!in your routine you should include the
!line:
!
!use constants, only : xxx, yyy, zzz
!
!where xxx,yyy,zzz are the constants you plan
!to use in your routine.
!This line needs to go before the
!implicit none declaration.
!
! Based on constants.h, a pre Fortran 90 version, by Scott Brozell
!   and Dave Case (TSRI, 2002)
! Converted into a Fortran 90 module by: Ross Walker (TSRI, 2005)
! Expanded by others including: Matthew Clark, Andreas Goetz,
!
!++++++++++++++++++++++++++++++++++++++++

module constants

  implicit none

  ! by default everything in this module is public
  public

  !------------------------------------------------------------
  ! Generic Floating Point Constants
  _REAL_, parameter :: TEN_TO_MINUS2  = 1.0d-2
  _REAL_, parameter :: TEN_TO_MINUS3  = 1.0d-3
  _REAL_, parameter :: TEN_TO_MINUS4  = 1.0d-4
  _REAL_, parameter :: TEN_TO_MINUS5  = 1.0d-5
  _REAL_, parameter :: TEN_TO_MINUS6  = 1.0d-6
  _REAL_, parameter :: TEN_TO_MINUS8  = 1.0d-8
  _REAL_, parameter :: TEN_TO_MINUS10 = 1.0d-10
  _REAL_, parameter :: TEN_TO_MINUS25 = 1.0d-25
  _REAL_, parameter :: TEN_TO_PLUS3   = 1.0d+3
  _REAL_, parameter :: TEN_TO_PLUS10  = 1.0d+10

  _REAL_, parameter :: zero      = 0.0d0
  _REAL_, parameter :: one       = 1.0d0
  _REAL_, parameter :: two       = 2.0d0
  _REAL_, parameter :: three     = 3.0d0
  _REAL_, parameter :: four      = 4.0d0
  _REAL_, parameter :: five      = 5.0d0
  _REAL_, parameter :: six       = 6.0d0
  _REAL_, parameter :: seven     = 7.0d0
  _REAL_, parameter :: eight     = 8.0d0
  _REAL_, parameter :: nine      = 9.0d0
  _REAL_, parameter :: ten       = 10.0d0
  _REAL_, parameter :: eleven    = 11.0d0
  _REAL_, parameter :: twelve    = 12.0d0
  _REAL_, parameter :: sixteen   = 16.0d0
  _REAL_, parameter :: twenty    = 20.0d0
  _REAL_, parameter :: thirtytwo = 32.0d0
  _REAL_, parameter :: sixtyfour = 64.0d0

  _REAL_, parameter :: half         = one/two
  _REAL_, parameter :: third        = one/three
  _REAL_, parameter :: fourth       = one/four
  _REAL_, parameter :: fifth        = one/five
  _REAL_, parameter :: sixth        = one/six
  _REAL_, parameter :: seventh      = one/seven
  _REAL_, parameter :: eighth       = one/eight
  _REAL_, parameter :: ninth        = one/nine
  _REAL_, parameter :: tenth        = one/ten
  _REAL_, parameter :: eleventh     = one/eleven
  _REAL_, parameter :: twelfth      = one/twelve
  _REAL_, parameter :: sixteenth    = one/sixteen
  _REAL_, parameter :: thirtysecond = one/thirtytwo
  _REAL_, parameter :: sixtyfourth  = one/sixtyfour
  
  _REAL_, parameter :: thirtieth    = one/30.0d0

  !------------------------------------------------------------
  !     THE ARRAY FC(I) CONTAINS THE FACTORIALS OF (I-1).

  _REAL_, parameter :: FC(1:17) =&
       (/ 1.0D0,1.0D0, 2.0D0, 6.0D0, 24.0D0, 120.0D0, 720.0D0, 5040.0D0, &
          40320.0D0, 362880.0D0, 3628800.0D0, 39916800.0D0,              &
          4.790016D+08, 6.2270208D+09, 8.71782912D+10,                   &
          1.307674368D+12, 2.092278989D+13 /)   
  
  _REAL_, parameter :: logFC(1:17) = (/ 0.0D0, 0.0D0, 0.6931471805599D0,       &
       &            1.7917594692281D0,  3.1780538303479D0,  4.7874917427820D0, &
       &            6.5792512120101D0,  8.5251613610654D0, 10.6046029027453D0, &
       &           12.8018274800815D0, 15.1044125730755D0, 17.5023078458739D0, &
       &           19.9872144956619D0, 22.5521638531234D0, 25.1912211827387D0, &
       &           27.8992713838409D0, 30.6718601061763D0 /)

  !     DEFINE C COEFFICIENTS FOR ASSOCIATE LEGENDRE POLYNOMIALS.         
  _REAL_, parameter::CC(1:21,1:3) = reshape ( (/   &      
       8.0D0,   8.0D0,   4.0D0,  -4.0D0,  4.0D0,   &
       4.0D0, -12.0D0,  -6.0D0,  20.0D0,  5.0D0,   & 
       3.0D0, -30.0D0, -10.0D0,  35.0D0,  7.0D0,   &     
       15.0D0,   7.5D0, -70.0D0, -17.5D0, 63.0D0,  &
       10.5D0,                                     &               
       0.0D0,   0.0D0,   0.0D0,  12.0D0,  0.0D0,   &
       0.0D0,  20.0D0,  30.0D0,   0.0D0,  0.0D0,   &       
       -30.0D0,  70.0D0,  70.0D0,   0.0D0,  0.0D0, &
       -70.0D0, -105.D0, 210.0D0, 157.5D0,  0.0D0, &
       0.0D0,                                      &     
       0.0D0,   0.0D0,   0.0D0,   0.0D0,  0.0D0,   &
       0.0D0,   0.0D0,   0.0D0,   0.0D0,  0.0D0,   &
       35.0D0,   0.0D0,   0.0D0,   0.0D0,  0.0D0,  &
       63.0D0, 157.5D0,   0.0D0,   0.0D0,  0.0D0,  &
       0.0D0/), (/ 21, 3 /) )      
          
  !------------------------------------------------------------
  ! Physical Constants
  _REAL_, parameter :: LIGHT_SPEED = 2.997924d08
  _REAL_, parameter :: HBAR = 627.509d0 * 0.0241888d-3 * 20.455d0 !Planck's constant in internal units
  _REAL_, parameter :: J_PER_CAL = 4.184d0 !  This is defined as the thermochemical calorie
  _REAL_, parameter :: JPKC = J_PER_CAL * 1000.0d0 !kilocalories per joule
  _REAL_, parameter :: BOLTZMANN = 1.380658d-23 !Boltzmann's constant in J/K
  _REAL_, parameter :: AVOGADRO = 6.0221367d+23 !Avogadro's number
  _REAL_, parameter :: KB = (BOLTZMANN * AVOGADRO) / JPKC !Boltzmann's constant in internal units
  _REAL_, parameter :: AMBER_ELECTROSTATIC = 18.2223d0
  _REAL_, parameter :: AMBER_ELECTROSTATIC2 = AMBER_ELECTROSTATIC * AMBER_ELECTROSTATIC
  !Ratio by which to scale amber charges to get electron charges - amberchg * oneqscale = electron charges
  ! = 1.0 / 18.2223d0
  _REAL_, parameter :: INV_AMBER_ELECTROSTATIC = 1.0d0/AMBER_ELECTROSTATIC
  _REAL_, parameter :: INV_AMBER_ELECTROSTATIC2 = 1.0d0/AMBER_ELECTROSTATIC2

  _REAL_, parameter :: CHARGE_ON_ELEC = 1.60217733d-19 !Charge on an electron in Coulombs
  _REAL_, parameter :: BOHR_RADIUS = 52.9177249d-12 ! in meter
  _REAL_, parameter :: BOHRS_TO_A = 0.529177249D0   ! Bohrs * this = angstroms - Same constants as used in dynamo v2.
  !_REAL_, parameter :: BOHRS_TO_A = 0.52917706D0   ! Bohrs * this = angstroms - Same constants as used in Gaussian 98
  !_REAL_, parameter :: BOHRS_TO_A = 0.529177D0     ! Bohrs * this = angstroms - Same constants as used in Mopac6 hcore.f
  !_REAL_, parameter :: BOHRS_TO_A = 0.529167D0     !                            as used in Mopac6 repp.f
  _REAL_, parameter :: A_TO_BOHRS = 1.0d0 / BOHRS_TO_A
  !_REAL_, parameter :: A_TO_BOHRS = 1.88976D0      !Same constants as used in Mopac6 gover.f
  _REAL_, parameter :: A2_TO_BOHRS2 = A_TO_BOHRS * A_TO_BOHRS !Mopac6 uses 3.5711928576D0 in gover.f for this.
  _REAL_, parameter :: A3_TO_BOHRS3 = A2_TO_BOHRS2 * A_TO_BOHRS
  _REAL_, parameter :: A4_TO_BOHRS4 = A2_TO_BOHRS2 * A2_TO_BOHRS2

  _REAL_, parameter :: ONE_AU = 27.2113962d0 !One atomic unit of energy in eV.
  _REAL_, parameter :: HARTREE_TO_JOULE = ONE_AU * CHARGE_ON_ELEC !conversion from hartrees to joules
  _REAL_, parameter :: HART_BOHR_TO_JOULE_A = HARTREE_TO_JOULE * BOHRS_TO_A !hartree*bohrs to joules*angstroms
  _REAL_, parameter :: COULOMB_CONST_E = HART_BOHR_TO_JOULE_A*AVOGADRO/JPKC
                                                !Coulomb's constant for charges in units of e
                                                !This is the same as AMBER_ELECTROSTATIC2 but to higher precision
  !_REAL_, parameter :: AU_TO_EV = ONE_AU  !Conversion from AU to EV - not used because we match dynamo v2 below.
  _REAL_, parameter :: AU_TO_EV = 27.21d0 !Conversion from AU to EV - Same as dynamo v2 uses and Gaussian 98
                                          !Note (RCW+MC): more precise would be: 1 a.u. 27.211396 eV
                                          !Mopac6 uses 27.21D0 in calpar.f, delri.f and repp.f but in
                                          !ffhpol.f it uses 27.2107 and in the manual it quotes 27.211
  _REAL_, parameter :: HALF_AU_TO_EV   = AU_TO_EV * half
  _REAL_, parameter :: FOURTH_AU_TO_EV = AU_TO_EV * fourth
  _REAL_, parameter :: EIGHTH_AU_TO_EV = AU_TO_EV * eighth
  _REAL_, parameter :: SXNTH_AU_TO_EV  = EIGHTH_AU_TO_EV * half
  _REAL_, parameter :: A2_TO_BOHRS2xAU_TO_EV = A2_TO_BOHRS2*AU_TO_EV

  !_REAL_, parameter :: EV_TO_KCAL = 23.060362D0      !Conversion from EV to KCAL/MOL
  !Dynamo parameter
  _REAL_, parameter :: EV_TO_KCAL = 23.061d0  !Dynamo's conversion
                                              !Mopac6 uses 23.061 in ffhpol.f analyt.f compfg.f datin.f dcart.f
                                              !                      delri1.f delri2.f deritr.f interp.f iter.f
                                              !                      moldat.f mopac.f
  _REAL_, parameter :: KCAL_TO_EV = one / EV_TO_KCAL

  _REAL_, parameter :: AU_TO_KCAL = AU_TO_EV*EV_TO_KCAL !1 hartree.


  ! The following are updated constants from Mohr, Taylor, Newell, Rev. Mod. Phys. 80 (2008) 633-730.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !       NAME                    REFERENCE               VALUE               UNITS    !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Avogadro constant        (Table XLVIII, p.708)   6.0221417930E23         mol^-1    !
  ! Bohr radius              (Table XLIX, p. 710)    0.5291772085936E-10     m         !
  ! a.u. of energy           (Table LII, p. 717)     4.3597439422E-18        J         !
  ! speed of light (vacuum)  (Table I, p.637)        299 792 458             m*s^−1    !
  ! a.u. of charge           (Table LIII, p. 717)    1.60217648740E-19       C         !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  _REAL_, parameter :: CODATA08_AVOGADRO       = 6.0221417930d23      ! Avogadro's number
  _REAL_, parameter :: CODATA08_BOHR_RADIUS    = 0.5291772085936d-10
  _REAL_, parameter :: CODATA08_ONE_AU         = 4.3597439422d-18     ! Atomic unit of energy in joules
  _REAL_, parameter :: CODATA08_LIGHT_SPEED    = 2.99792458d08
  _REAL_, parameter :: CODATA08_CHARGE_ON_ELEC = 1.60217648740d-19
  ! Derived values
  _REAL_, parameter :: CODATA08_A_TO_BOHRS     = 1d-10 / CODATA08_BOHR_RADIUS
  _REAL_, parameter :: CODATA08_AU_TO_KCAL     = CODATA08_ONE_AU / J_PER_CAL / 1000 * CODATA08_AVOGADRO
  _REAL_, parameter :: CODATA08_AU_TO_DEBYE    = &
                         CODATA08_LIGHT_SPEED * CODATA08_CHARGE_ON_ELEC * CODATA08_BOHR_RADIUS / 1d-21
  !_REAL_, parameter :: AU_TO_DEBYE = 1.0d0/0.393430 ! from http://cccbdb.nist.gov/debye.asp (April 12 2011)

  !------------------------------------------------------------
  !Numeric Constants
  _REAL_, parameter :: PI      = 3.1415926535897932384626433832795d0

  !The BOOK says :
  !
  !2Chronicles 4:2 reads thus, 'Also he made a molten sea of ten cubits 
  !from brim to brim, round in compass, and five cubits the height thereof; 
  !and a line of thirty cubits did compass it round about.'
  !
  !Hence, Pi is exactly equal to three and there is nothing more to discuss!
  !
  !If you want to use the value of PI defined by 'the BOOK' then uncomment
  !the following line and comment out the definition above...
  !_REAL_, parameter :: PI = 3.0d0

  _REAL_, parameter :: PI2     = PI*PI
  _REAL_, parameter :: HALFPI = PI * 0.5d0
  _REAL_, parameter :: TWOPI  = 2.0d0 * PI
  _REAL_, parameter :: FOURPI = 4.0d0 * PI
  _REAL_, parameter :: INVPI  = 1.0d0 / PI
  _REAL_, parameter :: SQRTPI = 1.77245385090551602729816748334d0 !sqrt(PI)
  _REAL_, parameter :: INVSQRTPI = 1.0d0 / SQRTPI
  _REAL_, parameter :: DEG_TO_RAD = PI / 180.0d0
  _REAL_, parameter :: RAD_TO_DEG = 180.0d0 / PI 
  _REAL_, parameter :: LN_TO_LOG = 2.30258509299404568402d0  ! log(1.0d1)

  _REAL_, parameter :: SQRT2     = 1.4142135623730950488016887242097d0
  _REAL_, parameter :: INVSQRT2  = 1.0d0 / SQRT2

  !------------------------------------------------------------
  !Generalised Born Constants
  _REAL_, parameter :: alpb_alpha = 0.571412d0 !Alpha prefactor for alpb_alpha

  !------------------------------------------------------------
  ! Unusual Constants
  integer, parameter :: RETIRED_INPUT_OPTION = -10301 ! first 5 digit palindromic prime
  integer, parameter :: NO_INPUT_VALUE = 12344321  ! from Bob Duke
  _REAL_, parameter  :: NO_INPUT_VALUE_FLOAT = 12344321.d0

  integer :: plumed
  character(256) :: plumedfile

contains

  function BinomialCoefficient(m, n) result (bioCoeff)
        
    integer, intent(in)::m,n
    _REAL_::bioCoeff
        
    integer, parameter::size=30
    integer, save::bc(size,size)=0
    logical, save::initialized=.false.;
        
    integer::i,j,k
        
    if (.not.initialized) then
       do i=1,size
          bc(i,1)=one
          bc(i,2:size)=zero
       end do
       do i=2, size
          do j=2, i
             bc(i,j)=bc(i-1,j-1)+bc(i-1,j)
          end do
       end do
       
       initialized=.true.
    end if
        
    bioCoeff=one*bc(m,n)
    
  end function BinomialCoefficient

end module constants

