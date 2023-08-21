module MultiFluidDE
  use DarkEnergyInterface
  use results
  use constants
  use classes
  implicit none

  private
  ! JVR - can we use these quantities from results.f90 without needing to bring them to this file via the Init subroutine?
  real(dl) :: grhocrit
  real(dl) :: grhode_today
  real(dl) :: Omega_de
  real(dl), parameter :: Tpl= sqrt(kappa*hbar/c**5)  ! sqrt(8 pi G hbar/c^5), reduced planck time
  real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  ! Convert to units of 1/Mpc^2

  type, extends(TDarkEnergyModel) :: TMultiFluidDE
    integer :: DebugLevel
    integer, dimension(max_num_of_fluids) :: models
    ! JVR - TODO go back to a parameter array
    ! JVR - model parameters
    real(dl) :: w0, wa ! CPL parameters
    real(dl) :: w1, w2, w3, z1, z2, z3 ! Binned w
    real(dl) :: zc, fde_zc, theta_i, wn ! Fluid EDE parameters
    real(dl) :: n, grhonode_zc, freq ! Fluid EDE internal parameters
    ! JVR - Internal variables, scalar field variables
    logical :: use_zc
    integer :: which_potential ! 1 for rock 'n' roll, 2 for axion 
    real(dl) :: V0 = 1d-115, m = 5d-54, f = 0.05, initial_phi = 1._dl ! Scalar field EDE parameters
    real(dl) :: astart = 1e-6_dl ! Scalar field integration start
    real(dl), dimension(:), allocatable :: sampled_a, phi_a, phidot_a
    logical :: did_scalar_field_init = .false.
    ! Steps for log a and linear spacing, switching at max_a_log (set by Init)
    integer :: npoints_linear, npoints_log, npoints = 200
    real(dl), private :: dloga, da, log_astart, max_a_log
    real(dl) :: integrate_tol=1d-3, min_steps_per_osc=10, zcfdeprec=1d-3
    real(dl), dimension(:), private, allocatable :: ddphi_a, ddphidot_a
    real(dl), dimension(:), allocatable :: fde, ddfde
    class(CAMBdata), private, pointer :: State
  contains
  procedure :: Init => TMultiFluidDE_Init
  procedure :: w_de => TMultiFluidDE_w_de
  procedure :: grho_de => TMultiFluidDE_grho_de
  procedure :: BackgroundDensityAndPressure => TMultiFluidDE_density
  procedure :: Effective_w_wa => TMultiFluidDE_Effective_w_wa
  procedure :: PerturbedStressEnergy => TMultiFluidDE_PerturbedStressEnergy
  procedure :: PerturbationEvolve => TMultiFluidDE_PerturbationEvolve
  procedure :: PerturbationInitial => TMultiFluidDE_PerturbationInitial
  procedure :: diff_rhopi_Add_Term => TMultiFluidDE_diff_rhopi_Add_Term
  procedure, nopass :: PythonClass => TMultiFluidDE_PythonClass
  procedure, nopass :: SelfPointer => TMultiFluidDE_SelfPointer
  procedure :: ReadParams => TMultiFluidDE_ReadParams
  procedure :: w_derivative
  procedure :: KG_Equations
  procedure :: KG_Equations_Log
  procedure :: Vofphi
  procedure :: Field_Values_At_a
  procedure :: check_error
  procedure :: calc_zc_fde
  procedure :: fde_peak
  end type TMultiFluidDE

  procedure(TClassDverk) :: dverk

  public TMultiFluidDE

  contains

  ! --------------------- Background Equations ---------------------

  function TMultiFluidDE_w_de(this, a) result (w_de)
    class(TMultiFluidDE) :: this
    real(dl), dimension(max_num_of_fluids) :: w_de
    real(dl), intent(in) :: a
    integer :: i
    real(dl) :: zc, ac, fde_zc, wn, z, phi, phidot, K, V

    ! i - counts the component: 1 is the late component, 2 is the early component.
    w_de = 0
    do i = 1, this%num_of_components
      if (i == 1) then
        if (this%models(i) == 1) then
          ! w_constant model
          w_de(i) = this%w0
        else if (this%models(i) == 2) then
          ! w0wa with PPF - allows phantom divide crossing
          w_de(i) = this%w0 + this%wa * (1-a)
        else if (this%models(i) == 3) then
          ! Binned w model
          z = 1._dl/a - 1._dl
          if (z < this%z1) then
            w_de(i) = this%w0
          else if (z < this%z2) then
            w_de(i) = this%w1
          else if (z < this%z3) then
            w_de(i) = this%w2
          else
            w_de(i) = this%w3
          end if
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 1"
        end if

      else if (i == 2) then
        if (this%models(i) == 1) then
          ! Axion effective fluid from arXiv:1806.10608
          zc = this%zc
          ac = 1._dl/(1._dl + zc)
          fde_zc = this%fde_zc
          wn = this%wn
          w_de(i) = (1._dl + wn)/(1._dl + (ac/a)**(3._dl*(1._dl + wn))) - 1._dl
        else if (this%models(i) == 2) then
          ! Scalar field EDE
          call this%Field_Values_At_a(a, phi, phidot)
          K = phidot**2/(2*a**2)
          V = this%Vofphi(phi, 0)
          w_de(i) = (K - V)/(K + V)
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 2"
        end if
      else 
        stop "[Multifluid DE] Invalid dark energy fluid"
      end if
    end do
  end function TMultiFluidDE_w_de

  function TMultiFluidDE_grho_de(this, a) result(grho_de)
    class(TMultiFluidDE) :: this
    real(dl), dimension(max_num_of_fluids) :: grho_de
    real(dl), intent(in) :: a
    integer :: i
    real(dl) :: w0, wa, zc, ac, fde_zc, oma_zc, wn, z, phi, phidot, grho_late_today, grho_ede_today

    ! Returns 8*pi*G*rho_de, no scale factors
    ! grhode_today = grhocrit * (1 - Omega_baryons - Omega_cdm - Omega_neutrinos - Omega_photons)
    grho_de = 0
    do i = 1, this%num_of_components
      if (i == 1) then
        if (this%num_of_components == 2 .and. this%models(2) == 1 .and. this%zc < 100) then ! Discounting EDE contribution today from total DE energy
          zc = this%zc
          fde_zc = this%fde_zc
          wn = this%wn
          oma_zc = fde_zc * this%grhonode_zc / (grhocrit * (1 - fde_zc))
          grho_ede_today = grhocrit * (2._dl * oma_zc)/(1 + (1+zc)**(3*(1+wn)))
          grho_late_today = grhode_today - grho_ede_today
        else ! Don't want to do these calculations if EDE is indeed negligible
          grho_late_today = grhode_today
        end if
        if (this%models(i) == 1) then
          ! w_constant model
          !w0 = this%de_params(max_num_of_params*(i-1) + 1)
          w0 = this%w0
          grho_de(i) = grho_late_today * a**(-3*(1 + w0))
        else if (this%models(i) == 2) then
          w0 = this%w0
          wa = this%wa
          grho_de(i) = grho_late_today * a**(-3*(1 + w0 + wa)) * exp(-3 * wa * (1-a))
        else if (this%models(i) == 3) then
          ! Binned w model
          z = 1._dl/a - 1
          if (z < this%z1) then
            grho_de(i) = grho_late_today * a**(-3*(1 + this%w0))
          else if (z < this%z2) then
            grho_de(i) = grho_late_today * (1+this%z1)**(3*(1 + this%w0)) * (a * (1+this%z1))**(-3*(1+this%w1))
          else if (z < this%z3) then
            grho_de(i) = grho_late_today * (1+this%z1)**(3*(1 + this%w0)) * ((1+this%z2)/(1+this%z1))**(3*(1+this%w1)) * (a * (1+this%z2))**(-3*(1 + this%w2))
          else
            grho_de(i) = grho_late_today * (1+this%z1)**(3*(1 + this%w0)) * ((1+this%z2)/(1+this%z1))**(3*(1+this%w1)) * ((1+this%z3) / (1+this%z2))**(3*(1 + this%w2)) * (a * (1+this%z3))**(-3*(1 + this%w3))
          end if
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 1"
        end if

      else if (i == 2) then
        if (this%models(i) == 1) then
          ! Axion effective fluid from arXiv:1806.10608
          z = 1/a - 1
          !zc = this%de_params(max_num_of_params*(i-1) + 1)
          !fde_zc = this%de_params(max_num_of_params*(i-1) + 2)
          !wn = this%de_params(max_num_of_params*(i-1) + 3)
          zc = this%zc
          fde_zc = this%fde_zc
          wn = this%wn
          ac = 1._dl/(1._dl + zc)
          oma_zc = fde_zc * this%grhonode_zc / (grhocrit * (1 - fde_zc))
          grho_de(i) = grhocrit * (2._dl * oma_zc)/(1 + ( (1+zc)/(1+z) )**(3*(1+wn)))
          ! grho_ede(a = 1) = grhocrit * (2._dl * oma_zc)/(1 + (1+zc)**(3*(1+wn)))
        else if (this%models(i) == 2) then
          ! Scalar field EDE with either rock 'n' roll potential or Axion
          if (.not. this%did_scalar_field_init) then
            ! To evolve the scalar field, I need to call the energy density of the first fluid.
            ! This ensures that I can call this function without any errors
            grho_de(i) = 0
          else
            call this%Field_Values_At_a(a, phi, phidot)
            grho_de(i) = phidot**2/(2*a**2) + this%Vofphi(phi, 0)
          end if

        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 2"
        end if
      
      else 
        stop "[Multifluid DE] Invalid dark energy fluid"
      end if
    end do
  end function TMultiFluidDE_grho_de

  subroutine TMultiFluidDE_density(this, grhov, a, grhov_t, w)
    ! Get grhov_t = 8*pi*G*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TMultiFluidDE), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), dimension(max_num_of_fluids), intent(out) :: grhov_t
    real(dl), optional, dimension(max_num_of_fluids), intent(out) :: w
    integer :: i

    if (a > 1e-10) then
      grhov_t = this%grho_de(a) * a**2
    else
      grhov_t = 0
    end if
    if (present(w)) then
      w = this%w_de(a)
    end if
  end subroutine TMultiFluidDE_density

  ! --------------------- Perturbation Equations ---------------------

  subroutine TMultiFluidDE_PerturbedStressEnergy(this, dgrhoe, dgqe, &
    a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TMultiFluidDE), intent(inout) :: this
    real(dl), dimension(max_num_of_fluids), intent(out) :: dgrhoe, dgqe
    real(dl), dimension(max_num_of_fluids), intent(in) :: grhov_t
    real(dl), dimension(max_num_of_fluids), intent(in) :: w
    real(dl), intent(in) ::  a, dgq, dgrho, grho, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: phi, phidot, delta_phi, delta_phi_prime
    integer :: i, delta_index

    ! TODO clean variables
    dgrhoe = 0
    dgqe = 0
    do i=1, this%num_of_components
      delta_index = w_ix + 2*(i - 1)
      if (i == 1 .and. (this%models(1) == 2) .or. this%models(i) == 3) then ! w0wa requires PPF
        call PPF_Perturbations(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, &
        etak, adotoa, k, kf1, ay, ayprime, w_ix)
      else
        if (this%models(1) == 2) then
          delta_index = delta_index - 1 ! PPF has only one perturbation equation, so I need to decrease one index
        end if

        if (i == 2 .and. this%models(i)==2) then
          ! Scalar field energy perturbations
          ! Check equation (A6) in https://arxiv.org/pdf/astro-ph/9801234.pdf
          call this%Field_Values_At_a(a, phi, phidot)
          delta_phi = ay(delta_index)
          delta_phi_prime = ay(delta_index + 1)
          dgrhoe(i) = phidot * delta_phi_prime + delta_phi * a**2 * this%Vofphi(phi,1)
          dgqe(i) = k * phidot * delta_phi
        else
          ! Usual fluid definitions
          dgrhoe(i) = ay(delta_index) * grhov_t(i)
          dgqe(i) = ay(delta_index + 1) * grhov_t(i) * (1 + w(i))
        end if
      end if
    end do
  end subroutine TMultiFluidDE_PerturbedStressEnergy

  function w_derivative(this, a, i) result(dwda)
    ! Returns dw/da, necessary for calculating the adiabatic sound speed
    ! Check equations 22, 23, 24 in https://arxiv.org/pdf/1806.10608.pdf
    class(TMultiFluidDE), intent(in) :: this
    integer, intent(in) :: i
    real(dl), intent(in) :: a
    real(dl) :: dwda
    real(dl) :: ac, zc, fde_zc
    real(dl) :: wn

    select case (i)
      case(1)
        if (this%models(i) == 1) then
          ! w_constant model
          dwda = 0
        else
          ! w0wa and binW models use PPF equations, so we don't need to calculate this!
          stop "[Multifluid DE] Invalid dark energy model for fluid 1"
        end if

      case(2)
        if (this%models(i) == 1) then
          ! Axion effective fluid from arXiv:1806.10608
          zc = this%zc
          wn = this%wn
          ac = 1._dl/(1._dl + zc)
          fde_zc = this%fde_zc
          dwda = 3 * ac * (1+wn)**2 * (ac/a)**(-1 + 3*(1+wn)) &
          / (a**2 * (1 + (ac/a)**(3*(1+wn)) )**2) 
        else if (this%models(i)==2) then
          stop "[Multifluid DE] Scalar Field DE doesn't need to invoke this subroutine. Something is wrong"
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 2"
        end if
      case default
        stop "[Multifluid DE] Invalid dark energy fluid"
    end select
  end function w_derivative

  subroutine TMultiFluidDE_PerturbationEvolve(this, ayprime, w, w_ix, &
    a, adotoa, k, z, y)
    class(TMultiFluidDE), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: y(:)
    real(dl), dimension(max_num_of_fluids), intent(in) :: w
    real(dl), intent(in) :: a, adotoa, k, z
    integer, intent(in) :: w_ix
    real(dl) :: Hv3_over_k, v, delta, fac
    real(dl) :: cs2, phi, phidot, delta_phi, delta_phi_prime
    integer :: i, delta_index

    ! Fluid equations in synchronous gauge https://arxiv.org/pdf/1806.10608.pdf
    ! Assuming cs_2 = 1 for all fluids
    ! PPF has only one equation for Gamma and it is initialized in the PerturbedStressEnergy subroutine
    do i=1, this%num_of_components
      delta_index = w_ix + 2*(i-1)
      if (this%models(1) == 2 .or. this%models(1) == 3) then
        delta_index = delta_index - 1 ! PPF has one less equation
        if (i == 1) then
          cycle ! The PPF evolution equations are defined in another function
        end if
      end if

      if (i == 2 .and. this%models(2) == 1 .and. this%wn < 0.9999) then
        ! Simulating sound speed dynamics for Axion fluid EDE
        fac = 2*a**(2-6*this%wn)*this%freq**2
        cs2 = (fac*(this%n-1) + k**2)/(fac*(this%n+1) + k**2)
      else
        cs2 = 1._dl
      end if

      
      if (i == 2 .and. this%models(i)==2) then
        ! Scalar field perturbation equations
        ! Check equation (A3) in https://arxiv.org/pdf/astro-ph/9801234.pdf
        call this%Field_Values_At_a(a,phi,phidot) ! wasting time calling this again..
        delta_phi = y(delta_index)
        delta_phi_prime = y(delta_index + 1)
        
        ayprime(delta_index)= delta_phi_prime
        ayprime(delta_index + 1) = - 2 * adotoa * delta_phi_prime - k * z * phidot - k**2 * delta_phi - a**2*delta_phi * this%Vofphi(phi,2)

      else
        ! Fluid perturbation equations
        delta = y(delta_index)
        v = y(delta_index + 1)
        Hv3_over_k =  3 * adotoa * v / k
        ! Density perturbation evolution
        ayprime(delta_index) = - 3 * adotoa * (cs2 - w(i)) * delta &
                                  - 3 * adotoa * (1 + w(i)) * (cs2 - w(i)) * Hv3_over_k &
                                  - (1 + w(i)) * k * v - (1 + w(i)) * k * z
        ! The remaining term is the one involving the adiabatic sound speed
        ayprime(delta_index) = ayprime(delta_index) - adotoa*Hv3_over_k*a*this%w_derivative(a,i)
        ! JVR - checked the term involving 1/3 * (d(ln(1+w))/lna)
        ! Velocity evolution
        if (abs(w(i)+1) > 1e-6) then
            ayprime(delta_index + 1) = - adotoa * (1 - 3 * cs2) * v + &
                k * cs2 * delta / (1 + w(i))
        else
            ayprime(delta_index + 1) = 0
        end if
      end if
    end do
  end subroutine TMultiFluidDE_PerturbationEvolve

  subroutine TMultiFluidDE_Effective_w_wa(this, w, wa)
    class(TMultiFluidDE), intent(inout) :: this
    real(dl), intent(out) :: w, wa
    ! JVR - For w_const, effective w_wa = w
    ! Don't need effective wwa for EDE
    if (this%models(1) == 1) then
      !w = this%de_params(1)
      w = this%w0
      wa = 0
    else if (this%models(1) == 2) then
      w = this%w0
      wa = this%wa
    else if (this%models(1) == 3) then
      w = this%w0
      wa = 0
    else
      stop "[Multifluid DE] Effective w0 wa not implemented for this DE model"
    end if
  end subroutine TMultiFluidDE_Effective_w_wa

  subroutine TMultiFluidDE_PerturbationInitial(this, y, a, tau, k)
    class(TMultiFluidDE), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer :: delta_index
    !Get initial values for perturbations at a (or tau)
    !For standard adiabatic perturbations can usually just set to zero to good accuracy

    y = 0

  end subroutine TMultiFluidDE_PerturbationInitial

  ! --------------------- PPF specific subroutines ---------------------

  subroutine PPF_Perturbations(this, dgrhoe, dgqe, &
    a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TMultiFluidDE), intent(inout) :: this
    real(dl), dimension(max_num_of_fluids), intent(out) :: dgrhoe, dgqe
    real(dl), dimension(max_num_of_fluids), intent(in) :: grhov_t
    real(dl), dimension(max_num_of_fluids), intent(in) :: w
    real(dl), intent(in) ::  a, dgq, dgrho, grho, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    integer :: i, delta_index
    real(dl) :: grhoT, vT, k2, sigma, S_Gamma, ckH, Gamma, Gammadot, Fa, c_Gamma_ppf

    k2=k**2
    grhoT = grho - grhov_t(1)
    vT = dgq / (grhoT + gpres_noDE)
    Gamma = ay(w_ix)
    c_Gamma_ppf = 0.4_dl

    !sigma for ppf
    sigma = (etak + (dgrho + 3 * adotoa / k * dgq) / 2._dl / k) / kf1 - &
        k * Gamma
    sigma = sigma / adotoa

    S_Gamma = grhov_t(1) * (1 + w(1)) * (vT + sigma) * k / adotoa / 2._dl / k2
    ckH = c_Gamma_ppf * k / adotoa

    if (ckH * ckH .gt. 3.d1) then ! ckH^2 > 30 ?????????
        Gamma = 0
        Gammadot = 0.d0
    else
        Gammadot = S_Gamma / (1 + ckH * ckH) - Gamma - ckH * ckH * Gamma
        Gammadot = Gammadot * adotoa
    endif
    ayprime(w_ix) = Gammadot !Set this here, and don't use PerturbationEvolve

    Fa = 1 + 3 * (grhoT + gpres_noDE) / 2._dl / k2 / kf1
    dgqe(1) = S_Gamma - Gammadot / adotoa - Gamma
    dgqe(1) = -dgqe(1) / Fa * 2._dl * k * adotoa + vT * grhov_t(1) * (1 + w(1))
    dgrhoe(1) = -2 * k2 * kf1 * Gamma - 3 / k * adotoa * dgqe(1)
  end subroutine PPF_Perturbations

  function TMultiFluidDE_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    ! [PPF] Get derivative of anisotropic stress
    class(TMultifluidDE), intent(in) :: this
    real(dl), intent(in) :: dgrhoe(max_num_of_fluids), dgqe(max_num_of_fluids), grho, gpres, w(max_num_of_fluids), grhok, adotoa, &
        k, grhov_t(max_num_of_fluids), z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot, hdotoh

    if (this%models(1) /= 2) then
      ppiedot = 0
    else
      hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / 6._dl / adotoa
      ppiedot = 3._dl * dgrhoe(1) + dgqe(1) * &
            (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
            grhov_t(1) * (1 + w(1)) * k * z / adotoa - 2._dl * k2 * Kf1 * &
            (yprime(w_ix) / adotoa - 2._dl * y(w_ix))
      ppiedot = ppiedot * adotoa / Kf1
    end if
    
  end function TMultiFluidDE_diff_rhopi_Add_Term

 ! --------------------- Technical Subroutines ---------------------

  subroutine TMultiFluidDE_Init(this, State)
    use classes
    use results
    class(TMultiFluidDE), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: ac, grho_rad, grho_matter, xc, F, p, mu, n, zeq

    ! JVR - If we need any information from other components, Init is the place to get it

    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 2 * this%num_of_components
    if (this%models(1) == 2) then
      this%num_perturb_equations = this%num_perturb_equations - 1
    end if

    ! Have access to quantities from other components by accessing this%State
    select type(State)
    class is (CAMBdata)
        this%State => State
    end select

    ac = 1._dl/(1+this%zc)

    select type (State)
    type is (CAMBdata)
      Omega_de = State%Omega_de
      grhode_today = State%grhov
      grhocrit = State%grhocrit
      this%grhonode_zc = State%grho_no_de(ac) / ac**4
      grho_rad = (kappa/c**2*4*sigma_boltz/c**3*State%CP%tcmb**4*Mpc**2*(1+3.046*7._dl/8*(4._dl/11)**(4._dl/3)))
      grho_matter = (State%grhoc + State%grhob) * ac
      zeq = (State%CP%ombh2+State%CP%omch2)/((State%grhog+State%grhornomass)/grhocrit) / (State%CP%H0/100)**2
    end select

    ! If wn = 1 then cs2 = 1
    ! If not, effective sound speed is dynamical
    ! Check equations in https://arxiv.org/pdf/1806.10608.pdf
    if (this%models(2)==1 .and. this%wn < 0.9999) then
      ! Get (very) approximate result for sound speed parameter; arXiv:1806.10608  Eq 30 (but mu may not exactly agree with what they used)
      n = nint((1+this%wn)/(1-this%wn))
      ac = 1._dl / (1._dl + this%zc)
      !Assume radiation or matter domination, standard neutrino model; H0 factors cancel
      F=7./8
      if (this%zc > zeq) then
        p = 1./2
        xc = ac**2 * p / sqrt(grho_rad/3)
      else
        p = 2./3
        xc = ac**2 * p / sqrt(grho_matter*ac / 3)
      end if

      mu = 1/xc*(1-cos(this%theta_i))**((1-n)/2.)*sqrt((1-F)*(6*p+2)*this%theta_i/n/sin(this%theta_i)) ! eq 28
      this%freq =  mu*(1-cos(this%theta_i))**((n-1)/2.)* &
          sqrt(const_pi)*Gamma((n+1)/(2.*n))/Gamma(1+0.5/n)*2.**(-(n**2+1)/(2.*n))*3.**((1./n-1)/2)*ac**(-6./(n+1)+3) &
          *( ac**(6*n/(n+1.))+1)**(0.5*(1./n-1)) ! eq 30
      this%n = n
    end if

    if (this%models(2)==2) then
      call Init_ScalarField(this)
    end if
  end subroutine TMultiFluidDE_Init

  subroutine TMultiFluidDE_ReadParams(this, Ini)
    use IniObjects
    class(TMultiFluidDE) :: this
    class(TIniFile), intent(in) :: Ini
    integer :: i, j
    character(len=30) :: read_model

    do i = 1, this%num_of_components
      write(read_model, "(A12,I1)") "model_fluid_", i
      this%models(i) = Ini%Read_Int(trim(read_model), 1)
    end do

    do i = 1, this%num_of_components
      do j = 1, max_num_of_params
        write(read_model, "(A12,I1)") "de_fluid_", i, "_param_", j
        this%models(i) = Ini%Read_Int(trim(read_model), 0)
      end do
    end do
  end subroutine TMultiFluidDE_ReadParams

  function TMultiFluidDE_PythonClass()
    character(LEN=:), allocatable :: TMultiFluidDE_PythonClass
    TMultiFluidDE_PythonClass = 'MultiFluidDE'
  end function TMultiFluidDE_PythonClass

  subroutine TMultiFluidDE_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TMultiFluidDE), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
  end subroutine TMultiFluidDE_SelfPointer


  ! --------------------- Scalar Field subroutines ---------------------

  subroutine KG_Equations_Log(this,num,loga,y,yprime)
    ! Evolve the background equation in terms of loga.
    ! Variables are phi=y(1), a^2 phi' = y(2)
    ! Assume otherwise standard background components
    class(TMultiFluidDE) :: this
    integer num
    real(dl) y(num),yprime(num)
    real(dl) loga, a

    a = exp(loga)
    call this%KG_Equations(num, a, y, yprime)
    yprime = yprime*a
  end subroutine KG_Equations_Log

  subroutine KG_Equations(this,num,a,y,yprime)
      ! Evolve the field according to the Klein-Gordon equation in terms of a. 
      ! Variables are phi=y(1), a^2 phi' = y(2)
      ! Assume otherwise standard background components
      class(TMultiFluidDE) :: this
      integer num
      real(dl) y(num),yprime(num)
      real(dl) a, a2, tot
      real(dl) phi, grhode, phidot, adot, grho_other_de_fluids(max_num_of_fluids)

      a2=a**2
      phi = y(1)
      phidot = y(2)/a2
      grho_other_de_fluids = this%grho_de(a)
      grhode = a2 * (0.5d0*phidot**2 + a2 * this%Vofphi(phi,0)) + grho_other_de_fluids(1)*a**4
      tot = this%state%grho_no_de(a) + grhode

      adot = sqrt(tot/3.0d0)
      yprime(1)=phidot/adot !d phi /d a
      yprime(2)= -a2**2*this%Vofphi(phi,1)/adot
  end subroutine KG_Equations

  subroutine Field_Values_At_a(this,a,aphi,aphidot)
      class(TMultiFluidDE) :: this
      ! Do interpolation for background phi and phidot at a (precomputed in Init)
      real(dl) a, aphi, aphidot
      real(dl) a0,b0,ho2o6,delta,da
      integer ix

      if (a >= 0.9999999d0) then
          aphi= this%phi_a(this%npoints_linear+this%npoints_log)
          aphidot= this%phidot_a(this%npoints_linear+this%npoints_log)
          return
      elseif (a < this%astart) then
          aphi = this%phi_a(1)
          aphidot = 0
          return
      elseif (a > this%max_a_log) then
          delta= a-this%max_a_log
          ix = this%npoints_log + int(delta/this%da)
      else
          delta= log(a)-this%log_astart
          ix = int(delta/this%dloga)+1
      end if
      da = this%sampled_a(ix+1) - this%sampled_a(ix)
      a0 = (this%sampled_a(ix+1) - a)/da
      b0 = 1 - a0
      ho2o6 = da**2/6._dl
      aphi=b0*this%phi_a(ix+1) + a0*(this%phi_a(ix)-b0*((a0+1)*this%ddphi_a(ix)+(2-a0)*this%ddphi_a(ix+1))*ho2o6)
      aphidot=b0*this%phidot_a(ix+1) + a0*(this%phidot_a(ix)-b0*((a0+1)*this%ddphidot_a(ix)+(2-a0)*this%ddphidot_a(ix+1))*ho2o6)
  end subroutine Field_Values_At_a

  function Vofphi(this, phi, deriv)
      !The input variable phi is sqrt(8*Pi*G)*psi
      !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
      !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
      class(TMultiFluidDE), intent(in) :: this
      real(dl), intent(in) :: phi
      integer, intent(in) :: deriv
      real(dl) :: Vofphi
      real(dl) :: theta, costheta, V0

      ! JVR - Recipe for implementing potentials:
      ! 1. Write your potential and its derivatives in natural units (assuming all parameters
      ! are in natural units of Mpl)
      ! 2. Multiply the result by this units factor
      
      ! Assume f = sqrt(kappa)*f_theory = f_theory/M_pl
      ! m = m_theory/M_Pl
      select case (this%which_potential)
          case(1) ! Rock 'n' Roll potential, already including the cosmological constant
              V0 = this%V0
              if (deriv == 0) then
                  Vofphi = units * V0 * phi**(2*this%n)
              else if (deriv == 1) then
                  Vofphi = units * V0 * (2 * this%n) * phi**(2*this%n - 1)
              else if (deriv == 2) then
                  Vofphi = units * V0 * (2 * this%n) * (2 * this%n - 1) * phi**(2*this%n - 2)
              end if
          
          case(2) ! AxionEDE
              theta = phi/this%f
              if (deriv==0) then
                Vofphi = units * this%m**2 * this%f**2 * (1 - cos(theta))**this%n
              else if (deriv ==1) then
                Vofphi = units * this%m**2 * this%f*this%n*(1 - cos(theta))**(this%n-1)*sin(theta)
              else if (deriv ==2) then
                costheta = cos(theta)
                Vofphi = units * this%m**2 * this%n * (1 - costheta)**(this%n-1) * (this%n * (1 + costheta) -1)
              end if
      end select
  end function Vofphi

  subroutine Init_ScalarField(this)
      use Powell
      class(TMultiFluidDE), intent(inout) :: this
      real(dl) aend, afrom
      integer, parameter ::  NumEqs=2
      real(dl) c(24),w(NumEqs,9), y(NumEqs)
      integer ind, i, ix
      real(dl), parameter :: splZero = 0._dl
      real(dl) lastsign, da_osc, last_a, a_c
      real(dl) initial_phi, initial_phidot, a2
      real(dl), dimension(:), allocatable :: sampled_a, phi_a, phidot_a, fde
      integer npoints, tot_points, max_ix
      logical has_peak
      real(dl) fzero, xzero
      integer iflag, iter
      Type(TTimer) :: Timer
      Type(TNEWUOA) :: Minimize
      real(dl) log_params(2), param_min(2), param_max(2)
      real(dl) :: fmatter, frad, grhorad, grhomatter, grhoquint, w_phi, quint_to_lmbd_frac, grho_other_de_fluids(max_num_of_fluids)

      ! JVR - The Init subroutine is used to initialize all variables in the class that are
      ! needed to compute the background quantities and the perturbation equations
      ! Thus, we need to evolve numerically the field, since rho_phi(a) = phi'(a)^2 / 2a^2 + V(phi(a))
      ! All other components have been initialized at this point

      this%log_astart = log(this%astart)

      if (this%use_zc) then
          call zc_fde_adjust(this)
      end if

      ! Setting initial field value
      if (this%which_potential == 1) then
          initial_phi = this%initial_phi
      else
          initial_phi = this%theta_i*this%f
      end if

      ! Use log spacing in a up to max_a_log, then linear. Switch where step matches
      this%dloga = (-this%log_astart)/(this%npoints-1)
      this%max_a_log = 1.d0/this%npoints/(exp(this%dloga)-1)
      npoints = (log(this%max_a_log)-this%log_astart)/this%dloga + 1

      if (allocated(this%phi_a)) then
          deallocate(this%phi_a,this%phidot_a)
          deallocate(this%ddphi_a,this%ddphidot_a, this%sampled_a)
      end if
      allocate(phi_a(npoints),phidot_a(npoints), sampled_a(npoints), fde(npoints))

      y(1) = initial_phi
      initial_phidot = 0
      y(2) = initial_phidot*this%astart**2

      phi_a(1) = y(1)
      phidot_a(1) = y(2)/this%astart**2
      sampled_a(1) = this%astart
      da_osc = 1
      last_a = this%astart
      max_ix = 0


      ! ----------------- Background Integration -----------------------

      ind=1
      afrom=this%log_astart
      do i=1, npoints-1
          aend = this%log_astart + this%dloga*i
          ix = i+1
          sampled_a(ix)=exp(aend)
          a2 = sampled_a(ix)**2
          grho_other_de_fluids = this%grho_de(sampled_a(ix))
          quint_to_lmbd_frac = (phidot_a(i)**2/(2*a2) + this%Vofphi(phi_a(i), 0))/grho_other_de_fluids(1)
          if (quint_to_lmbd_frac < 1e-4) then
              phi_a(ix) = 0
              phidot_a(ix) = 0
              cycle
          end if 
          call dverk(this,NumEqs,KG_Equations_Log,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
          if (.not. this%check_error(exp(afrom), exp(aend))) return
          call KG_Equations_Log(this,NumEqs,aend,y,w(:,1))
          phi_a(ix)=y(1)
          phidot_a(ix)=y(2)/a2
          
          ! probing oscillations
          !if (i==1) then
          !    lastsign = y(2)
          !elseif (y(2)*lastsign < 0) then
          !    !derivative has changed sign. Use to probe any oscillation scale:
          !    da_osc = min(da_osc, exp(aend) - last_a)
          !    last_a = exp(aend)
          !    lastsign= y(2)
          !end if

          !Define fde as ratio of early dark energy density to total
          
          !if (max_ix==0 .and. ix > 2 .and. fde(ix)< fde(ix-1)) then
          !    max_ix = ix-1
          !end if
          !if (sampled_a(ix)*(exp(this%dloga)-1)*this%min_steps_per_osc > da_osc) then
          !    !Step size getting too big to sample oscillations well
          !    exit
          !end if
      end do
      ! Do remaining steps with linear spacing in a, trying to be small enough
      this%npoints_log = ix
      this%max_a_log = sampled_a(ix)
      this%da = min(this%max_a_log *(exp(this%dloga)-1), &
          da_osc/this%min_steps_per_osc, (1- this%max_a_log)/(this%npoints-this%npoints_log))
      this%npoints_linear = int((1- this%max_a_log)/ this%da)+1
      this%da = (1- this%max_a_log)/this%npoints_linear

      tot_points = this%npoints_log+this%npoints_linear
      allocate(this%phi_a(tot_points),this%phidot_a(tot_points))
      allocate(this%ddphi_a(tot_points),this%ddphidot_a(tot_points))
      allocate(this%sampled_a(tot_points), this%fde(tot_points), this%ddfde(tot_points))
      this%sampled_a(1:ix) = sampled_a(1:ix)
      this%phi_a(1:ix) = phi_a(1:ix)
      this%phidot_a(1:ix) = phidot_a(1:ix)
      this%sampled_a(1:ix) = sampled_a(1:ix)
      this%fde(1:ix) = fde(1:ix)

      ind=1
      afrom = this%max_a_log
      do i=1, this%npoints_linear
          ix = this%npoints_log + i
          aend = this%max_a_log + this%da*i
          a2 =aend**2
          this%sampled_a(ix)=aend
          grho_other_de_fluids = this%grho_de(aend)
          quint_to_lmbd_frac = (this%phidot_a(ix-1)**2/(2*a2) + this%Vofphi(this%phi_a(ix-1), 0))/grho_other_de_fluids(1)
          if (quint_to_lmbd_frac < 1e-4) then
              this%phi_a(ix) = 0
              this%phidot_a(ix) = 0
              cycle
          end if
          call dverk(this,NumEqs,KG_Equations,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
          if (.not. this%check_error(afrom, aend)) return
          call KG_Equations(this,NumEqs,aend,y,w(:,1))
          this%phi_a(ix)=y(1)
          this%phidot_a(ix)=y(2)/a2

          this%fde(ix) = 1/((this%state%grho_no_de(aend) +  grho_other_de_fluids(1) * a2**2) &
              /(a2*(0.5d0* this%phidot_a(ix)**2 + a2*this%Vofphi(y(1),0))) + 1)
          if (max_ix==0 .and. this%fde(ix)< this%fde(ix-1)) then
              max_ix = ix-1
          end if
      end do

      call spline(this%sampled_a,this%phi_a,tot_points,splZero,splZero,this%ddphi_a)
      call spline(this%sampled_a,this%phidot_a,tot_points,splZero,splZero,this%ddphidot_a)
      call spline(this%sampled_a,this%fde,tot_points,splZero,splZero,this%ddfde)
      this%did_scalar_field_init = .true.
  end subroutine Init_ScalarField

  logical function check_error(this, afrom, aend)
      class(TMultiFluidDE) :: this
      real(dl) afrom, aend

      if (global_error_flag/=0) then
          write(*,*) 'Scalar field error integrating', afrom, aend
          write(*,*) this%n, this%f, this%m, this%theta_i
          !stop
          check_error = .false.
          return
      end if
      check_error= .true.
  end function check_error

  logical function fde_peak(this, peak, xlo, xhi, Flo, Fhi, ddFlo, ddFhi)
      class(TMultiFluidDE) :: this
      real(dl), intent(out) :: peak
      real(dl) Delta
      real(dl), intent(in) :: xlo, xhi, ddFlo, ddFhi,Flo, Fhi
      real(dl) a, b, c, fac

      !See if derivative has zero in spline interval xlo .. xhi

      Delta = xhi - xlo

      a = 0.5_dl*(ddFhi-ddFlo)/Delta
      b = (xhi*ddFlo-xlo*ddFhi)/Delta
      c = (Fhi-Flo)/Delta+ Delta/6._dl*((1-3*xhi**2/Delta**2)*ddFlo+(3*xlo**2/Delta**2-1)*ddFhi)
      fac = b**2-4*a*c
      if (fac>=0) then
          fac = sqrt(fac)
          peak = (-b + fac)/2/a
          if (peak >= xlo .and. peak <= xhi) then
              fde_peak = .true.
              return
          else
              peak = (-b - fac)/2/a
              if (peak >= xlo .and. peak <= xhi) then
                  fde_peak = .true.
                  return
              end if
          end if
      end if
      fde_peak = .false.
  end function fde_peak

  function match_zc(this, logm)
      class(TMultiFluidDE), intent(inout) :: this
      real(dl), intent(in) :: logm
      real(dl) match_zc, zc, fde_zc

      this%m = exp(logm)
      call this%calc_zc_fde(zc, fde_zc)
      match_zc = zc - this%zc

      end function match_zc

  function match_fde(this, logf)
      class(TMultiFluidDE), intent(inout) :: this
      real(dl), intent(in) :: logf
      real(dl) match_fde, zc, fde_zc

      this%f = exp(logf)
      call this%calc_zc_fde(zc, fde_zc)
      match_fde = fde_zc - this%fde_zc
  end function match_fde

  function match_fde_zc(this, x)
      class(TMultiFluidDE) :: this
      real(dl), intent(in) :: x(:)
      real(dl) match_fde_zc, zc, fde_zc

      if (this%which_potential == 2) then
          this%f = exp(x(1))
          this%m = exp(x(2))
      else
          this%V0 = exp(x(1))
          this%initial_phi = exp(x(2))
      end if
      
      call this%calc_zc_fde(zc, fde_zc)

      match_fde_zc = (log(this%fde_zc)-log(fde_zc))**2 + (log(zc)-log(this%zc))**2
      if (this%DebugLevel>1) then
          write(*,*) 'search f, m, zc, fde_zc, chi2', this%f, this%m, zc, fde_zc, match_fde_zc
      end if
  end function match_fde_zc

  subroutine calc_zc_fde(this, z_c, fde_zc)
      class(TMultiFluidDE), intent(inout) :: this
      real(dl), intent(out) :: z_c, fde_zc
      real(dl) aend, afrom
      integer, parameter ::  NumEqs=2
      real(dl) c(24),w(NumEqs,9), y(NumEqs)
      integer ind, i, ix
      real(dl), parameter :: splZero = 0._dl
      real(dl) a_c
      real(dl) initial_phi, initial_phidot, a2
      real(dl), dimension(:), allocatable :: sampled_a, fde, ddfde
      integer npoints, max_ix
      logical has_peak
      real(dl) a0, b0, da, grho_other_de_fluids(max_num_of_fluids)

      ! Get z_c and f_de(z_c) where z_c is the redshift of (first) peak of f_de (de energy fraction)
      ! Do this by forward propagating until peak, then get peak values by cubic interpolation

      if (this%which_potential == 1) then
          initial_phi = this%initial_phi
      else
          initial_phi = this%theta_i*this%f
      end if

      this%log_astart = log(this%astart)
      this%dloga = (-this%log_astart)/(this%npoints-1)

      npoints = this%npoints

      allocate(sampled_a(npoints), fde(npoints), ddfde(npoints))

      y(1) = initial_phi
      initial_phidot = 0
      y(2) = initial_phidot*this%astart**2
      sampled_a(1) = this%astart
      max_ix = 0
      ind = 1
      afrom = this%log_astart
      do i = 1, npoints-1
          aend = this%log_astart + this%dloga*i
          ix = i+1
          sampled_a(ix)=exp(aend)
          a2 = sampled_a(ix)**2
          grho_other_de_fluids = this%grho_de(sampled_a(ix))
          call dverk(this,NumEqs,KG_Equations_Log,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
          if (.not. this%check_error(exp(afrom), exp(aend))) return
          call KG_Equations_Log(this,NumEqs,aend,y,w(:,1))
          fde(ix) = 1/((this%state%grho_no_de(sampled_a(ix)) +  grho_other_de_fluids(1) * a2**2) &
              /(0.5d0*y(2)**2/a2 + a2**2*this%Vofphi(y(1),0)) + 1)
          if (max_ix == 0 .and. ix > 2 .and. fde(ix) < fde(ix-1)) then
              max_ix = ix-1
          end if
          if (max_ix/=0 .and. ix > max_ix+4) exit
      end do

      call spline(sampled_a,fde,ix,splZero,splZero,ddfde)

      has_peak = .false.
      if (max_ix >0) then
          has_peak = this%fde_peak(a_c, sampled_a(max_ix), sampled_a(max_ix+1), fde(max_ix), &
              fde(max_ix+1), ddfde(max_ix), ddfde(max_ix+1))
          if (.not. has_peak) then
              has_peak = this%fde_peak(a_c, sampled_a(max_ix-1), sampled_a(max_ix), &
                  fde(max_ix-1), fde(max_ix), ddfde(max_ix-1), ddfde(max_ix))
          end if
      end if

      if (has_peak) then
          z_c = 1/a_c-1
          ix = int((log(a_c)-this%log_astart)/this%dloga)+1
          da = sampled_a(ix+1) - sampled_a(ix)
          a0 = (sampled_a(ix+1) - a_c)/da
          b0 = 1 - a0
          fde_zc=b0*fde(ix+1) + a0*(fde(ix)-b0*((a0+1)*ddfde(ix)+(2-a0)*ddfde(ix+1))*da**2/6._dl)
      else
          write(*,*) '[EDE] calc_zc_fde: seems like there is no energy density peak'
          if (this%which_potential == 1) then
              write(*,*) '[EDE] (V0,initial_phi) = ', this%V0, this%initial_phi
          else
              write(*,*) '[EDE] (m,f) = ', this%V0, this%initial_phi
          end if
          z_c = -1
          fde_zc = 0
      end if
  end subroutine calc_zc_fde

  function fdeAta(this,a)
      class(TMultiFluidDE) :: this
      real(dl), intent(in) :: a
      real(dl) fdeAta, aphi, aphidot, a2, grho_other_de_fluids(max_num_of_fluids)

      call this%Field_Values_At_a(a, aphi, aphidot)
      a2 = a**2
      grho_other_de_fluids = this%grho_de(a)
      fdeAta = 1/((this%state%grho_no_de(a) +  grho_other_de_fluids(1) * a2**2) &
          /(a2*(0.5d0* aphidot**2 + a2*this%Vofphi(aphi,0))) + 1)
  end function fdeAta

  subroutine zc_fde_adjust(this)
      ! Find underlying parameters (V0, initial_phi) for Monomial or (m,f) for AxionEDE
      ! to give specified zc and fde_zc (peak early dark energy fraction)
      ! Input parameters are used as starting values for search, which is done by brute force
      ! (so should generalize easily, but not optimized for this specific potential)
      use Powell
      Type(TTimer) :: Timer
      Type(TNEWUOA) :: Minimize
      real(dl) fzero, xzero
      integer iflag, iter
      class(TMultiFluidDE), intent(inout) :: this
      real(dl) :: log_params(2)
      real(dl) :: ac

      ! 1 - log_params(*) will be the initial guess for the algorithm
      if (this%which_potential == 1) then
          ! log_params(1) = log(this%V0)]
          ac = 1._dl/(1 + this%zc)
          log_params(1) = log(this%V0)
          log_params(2) = log(this%initial_phi)
      else
          log_params(1) = log(this%f)
          log_params(2) = log(this%m)
      end if

      if (.false.) then
          ! Can just iterate linear optimizations when nearly orthogonal
          call Timer%Start()
          do iter = 1, 2
              call brentq(this,match_fde,log(0.01_dl),log(10._dl), 1d-3,xzero,fzero,iflag)
              if (iflag/=0) print *, 'BRENTQ FAILED f'
              this%f = exp(xzero)
              print *, 'match to m, f =', this%m, this%f, fzero
              call brentq(this,match_zc,log(1d-55),log(1d-52), 1d-3,xzero,fzero,iflag)
              if (iflag/=0) print *, 'BRENTQ FAILED m'
              this%m = exp(xzero)
              print *, 'match to m, f =', this%m, this%f, fzero
              call this%calc_zc_fde(fzero, xzero)
              print *, 'matched outputs', fzero, xzero
          end do
          call Timer%WriteTime('Timing for fitting')
      end if
      
      ! 2 - Could be useful to time the fitting process
      if (this%DebugLevel>0) call Timer%Start()
      !Minimize in log f, log m
      ! param_min(1) = log(0.001_dl)
      ! param_min(2) = log(1d-58)
      ! param_max(1) = log(1e5_dl)
      ! param_max(2) = log(1d-50)
      ! if (Minimize%BOBYQA(this, match_fde_zc, 2, 5, log_params,param_min, &
      !           param_max, 0.8_dl,1e-4_dl,this%DebugLevel,2000)) then

      ! 3 - NEWUOA is an Unconstrained Optimization Algorithm.
      ! Check its implementation in PowellMinimize.f90
      ! log_params will be updated with values that give desired zc/fde
      ! Arguments: NEWUOA(this - class instance where values are stored;
      !                   match_fde_zc - the function to be optimized;
      !                   2 - the number of variables
      !                   5 - interpolation conditions
      !                   log_params - the initial values)
      if (Minimize%NEWUOA(this, match_fde_zc, 2, 5, log_params,&
          0.8_dl,this%zcfdeprec,this%DebugLevel,500)) then

          if (Minimize%Last_bestfit > 1e-3) then
              global_error_flag = error_darkenergy
              global_error_message = '[EDE] ERROR converging solution for fde, zc'
              write(*,*) 'last-bestfit = ', Minimize%Last_bestfit
              return
          end if

          if (this%which_potential == 1) then
              this%V0 = exp(log_params(1))
              this%initial_phi = exp(log_params(2))
          else
              this%f = exp(log_params(1))
              this%m = exp(log_params(2))
          end if
          
          if (this%DebugLevel>0) then
              call this%calc_zc_fde(fzero, xzero)
              write(*,*) 'matched outputs Bobyqa zc, fde = ', fzero, xzero
          end if
      else
          global_error_flag = error_darkenergy
          global_error_message = '[EDE] ERROR finding solution for fde, zc'
          return
      end if
      
      if (this%DebugLevel>0) call Timer%WriteTime('Timing for parameter fitting')
  end subroutine zc_fde_adjust

end module MultiFluidDE