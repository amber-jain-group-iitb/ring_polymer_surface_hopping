Module mod_iesh
!! Hammes-Schiffer, Tully, JCP 101, 4657 (1994)
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28,au2mom=1.9928191410d-24
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Potential
integer nquant,n_el
real*8 g_coup,epsilon
real*8 band_width,gama_coup
real*8 V_exothermicity,omg_B,gamma_B,temperature
real*8,allocatable :: Vc(:),e_metal(:)
real*8 beta,gamma_D,lambda_B,V_reorg,V_barrier
real*8 omg_c,omg_scaled
real*8 s01,s02,x_cr
real*8,allocatable :: mass(:),omg(:)
real*8,allocatable :: knots(:),weights(:)

!! Paramter for NO-Au system using two degree of freedom. r_N_O, r_N_Au
real*8::a00,b00,c00,d00,x0_00,y0_00
real*8::a01,b01,c01,d01,x0_01,y0_01
real*8::a11,b11,c11,d11,x0_11,y0_11,e11
real*8::conv
integer::step_no

!! Input/Output
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:)
complex*16,allocatable :: rho(:,:,:)

!! Classical
integer nclass,idistribution
real*8,allocatable :: x(:),v(:),acc(:)
real*8,allocatable :: x_old(:),v_old(:),acc_old(:)
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_hop,flag_ortho
complex*16,allocatable :: delr(:,:,:),delp(:,:,:),delacc(:,:,:)
complex*16,allocatable :: delr_old(:,:,:),delp_old(:,:,:)

!! Quantum
integer nbasis
integer,allocatable :: state(:),state_tentative(:),state_old(:)
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:),vdotd(:,:),V_k_old(:)
real*8,allocatable :: Hamil_site(:,:),Hamil_diab(:,:),delH_dels(:,:,:),delH_dels_ad(:,:,:)
real*8,allocatable :: pot(:,:),force(:,:,:),force_old(:,:,:),delF(:,:,:)
complex*16,allocatable :: ci(:,:),ci_old(:,:),sigma(:,:,:)
real*8,allocatable :: si_adiab_prev(:,:)
complex*16,allocatable :: mat(:,:),mat_adiab(:,:)
real*8,allocatable :: hop_prob(:,:),W_overlap(:,:),hop_prob_net(:)
integer ielec_hop

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,dtq,total_time,curr_time,traj_num,tim_eq
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!!ring polymer beads
real*8,allocatable::x_rp(:,:),v_rp(:,:),acc_rp(:,:),x_old_rp(:,:),v_old_rp(:,:),acc_old_rp(:,:)
real*8,allocatable::delH_delr_rp(:,:,:,:),mass_rp(:)
real*8::omg_n,beta_n,init_mom
integer::n_bead,trivial_cross,mom_num
real*8,allocatable::si_adiab_rp(:,:,:),si_adiab_rp_prev(:,:,:)

!! Misc
integer nold,cnt_rate
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2
real*8 tim_T_jk
integer,allocatable:: seed(:)
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

!!tully models
real*8::a,b,c,d,e

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)
  real*8 rnd,c_0,c_e,kt

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="iesh.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) dtq
  read(10,*) total_time
  read(10,*) iwrite
  read(10,*) nstep_write
  read(10,*) nstep_avg
  read(10,*) idistribution
  read(10,*) flag_frust
  read(10,*) flag_ortho
  read(10,*) energy_cutoff
  read(10,*) nclass
  read(10,*) nquant
  read(10,*) n_el
  read(10,*) n_bead
  read(10,*) temperature
  read(10,*) band_width
  read(10,*) gama_coup
  read(10,*) iforward
  read(10,*) icollapse
  read(10,*) ifriction
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  !nquant=nquant*nb_vib
  nbasis=nquant

  energy_cutoff=energy_cutoff*wave_to_J
  temperature=temperature*au2J/kb
  band_width=band_width*au2J
  gama_coup=gama_coup*au2J
!temperature=7894.d0
  kt=kb*temperature
  beta=1.d0/(kb*temperature)
  beta_n=beta/dfloat(n_bead)
  omg_n=n_bead/(beta*hbar)
  nsteps=nint(total_time/dtc)+1
  !lambda_D=V_reorg/4.d0
  
!-----------------------------------------------------------------  
  i=nsteps/nstep_avg+1
  allocate(pop(nquant,i),pop_surf(nquant,i),pop_amp(nquant,i))
  allocate(knots(nquant/2),weights(nquant/2))
  allocate(rho(nquant,nquant,i))
  allocate(x(nclass),v(nclass),acc(nclass))
  allocate(x_old(nclass),v_old(nclass),acc_old(nclass))
  allocate(state(n_el),state_tentative(n_el),state_old(n_el))
  allocate(mass(nclass),omg(nclass))
  allocate(delr(nquant,nquant,nclass),delp(nquant,nquant,nclass),delacc(nquant,nquant,nclass))
  allocate(delr_old(nquant,nquant,nclass),delp_old(nquant,nquant,nclass))
  allocate(si_adiab(nbasis,nquant),ci(nquant,n_el),V_k(nquant),V_k_old(nquant),sigma(nquant,nquant,n_el))
  allocate(Hamil_site(nbasis,nbasis),Hamil_diab(nbasis,nbasis),delH_dels(nbasis,nbasis,nclass),delH_dels_ad(nquant,nquant,nclass))
  allocate(pot(nquant,nquant),force(nquant,nquant,nclass),force_old(nquant,nquant,nclass),delf(nquant,nquant,nclass))
  allocate(mat(nbasis,nbasis),mat_adiab(nquant,nquant))
  allocate(d_ij(nquant,nquant,nclass),vdotd(nquant,nquant),hop_prob(nquant,n_el),W_overlap(nquant,nquant))
  allocate(hop_prob_net(nquant))
  allocate(ci_old(nquant,n_el),si_adiab_prev(nbasis,nquant))
  allocate(Vc(nquant),e_metal(nquant))

  allocate(x_rp(nclass,n_bead),v_rp(nclass,n_bead),acc_rp(nclass,n_bead),x_old_rp(nclass,n_bead),v_old_rp(nclass,n_bead),acc_old_rp(nclass,n_bead))
  allocate(delH_delr_rp(nquant,nquant,nclass,n_bead),mass_rp(nclass),si_adiab_rp(nquant,nquant,n_bead),si_adiab_rp_prev(nquant,nquant,n_bead))
  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
        call init_cond
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j,k,n1
  real*8 t1,t2,trans_prob1,reflec_prob,trans_prob2,kt

  call files(0)

  call cpu_time(t1)

  call setup_parameters
  call initialize_averages


do mom_num=1,40

init_mom=0.25d0*(mom_num-1)+9.d0
temperature=hbar*init_mom*au2mom/kb/mass(1)
kt=kb*temperature
beta=1.d0/kt
beta_n=real(n_bead)/kt
!write(6,*)hbar*init_mom*au2mom/kb/mass(1)
!stop

reflec_prob=0.d0
trans_prob1=0.d0
trans_prob2=0.d0

  do i=1,N_traj
    traj_num=i
    call init_cond
    call evolve(nsteps)
!write(510,*)traj_num,state
if(x(1)<0.d0.and.state(1)==1) reflec_prob=reflec_prob+1.d0
if(x(1)>0.d0.and.state(1)==1) trans_prob1=trans_prob1+1.d0
if(x(1)>0.d0.and.state(1)==2) trans_prob2=trans_prob2+1.d0
!    call average_end

  enddo
  call write_average

!write(8,'(4es15.7)')init_mom,trans_prob1/real(N_traj),reflec_prob/real(N_traj),trans_prob2/real(N_traj)
write(314,'(5es15.7)')init_mom,log(init_mom**2/4000.d0),trans_prob1/real(N_traj),reflec_prob/real(N_traj),trans_prob2/real(N_traj)
enddo

  call cpu_time(t2);tim_tot=tim_tot+t2-t1
  call files(1)

end subroutine main
!---------------------------------------------------------- 

subroutine files(flag)
  implicit none
  integer,intent(in)::flag

  if(flag==0) then
    open(10,file="output")
    open(11,file="output_cl")
    open(12,file="output_state")
    open(13,file="output_hop")
    open(14,file="output_overlap")
    open(15,file="output_dec")

    open(100,file="pop.out")
    open(101,file="cnts.out")
  else
    write(10,*)
    write(10,*)"Total time=",tim_tot
    close(10);close(11);close(12);close(13);close(14);close(15)
    close(100);close(101)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine initialize_averages
  implicit none

  cnt_frust=0.d0
  cnt_collapse=0.d0
  cnt_init=0.d0
  cnt_term=0.d0
  pop=0.d0
  rho=0.d0
  pop_surf=0.d0
  pop_amp=0.d0

end subroutine initialize_averages
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i,j
  real*8 sig_x,sig_p,rnd,ak,su
  real*8 energy_0
  real*8::sum_ci,si_coeff(nquant,n_el)


  ci=0.d0
    state(1)=1
    ci(1,1)=1.d0
!initialize centroid poirtion and momentum and generate bead pos and momentum
!and equilibrate

call initial_cond

  call evaluate_variables(0)
  call evaluate_variables(1)

  delr=0.d0
  delp=0.d0

  ihop=1
  iaverage=1
  iterminate=0
  flag_terminate=0

  curr_time=0.d0
  call evaluate_variables(0)
  call evaluate_variables(1)
  call compute_mat_diab

  !! to compute the standard deviation of the energy of the trajectory
  en_avg=0.d0;ensq_avg=0.d0
  nst_av=0

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine initial_cond
implicit none
integer::i,j
real*8::rnd,sig_x,sig_p,sigma,fic_temp,fic_kt

sigma=1.414d0
call gaussian_random_number(rnd)
x(1)=(-15.d0)*au2m

v(1)=au2mom/mass(1)*init_mom!(0.25d0*mom_num+6.d0)

sig_p=dsqrt(mass_rp(1)/beta_n)
sig_x=sigma

do i=1,n_bead
call gaussian_random_number(rnd)
x_rp(1,i)=rnd*sig_x*au2m
call gaussian_random_number(rnd)
v_rp(1,i)=(1/mass_rp(1))*(rnd*sig_p)
enddo
v_rp(1,:)=v_rp(1,:)-sum(v_rp(1,:))/real(n_bead)

call equilibrate_ring_polymer
!v_rp=0.d0
x_rp(1,:)=x_rp(1,:)+x(1)
v_rp(1,:)=v_rp(1,:)+v(1)

call evaluate_variables(0)


end subroutine initial_cond

!-----------------------------------------------------------------  
subroutine equilibrate_ring_polymer
implicit none

integer::i,j,k,steps
real*8::pos(n_bead,1000)

do steps=1,1000
x_rp=x_rp+v_rp*dtc+0.5d0*acc_rp*dtc*dtc
v_rp=v_rp+0.5d0*acc_rp*dtc
!call accelaration_rp
call rp_acc
v_rp=v_rp+0.5d0*dtc*acc_rp
!write(200,'(10es14.3)')v_rp(1,:)
!write(201,'(10es14.3)')x_rp(1,:)/au2m
enddo

end subroutine
!-----------------------------------------------------------------  
subroutine rp_acc
implicit none
integer::i,j,k

do i=1,n_bead
j=i-1
k=i+1
if(j<1) j=n_bead
if(k>n_bead) k=1

acc_rp(1,i)=-omg_n**2*(2*x_rp(1,i)-x_rp(1,j)-x_rp(1,k))
enddo

end subroutine
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in) :: nsteps
  integer i,j,nstep_sm,iflag_coll,i_do_something
  real*8 t1,t2
  integer iterm

  !call cpu_time(t1)

  call write_output(1,1)
  iterm=0
step_no=0
do while(x(1)>-18.d0*au2m.and.x(1)<10.d0*au2m)
step_no=step_no+1
!    call write_output(step_no,0)
!    call average(step_no)
    call save_old_state
    call evolve_rp_classical(dtc)
    call evolve_quantum_small_dtq
    if(ihop==1)call hop
!    if(icollapse==1)call collapse(dtc,iflag_coll)
    if(flag_terminate==1) call traj_terminate(iterm)
!      if(iterm==1)exit
    curr_time=curr_time+dtc
enddo
  call write_output(1,1)

  !call cpu_time(t2)
  !tim_evolve=tim_evolve+t2-t1

end subroutine evolve
!-----------------------------------------------------------------  

subroutine average(i)
  implicit none
  integer,intent(in) :: i
  integer j,i1,j1,k,kp
  complex*16 ci_diab(nquant,n_el),rho_ad(nquant,nquant)
  complex*16 rho_tmp(nquant,nquant)
  real*8 r_avg,U(nquant,nquant),U_exc(nquant,nquant)
  real*8 pd
  integer if_reactant
  real*8 t1,t2

  !call cpu_time(t1)

  if(iwrite==1) then
    en_avg=en_avg+energy
    ensq_avg=ensq_avg+energy*energy
    nst_av=nst_av+1
  endif

  if(iaverage==1.and.(mod(i,nstep_avg)==1.or.nstep_avg==1)) then
    if(nstep_avg==1) then
      j=i
    else
      j=i/nstep_avg+1
    endif


    ci_diab=matmul(si_adiab,ci)
    pd=0.d0
    do k=1,n_el
      pd=pd+abs(ci_diab(2,k))**2
    enddo
!if(state(1)==2) pd=1.d0
    pop(1,j)=pop(1,j)+pd


  endif

  !call cpu_time(t2)
  !tim_coll=tim_coll+t2-t1

end subroutine average
!-----------------------------------------------------------------  

subroutine average_end
  implicit none

end subroutine average_end
!-----------------------------------------------------------------  

subroutine save_old_state
  implicit none

  x_old=x
  v_old=v
  acc_old=acc
  ci_old=ci
  state_old=state
  !ci2_old=ci2
  si_adiab_prev=si_adiab
  V_k_old=V_k
  force_old=force
  energy_old=energy
  delr_old=delr
  delp_old=delp

x_old_rp=x_rp
v_old_rp=v_rp
acc_old_rp=acc_rp

end subroutine save_old_state
!-----------------------------------------------------------------  

subroutine revert_state
  implicit none

  x=x_old
  v=v_old
  state=state_old
  ci=ci_old
  delr=delr_old
  delp=delp_old
  force=force_old
  !ci2=ci2_old
  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine revert_state
!-----------------------------------------------------------------  

subroutine evolve_quantum_small_dtq
  implicit none
  integer i,nstep_qm
  real*8 dtq1,dtq2
  real*8 V_k_hold(nquant),dVk_dt(nquant)
  real*8 dforce_dt(nquant,nclass)
  complex*16 ci_prev(nquant),dci_dt(nquant)

  call compute_vdotd
  dVk_dt=(V_k-V_k_old)/dtc
!  if(icollapse==1) then
!    call compute_delH_dels_ad
!  endif

!  dtq1=0.02/maxval(vdotd)
!  dtq2=0.02*hbar/maxval(V_k-sum(V_k)/real(nquant))
!  dtq=dtq1
!  if(dtq>dtq2)dtq=dtq2

!  if(dtq>dtc)dtq=dtc

  nstep_qm=nint(dtc/dtq)
!  dtq=dtc/real(nstep_qm)
  hop_prob=0.d0
  hop_prob_net=0.d0
  V_k_hold=V_k
  V_k=V_k_old
  call compute_mat_adiab

  flag_hop=0
  do i=1,nstep_qm
    call compute_hop_prob(dtq)
    if(flag_hop==0)call check_hop(i*dtq)
    call rk4(ci,dtq,dVk_dt)
    !if(icollapse==1)call rk4_decoherence(dtc)
  enddo
  !if(icollapse==1)call vv_decoherence(dtc)

  !if(icollapse==1) then
  !  call verlet_decoherence(dtc,W_overlap,V_k_old,dvk_dt)
  !endif

  do i=1,nquant
    if(hop_prob_net(i)<0.d0)hop_prob_net=0.d0
    hop_prob_net(i)=1.d0-dexp(-hop_prob_net(i))
  enddo

end subroutine evolve_quantum_small_dtq
!-----------------------------------------------------------------  

subroutine compute_hop_prob(dtq)
  implicit none
  real*8,intent(in)::dtq
  integer i,j,state_tent(n_el)
  complex*16 AA,AA_KK
  real*8 pr

  call calculate_Akj(AA_KK,state,state)
  hop_prob=0.d0
  do j=1,n_el
    do i=1,nquant
      if(.not.(any(i==state))) then
        state_tent=state
        state_tent(j)=i
        call calculate_Akj(AA,state,state_tent)
        pr=-2*real(conjg(AA))*vdotd(i,state(j))
        pr=pr*dtq/cdabs(AA_KK)
        if(pr<0.d0)pr=0.d0     !!!! CAUTION AMBER CHECK !!!!
        hop_prob(i,j)=pr
        hop_prob_net(i)=hop_prob_net(i)+pr
      endif
    enddo
  enddo
end subroutine compute_hop_prob
!-----------------------------------------------------------------  

subroutine calculate_Akj(AA,state_k,state_j)
  implicit none
  complex*16,intent(out) :: AA
  integer,intent(in) :: state_k(n_el),state_j(n_el)
  complex*16 Sij(n_el,n_el),det_k,det_j
  integer i,j
  real*8::si_pop_mat(nbasis,n_el)

  sij=0.d0
  do i=1,n_el
    do j=1,n_el
      sij(i,j)=ci(state_k(i),j)
    enddo
  enddo
  det_k=determinant(sij,n_el)


  sij=0.d0
  do i=1,n_el
    do j=1,n_el
      sij(i,j)=ci(state_j(i),j)
    enddo
  enddo
  det_j=determinant(sij,n_el)

  AA=det_k*conjg(det_j)

end subroutine calculate_Akj
!-----------------------------------------------------------------  


subroutine check_hop(tim)
  implicit none
  real*8,intent(in)::tim
  integer i,j
  real*8 rnd,pr

  call random_number(rnd)
  pr=0.d0
  flag_hop=0
outer:  do j=1,n_el
    do i=1,nquant
      if(.not.(any(i==state))) then
        pr=pr+hop_prob(i,j)                     !pr=pr+hop_prob(i,j)
        if(rnd<pr) then
          state_tentative=state
          state_tentative(j)=i
          ielec_hop=j
          flag_hop=1
          exit outer
        endif
      endif
    enddo
  enddo outer

end subroutine check_hop
!-----------------------------------------------------------------  

subroutine rk4(ci,dtq,dVk_dt)
  implicit none
  complex*16,intent(inout)::ci(nquant,n_el)
  real*8,intent(in) :: dtq,dVk_dt(nquant)
  complex*16,dimension(nquant,n_el):: k1,k2,k3,k4

  k1=matmul(mat_adiab,ci)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k2=matmul(mat_adiab,ci+0.5*dtq*k1)
  k3=matmul(mat_adiab,ci+0.5*dtq*k2)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k4=matmul(mat_adiab,ci+dtq*k3)

  ci=ci+dtq/6.d0*(k1+2*k2+2*k3+k4)


end subroutine rk4
!-----------------------------------------------------------------  

subroutine evolve_classical(dt)
  !! Velocity Verlet
  implicit none
  integer i
  real*8,intent(in) :: dt
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(nclass),delta_v(nclass),acc_sav(nclass)
  real*8 t1,t2

  !call cpu_time(t1)
  

  if(ifriction==0) then
    !! Step 1
    x=x+v*dt+0.5*acc*dt*dt
    v=v+0.5*acc*dt
    acc_old=acc
    call evaluate_variables(0)
    v=v+0.5*dt*acc
    call evaluate_variables(1)

  endif

  if(ifriction==1) then
    gama_dt=gamma_B*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)
     call stochastic_force(delta_r,delta_v,dt)
     x=x+c1*dt*v+c2*dt*dt*acc+delta_r
     acc_old=acc
     call evaluate_variables(0)
     v=c0*v+(c1-c2)*dt*acc_old+c2*dt*acc+delta_v
     call evaluate_variables(1)
  endif

!      write(500,'(5es17.9)')curr_time*1.d15,pot_en/conv,KE_en/conv
  !call cpu_time(t2);tim_ev_cl=tim_ev_cl+t2-t1

end subroutine evolve_classical
!-----------------------------------------------------------------  
subroutine evolve_rp_classical(dt)
implicit none
real*8,intent(in)::dt
integer::i


x_rp=x_rp+v_rp*dt+0.5d0*acc_rp*dt*dt
v_rp=v_rp+0.5d0*acc_rp*dt
acc_old_rp=acc_rp
call accelaration_rp
v_rp=v_rp+0.5d0*acc_rp*dt
call coordinate_change

     call evaluate_variables(0)
     call evaluate_variables(1)
end subroutine evolve_rp_classical
!-----------------------------------------------------------------  
subroutine coordinate_change
implicit none

integer::i,j

do i=1,nclass
x(i)=sum(x_rp(i,:))/real(n_bead)
v(i)=sum(v_rp(i,:))/real(n_bead)
enddo

end subroutine
!-----------------------------------------------------------------  

subroutine traj_terminate(iterm)
  implicit none
  integer,intent(out) :: iterm

  iterm=0

end subroutine traj_terminate
!-----------------------------------------------------------------  

subroutine compute_mat_diab
  implicit none
  integer i,j
  real*8 t1,t2

  !call cpu_time(t1)

  mat=0.d0
  do i=1,nbasis
    do j=1,nbasis
      mat(i,j)=-iota/hbar*sum(si_adiab(i,:)*si_adiab(j,:)*V_k(1:nquant))
    enddo
  enddo

  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1

end subroutine compute_mat_diab
!-----------------------------------------------------------------  

subroutine compute_mat_adiab
  implicit none
  integer i,j
  real*8 t1,t2
  real*8 V_avg
  
  !call cpu_time(t1)

  mat_adiab=-vdotd
  V_avg=sum(V_k)/real(nquant)
  do i=1,nquant
    mat_adiab(i,i)=mat_adiab(i,i)-iota/hbar*(V_k(i)-V_avg)
  enddo
      
  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1
  
end subroutine compute_mat_adiab
!-----------------------------------------------------------------  

subroutine hop
  implicit none
  integer ifrust

  if(flag_hop==1) then
    call velocity_adjust1(state_tentative,ifrust)
!state=state_tentative
  endif

end subroutine hop
!-----------------------------------------------------------------  
subroutine velocity_adjust1(state_tentative,ifrust)
  implicit none
  integer,intent(in)::state_tentative(n_el)
  integer,intent(out)::ifrust
  real*8 gij,gama,aa,bb,cc,discr,dp(nclass),vd,f1,f2
  integer i,j,k,kp

  k=state(ielec_hop);kp=state_tentative(ielec_hop)
  cc=V_k(state(ielec_hop))-V_k(state_tentative(ielec_hop))

  call compute_dij_2state(x,k,kp,dp)
  dp=dp/dsqrt(sum(dp*dp))
  !dp=0.d0;dp(1)=1.d0

gama=0.d0

  aa=0.d0
  bb=0.d0
  do i=1,nclass

    aa=aa+0.5/mass(i)*(dp(i)*dp(i))
    bb=bb+(v(i)*dp(i))

  enddo

  discr=bb**2+4*aa*cc
  if(discr<0.d0) then
    ifrust=1
    cnt_frust=cnt_frust+1.d0
    if(flag_frust==0)then
      gama=bb/aa
    endif
    if(flag_frust>0)gama=0.d0
  else
    ifrust=0
    if(bb>=0.d0) gama=(bb-dsqrt(discr))/(2*aa)
    if(bb<0.d0)  gama=(bb+dsqrt(discr))/(2*aa)
    state=state_tentative
    delr=0.d0
    delp=0.d0
  endif
do j=1,nclass
  do i=1,n_bead
    v_rp(j,i)=v_rp(j,i)-gama*dp(j)/mass(j)
  enddo
enddo

  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine velocity_adjust1
!-----------------------------------------------------------------  

subroutine reverse_velocity
  implicit none
  

end subroutine reverse_velocity
!-----------------------------------------------------------------  

subroutine write_output(n,nflag)
  !! nflag=0: Writes various variables as a function of time
  !! nflag=1: writes minimal useful information at the start and end of trajectory
  implicit none
  integer,intent(in)::nflag,n
  integer i
  real*8 t1,t2
  real*8 phase

  !call cpu_time(t1)

  if(nflag==0) then
    if(iwrite==1) then
      if(mod(n,nstep_write)==1.or.nstep_write==1) then
        write(10,'(4es17.7,i5)')curr_time*1.d15,energy/wave_to_J,sum(cdabs(ci)**2),temperature!,state
        write(11,'(es15.5$)')curr_time*1.d15
        write(12,'(es15.5$)')curr_time*1.d15
    !    write(13,'(5es15.5)')curr_time*1.d15,vdotd(1,2),dasin(W_overlap(1,2))/dtc,hop_prob_net(3-state),state*1.d0
    !    write(14,'(6f15.5)')curr_time*1.d15,W_overlap(1,1:2),W_overlap(2,1:2),determinant(W_overlap,nquant)
    !    write(15,'(6es15.5)')curr_time*1.d15,delr(1,1,1)*1.d10,delr(2,2,1)*1.d10
        do i=1,nclass
          write(11,'(2es15.5$)')x(i)*1.d10,v(i)
        enddo
        write(11,*)
        do i=1,nclass
          write(12,'(2i10$)')state(i)
        enddo
        write(12,*)
        !write(500,'(5es17.9)')curr_time*1.d15,pot_en/conv,KE_en/conv
        !write(501,*)curr_time*1.d15,state
      endif
    endif
  endif

  if(nflag==1) then
    if(iwrite==0)then
      write(10,'(5es15.5)')traj_num,energy/wave_to_J,sum(cdabs(ci)**2),temperature
      write(11,*) traj_num
      write(11,'(es15.5$)')curr_time*1.d15
      do i=1,nclass
        write(11,'(2es15.5$)')x(i)*1.d10,v(i)
      enddo
      write(11,*)
      write(11,*)
    endif
    if(iwrite==1) then
      write(10,*)"traj num=",traj_num
      write(10,*)"standard deviation=",dsqrt((ensq_avg-en_avg**2/dfloat(nst_av))/dfloat(nst_av))/wave_to_J
      write(10,*)"ci**2=",sum(cdabs(ci)**2)
      write(10,*);write(10,*)
      write(11,*);write(11,*)
      write(12,*);write(12,*)
      write(13,*);write(13,*)
      write(14,*);write(14,*)
      write(15,*);write(15,*)
    endif
  endif

  !call cpu_time(t2)
  !tim_wr_out=tim_wr_out+t2-t1

end subroutine write_output
!-----------------------------------------------------------------  

subroutine write_average
  !! Writes the final useful output
  implicit none
  integer i,j,i1,k
  real*8 nf,pop_el(2)

  nf=dfloat(n_traj)
  cnt_frust=cnt_frust/nf
  cnt_collapse=cnt_collapse/nf

  pop=pop/nf
  rho=rho/nf
  pop_surf=pop_surf/nf
  pop_amp=pop_amp/nf

  do i=1,nsteps/nstep_avg
    pop_el=0.d0
   ! do i1=1,2
   !   j=(i1-1)*nb_vib
   !   do k=1,nb_vib
   !     pop_el(i1)=pop_el(i1)+(rho(j+k,j+k,i))
   !   enddo
   !   write(100,'(21f15.7)')(i-1)*nstep_avg*dtc*1.d15,pop_el
   ! enddo
   write(100,'(21f15.7)')(i-1)*nstep_avg*dtc*omg_B,pop(1,i)!rho(1,1,i)
  enddo

  write(101,*) cnt_frust,cnt_collapse

end subroutine write_average
!-----------------------------------------------------------------  

subroutine evaluate_variables(flag)
  implicit none
  integer,intent(in):: flag
  integer i,j,k

  if(flag==0) then
    !! position dependant variables only
    call tise
  endif

  if(flag==1) then
    KE_en=0.d0
    do i=1,nclass
      KE_en=KE_en+0.5*mass(i)*v(i)*v(i)
    enddo

    energy=pot_en+KE_en
!write(210,*)pot_en/au2J,KE_en/au2J,energy/au2J
  endif

end subroutine evaluate_variables
!-----------------------------------------------------------------  

subroutine tise
  !! time independent schrodinger equation
  !! Output - pot_en,acc
  !! Output - V_k,d_ij
  implicit none
  integer i,j,k
  real*8 Hamil(nbasis,nbasis),ens(nbasis),vect(nbasis,nquant)
  real*8 pot_cl,acc_cl(nclass),acc_qm(nclass),dpotcl_dx(nclass)
  real*8 si_adiab_old(nquant,nbasis)
  real*8 t1,t2,si_coeff(nquant,n_el)
  real*8::pot_cl_rp,dpotcl_dx_rp(n_bead),acc_qm_rp(n_bead),acc_cl_rp(n_bead)

  !call cpu_time(t1)

  call compute_potential(Hamil,delH_dels)

  Hamil_diab=Hamil
  call diag(Hamil,nbasis,ens,vect,nquant)
  V_k=ens

  do i=1,nquant
    si_adiab(:,i)=vect(:,i)
    if(sum(si_adiab(:,i)*si_adiab_prev(:,i))<0.d0)si_adiab(:,i)=-si_adiab(:,i)
    if(sum(si_adiab(:,i)*si_adiab_prev(:,i))==0.d0)then
        call trivial_crossing
    endif
  enddo

pot_en=V_k(state(1))


!  acc_qm=0.d0
!  pot_en=0.d0
!  do j=1,n_el
!    do i=1,nclass
!     acc_qm(i)=acc_qm(i)-sum(si_adiab(:,state(j))*matmul(delH_dels(:,:,i),si_adiab(:,state(j))))/mass(i)
!    enddo
!      pot_en=pot_en+V_k(state(j))
!  enddo


!  call potential_classical(pot_cl,dpotcl_dx)

!  acc_cl=-1.d0/mass*dpotcl_dx

  !V_k=pot_cl+ens
!  pot_en=pot_en+pot_cl
!  acc=acc_cl+acc_qm
!  pot_en=pot_en+pot_cl_rp
!  acc_rp=acc_cl_rp+acc_qm_rp

  !call cpu_time(t2);tim_cl=tim_cl+(t2-t1)

end subroutine tise
!-----------------------------------------------------------------  

subroutine trivial_crossing
  implicit none
  integer i,j
  real*8 si_adiab_e(nbasis,nquant)
  real*8 Hamil_e(nbasis,nbasis),ens_e(nbasis),vect_e(nbasis,nquant)
  real*8 epsilon

trivial_cross=trivial_cross+1
  epsilon=1.d-10*wave_to_J

  Hamil_e=Hamil_diab+epsilon
  call diag(Hamil_e,nbasis,ens_e,vect_e,nquant)
  do i=1,nquant
    si_adiab_e(:,i)=vect_e(:,i)
    if(sum(si_adiab_e(:,i)*si_adiab_prev(:,i))<0.d0)si_adiab_e(:,i)=-si_adiab_e(:,i)
  enddo

  do i=1,nquant
    if(sum(si_adiab_e(:,i)*si_adiab(:,i))<0.d0)si_adiab(:,i)=-si_adiab(:,i)
  enddo

end subroutine trivial_crossing

!-----------------------------------------------------------------  

subroutine compute_dij
  implicit none
  integer i,k,kp

  do k=1,nquant-1
    do kp=k+1,nquant
      do i=1,nclass
        d_ij(k,kp,i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,kp)))
      enddo
      d_ij(k,kp,:)=d_ij(k,kp,:)/(V_k(kp)-V_k(k))
      d_ij(kp,k,:)=-d_ij(k,kp,:)
    enddo
  enddo

end subroutine compute_dij
!-----------------------------------------------------------------  

subroutine compute_dij_2state(x_hop,k,kp,dp)
  implicit none
  integer,intent(in):: k,kp
  real*8,intent(in):: x_hop(nclass)
  real*8,intent(out):: dp(nclass)
  real*8 x_sav(nclass)
  integer i

  x_sav=x
  x=x_hop
  call evaluate_variables(0)

  do i=1,nclass
    dp(i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,kp)))
  enddo
  dp=dp/(V_k(kp)-V_k(k))
  x=x_sav
  call evaluate_variables(0)

end subroutine compute_dij_2state
!-----------------------------------------------------------------  

subroutine compute_vdotd
  !! T matrix computation
  implicit none
  integer i,j,k
  real*8,dimension(nquant,nquant) :: W,ci_W,si_W
  real*8 A,B,C,D,E
  real*8 Wlj,Wlk

  !Method 3
  do i=1,nquant
    do j=1,nquant
      W_overlap(i,j)=sum(si_adiab_prev(:,i)*si_adiab(:,j))
    enddo
  enddo

!  if(flag_ortho==1)call orthoganalize(W_overlap,nquant)
  call logm(W_overlap,vdotd,nquant)
  vdotd=vdotd/dtc
end subroutine compute_vdotd
!-----------------------------------------------------------------  

subroutine orthoganalize(mat,n)
  integer,intent(in)::n
  real*8,intent(inout)::mat(n,n)
  real*8 S_mat(n,n)

  S_mat=matmul(transpose(mat),mat)
  call inverse_squareroot(S_mat,n)
  mat=matmul(mat,S_mat)

end subroutine orthoganalize
!-----------------------------------------------------------------  

subroutine setup_parameters
  implicit none
  integer i
  real*8 si_diab(nbasis,2),Vb
  real*8 c_0,c_e
  real*8 omg_max,delw
  real*8 rho

conv=1000/av

mass(1)=2000.d0*au2kg
mass_rp=mass(1)

!!model 1
!a=0.01d0*au2J          !energy
!b=1.6d0/au2m         !length inverse
!c=0.005d0*au2J       !energy
!d=1.d0/(au2m**2)          !length^-2

!!model 2
a=0.1d0*au2J
b=0.28d0/(au2m**2)
c=0.015d0*au2J
d=0.06d0/(au2m**2)
e=0.05d0*au2J

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine compute_potential(H_diab,delV_dels)
  implicit none
  real*8,intent(out) :: H_diab(nquant,nquant),delV_dels(nquant,nquant,nclass)
  integer i,j,k
  real*8 h1,dh1(nclass),der_pos(nquant,nquant,nclass)
  real*8 coup,gr_pot,ex_pot

  H_diab=0.d0
  delV_dels=0.d0

call coordinate_change
!!model 1

!if(x(1)<0.d0) H_diab(1,1)=-a*(1.d0-dexp(b*x(1)))
!if(x(1)>0.d0) H_diab(1,1)=a*(1.d0-dexp(-b*x(1)))
!H_diab(1,2)=c*dexp(-d*x(1)**2)

!H_diab(2,2)=-H_diab(1,1)
!H_diab(2,1)=H_diab(1,2)

!if(x(1)<0.d0) delV_dels(1,1,1)=a*b*dexp(b*x(1))
!if(x(1)>0.d0) delV_dels(1,1,1)=a*b*dexp(-b*x(1))

!delV_dels(2,2,1)=-delV_dels(1,1,1)
!delV_dels(1,2,1)=-2.d0*c*d*x(1)*dexp(-d*x(1)**2)
!delV_dels(2,1,1)=delV_dels(1,2,1)

!!model2
H_diab(1,1)=0.d0
H_diab(2,2)=-a*dexp(-b*x(1)**2)+e
H_diab(1,2)=c*dexp(-d*x(1)**2)
H_diab(2,1)=H_diab(1,2)

delV_dels(1,1,1)=0.d0
delV_dels(2,2,1)=2*a*b*x(1)*dexp(-b*x(1)**2)
delV_dels(1,2,1)=-2*c*d*x(1)*dexp(-d*x(1)**2)
delV_dels(2,1,1)=delV_dels(1,2,1)
end subroutine compute_potential
!-----------------------------------------------------------------  

subroutine accelaration_rp
implicit none
integer::i,j,k
real*8::bead,mass_rp_new,hamil_rp(2,2,n_bead),gr_vect(1,2),ex_vect(1,2),gr_eig,ex_eig
real*8::vect(2,2),ens(2)

delH_delr_rp=0.d0
acc_rp=0.d0
bead=1.d0!real(n_bead)
mass_rp_new=mass_rp(1)

do i=1,n_bead

j=i-1
k=i+1
if(j<1) j=n_bead
if(k>n_bead) k=1

!if(x_rp(1,i)<0.d0)then
!        delH_delr_rp(1,1,1,i)=a*b*dexp(b*x_rp(1,i))/bead+mass_rp_new*omg_n**2*(2*x_rp(1,i)-x_rp(1,j)-x_rp(1,k))
!        delH_delr_rp(2,2,1,i)=-a*b*dexp(b*x_rp(1,i))/bead+mass_rp_new*omg_n**2*(2*x_rp(1,i)-x_rp(1,j)-x_rp(1,k))
!endif

!if(x_rp(1,i)>0.d0)then
!        delH_delr_rp(1,1,1,i)=a*b*dexp(-b*x_rp(1,i))/bead+mass_rp_new*omg_n**2*(2*x_rp(1,i)-x_rp(1,j)-x_rp(1,k))
!        delH_delr_rp(2,2,1,i)=-a*b*dexp(-b*x_rp(1,i))/bead+mass_rp_new*omg_n**2*(2*x_rp(1,i)-x_rp(1,j)-x_rp(1,k))
!endif
hamil_rp(1,1,i)=0.5d0*mass_rp_new*omg_n**2*((x_rp(1,i)-x_rp(1,j))**2+(x_rp(1,i)-x_rp(1,k))**2)
hamil_rp(2,2,i)=0.5d0*mass_rp_new*omg_n**2*((x_rp(1,i)-x_rp(1,j))**2+(x_rp(1,i)-x_rp(1,k))**2)-a*dexp(-b*x_rp(1,i)**2)/bead+e/bead
hamil_rp(1,2,i)=c*dexp(-d*x_rp(1,i)**2)/bead
hamil_rp(2,1,i)=hamil_rp(1,2,i)

!call diagonalization(hamil_rp(:,:,i),gr_eig,ex_eig,gr_vect,ex_vect)
!si_adiab_rp(:,1,i)=gr_vect(1,:)
!si_adiab_rp(:,2,i)=ex_vect(1,:)

call diag(hamil_rp(:,:,i),nbasis,ens,vect,nquant)
si_adiab_rp(:,:,i)=vect(:,:)

delH_delr_rp(1,1,1,i)=mass_rp_new*omg_n**2*(2*x_rp(1,i)-x_rp(1,j)-x_rp(1,k))
delH_delr_rp(2,2,1,i)=2*a*b*x_rp(1,i)*dexp(-b*x_rp(1,i)**2)/bead+mass_rp_new*omg_n**2*(2*x_rp(1,i)-x_rp(1,j)-x_rp(1,k))
enddo

!acc_rp(1,:)=-delH_delr_rp(state(1),state(1),1,:)/mass_rp_new
do i=1,n_bead
acc_rp(1,i)=-sum(si_adiab_rp(:,state(1),i)*matmul(delH_delr_rp(:,:,1,i),si_adiab_rp(:,state(1),i)))/mass_rp_new
!acc_rp(1,i)=-sum(si_adiab(:,state(1))*matmul(delH_delr_rp(:,:,1,i),si_adiab(:,state(1))))/mass_rp_new
enddo

end subroutine accelaration_rp
!-----------------------------------------------------------------  

subroutine diagonalization(mat,gr_eig,ex_eig,gr_vec,ex_vec)
        implicit none
        real*8,intent(in),dimension(2,2)::mat
        real*8,intent(out)::gr_eig,ex_eig,gr_vec(1,2),ex_vec(1,2)
        gr_eig=0.5d0*((mat(1,1)+mat(2,2))-sqrt((mat(1,1)+mat(2,2))**2-4.d0*(mat(1,1)*mat(2,2)-(mat(1,2)*mat(2,1)))))
        ex_eig=0.5d0*((mat(1,1)+mat(2,2))+sqrt((mat(1,1)+mat(2,2))**2-4.d0*(mat(1,1)*mat(2,2)-(mat(1,2)*mat(2,1)))))
        gr_vec(1,1)=mat(1,2)/sqrt(((mat(1,1)-gr_eig)**2)+((mat(1,2))**2))
        ex_vec(1,1)=mat(1,2)/sqrt(((mat(1,1)-ex_eig)**2)+((mat(1,2))**2))
        gr_vec(1,2)=(mat(1,1)-gr_eig)/sqrt(((mat(1,1)-gr_eig)**2)+((mat(1,2))**2))
        ex_vec(1,2)=(mat(1,1)-ex_eig)/sqrt(((mat(1,1)-ex_eig)**2)+((mat(1,2))**2))
end subroutine diagonalization

!-----------------------------------------------------------------  

subroutine potential_classical(pot_cl,acc_cl)
  implicit none
  real*8,intent(out) :: pot_cl,acc_cl(nclass)
  integer i
  real*8 q1,q3

pot_cl=0.d0
acc_cl=0.d0

if(x(1)<0.d0) pot_cl=-a*(1.d0-dexp(b*x(1)))
if(x(1)>0.d0) pot_cl=a*(1.d0-dexp(-b*x(1)))
end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine check_acceleration
  !! A test subroutine that compares analytical accelerations with numerical
  !accelerations
  implicit none
  integer i,nflag
  real*8 delx,en_old,acc_sav(nclass)
  real*8 q0,rnd

  delx=1.d-17
  state=1

  do i=1,nclass
    call random_number(rnd)
    x(i)=(rnd*2-1.d0)*1.d-10
  enddo

  call init_cond

  call evaluate_variables(0)
  en_old=pot_en;acc_sav=acc

  write(6,*) "delx=",delx
  write(6,*)

  do i=1,nclass
      x(i)=x(i)+delx
      call evaluate_variables(0)
      acc(i)=-(pot_en-en_old)/delx/mass(i)
      write(6,*)"Analytical acceleration =",acc_sav(i)
      write(6,*)"Numerical acceleration  =",acc(i)
      write(6,*)"Error =",(acc(i)-acc_sav(i))/acc(i)*100.d0
      write(6,*)
      x(i)=x(i)-delx
  enddo

  stop

end subroutine check_acceleration
!---------------------------------------------------------- 

subroutine stochastic_force(delr,delv,dt)
  !! stoachastic forces for langevin equation
  !! Not used for the Holstein model results 
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass),delv(nclass)!f(nclass)
  integer i
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma_B*dt

  do i=1,nclass

    sig_r=dt*dsqrt(kb*temperature/mass(i) *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
    sig_v=dsqrt(kb*temperature/mass(i)*(1-dexp(-2*gdt)))
    sig_rv=(dt*kb*temperature/mass(i)* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

    call gaussian_random_number(rnd1)
    call gaussian_random_number(rnd2)
    delr(i)=sig_r*rnd1
    delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------  

function commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  complex*16 tmp
  integer j,k

  if(iflag==0) commute=matmul(A,B)-matmul(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do j=1,nquant
      do k=1,nquant
        commute(j,k)=B(j,k)*(A(j,j)-A(k,k))
      enddo
    enddo
  endif

  if(iflag==2) then
    !! Assume A is tridiagonal, with a_ii=0, and a_ij=-a_ji (a is assumed to be d_ij)
    do j=1,nquant
      do k=1,nquant
        tmp=0.d0
        if(j<nquant) tmp=tmp+A(j,j+1)*B(j+1,k)
        if(j>1) tmp=tmp-A(j-1,j)*B(j-1,k)
        if(k>1) tmp=tmp-A(k-1,k)*B(j,k-1)
        if(k<nquant) tmp=tmp+A(k,k+1)*B(j,k+1)
      enddo
    enddo
  endif

end function commute
!-----------------------------------------------------------------  

function anti_commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 anti_commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  integer i,j

  if(iflag==0) anti_commute=matmul(A,B)+matmul(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do i=1,nquant
      do j=1,nquant
       anti_commute(i,j)=B(i,j)*(A(i,i)+A(j,j))
      enddo
    enddo
  endif

end function anti_commute
!-----------------------------------------------------------------  

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 

subroutine logm(mat,log_mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(in):: mat(n,n)
  real*8,intent(out):: log_mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
  enddo

  log_mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine logm
!-----------------------------------------------------------------  

subroutine inverse_squareroot(mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(inout):: mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=1.d0/t(i,i)**0.5d0
  enddo

  mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine inverse_squareroot
!-----------------------------------------------------------------  

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n
  integer,intent(inout) :: nold
  complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
  real*8,intent(in) :: mat(n,n)
  complex*16,intent(out) :: T(n,n)
  complex*16,allocatable,intent(inout):: cwork(:)
  real*8 rwork(n)
  complex*16 mat_c(n,n)

  integer lwork
  logical:: select
  logical bwork(n)
  integer sdim,info,AllocateStatus

  T=mat

  info=0
  sdim=0

  if(nold.ne.n .or. .not.allocated(cwork)) then
  !if(nold.ne.n) then
    lwork=-1
    if(allocated(cwork))deallocate(cwork)
    allocate(cwork(n))
    call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
    lwork=int(cwork(1))
    deallocate(cwork)
    allocate(cwork(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in schur, allocation"
    nold=n
  endif

  lwork=size(cwork)
  call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
  if(info.ne.0) then
    write(6,*) "problem in scur",info
        write(6,*)T
    stop
  endif

end subroutine schur
!---------------------------------------------------------- 

complex*16 FUNCTION determinant(matrix, n)
    !!http://web.hku.hk/~gdli/UsefulFiles/Example-Fortran-program.html
    IMPLICIT NONE
    complex*16, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    complex*16 :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                determinant= 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    determinant= l
    DO i = 1, n
        determinant= determinant* matrix(i,i)
    END DO
    
END FUNCTION determinant
!-----------------------------------------------------------------  

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!---------------------------------------------------------- 

subroutine compute_KE_matrix_dvr(KE_matrix,ngrid,delq,mass)
  !! computes KE matrix in DVR basis
  !! Appendix A of JCP 96, 1982 (1991)
  implicit none
  integer,intent(in) :: ngrid       !! size of DVR grid
  real*8,intent(inout) :: KE_matrix(ngrid,ngrid)
  real*8,intent(in) :: delq         !! step size in position of DVR basis
  real*8,intent(in) :: mass         !! mass
  integer i,j
  real*8 pi,hbar

  pi=dacos(-1.d0);hbar=1.05457266D-34

  KE_matrix=hbar**2/(2*mass*delq**2)
  do i=1,ngrid
    do j=1,ngrid
      KE_matrix(i,j)=KE_matrix(i,j)*(-1.d0)**(i-j)
      if(i==j) then
        KE_matrix(i,j)=KE_matrix(i,j)*pi**2/3.d0
      else
        KE_matrix(i,j)=KE_matrix(i,j)*2.d0/real(i-j)**2
      endif
    enddo
  enddo
end subroutine compute_KE_matrix_dvr
!---------------------------------------------------------- 

subroutine draw_pes
  implicit none
  integer i

  call init_cond

  do i=1,n_el
    state(i)=i
  enddo
  !state(n_el-1)=state(n_el-1)+2
  state(n_el)=state(n_el)+2

  do i=1,1000
    x(1)=-5d-10+20.d-10*i/999.d0
    call evaluate_variables(0)
    write(20,*) x(1)*1.d10,V_k(1:6)/au2J
    write(21,*) x(1)*1.d10,pot_en/au2J
  enddo
  stop

end subroutine draw_pes
!-----------------------------------------------------------------  
!-----------------------------------------------------------------  

End Module mod_iesh
