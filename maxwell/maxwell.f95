program MaxwellModel
    implicit none
    integer :: N, i, j, k, l, m,p
    real, allocatable :: G(:), tau(:), ni(:), kl(:), kd(:)
    real ::  ni0, produkt, vsota, dt
    real ::  t, e, e1, e2, e3, e4, e5, de,de1,de2,de3,de4,de5, sigma, napetost
    real :: k1e,k2e,k3e,k4e
    real :: k1e1,k2e1,k3e1,k4e1
    real :: k1e2,k2e2,k3e2,k4e2
    real :: k1e3,k2e3,k3e3,k4e3
    real :: k1e4,k2e4,k3e4,k4e4
    real :: k1e5,k2e5,k3e5,k4e5

    N = 6
    allocate(G(N))
    allocate(tau(N))
    allocate(ni(N))
    allocate(kl(N))
    allocate(kd(N))


    ni0 = 1.0
    dt = 0.01

    G = (/0.01,0.005,0.0005,1.0,3.0,0.00069/)
    tau = (/1.1,1.2,1.3,1.4,1.5,1.1/)
    ni = tau*G
    kl = 0.0
    kd = 0.0


    do i = 1, N
        kd(1) = kd(1) + G(i)*tau(i)
    end do

    do i = 1, N-1
        do j = i+1, N
            vsota = G(i)+G(j)
            produkt = tau(i)*tau(j)
            kd(2) = kd(2) + vsota*produkt
        end do
    end do

    do i = 1, N-2
        do j = i+1, N-2
            do k = j+1, N
                vsota = G(i)+G(j)+G(k)
                produkt = tau(i)*tau(j)*tau(k)
                kd(3) = kd(3) + vsota*produkt
            end do
        end do
    end do

    do i = 1, N-3
        do j = i+1, N-2
            do k = j+1, N-1
                do l = k+1, N
                    vsota = G(i)+G(j)+G(k)+G(l)
                    produkt = tau(i)*tau(j)*tau(k)*tau(l)
                    kd(4) = kd(4) + vsota*produkt
                end do
            end do
        end do
    end do

    do i = 1, N-4
        do j = i+1, N-3
            do k = j+1, N-2
                do l = k+1, N-1
                    do m = l+1, N
                        vsota = G(i)+G(j)+G(k)+G(l)+G(m)
                        produkt = tau(i)*tau(j)*tau(k)*tau(l)*tau(m)
                        kd(5) = kd(5) + vsota*produkt
                    end do
                end do
            end do
        end do
    end do

    do i = 1, N-5
        do j = i+1, N-4
            do k = j+1, N-3
                do l = k+1, N-2
                    do m = l+1, N-1
                        do p = m+1, N
                            vsota = G(i)+G(j)+G(k)+G(l)+G(m)+G(p)
                            produkt = tau(i)*tau(j)*tau(k)*tau(l)*tau(m)*tau(p)
                            kd(6) = kd(6) + vsota*produkt
                        end do
                    end do
                end do
            end do
        end do
    end do

    print*, "kd", kd

    t = 0.0



    sigma = napetost(t)
    e  = 0.0
    e1 = 0.0
    e2 = 0.0
    e3 = 0.0
    e4 = 0.0
    e5 = 0.0

    open(15, file="rungekutta.dat")
    write(15,*) t, sigma, e, e1, e2, e3, e4, e5


    do while (t <= 120.0)

        k1e  = dt* de(e1)
        k1e1 = dt*de1(e2)
        k1e2 = dt*de2(e3)
        k1e3 = dt*de3(e4)
        k1e4 = dt*de4(e5)
        k1e5 = dt*de5(kd,t,e1,e2,e3,e4,e5)

        k2e  = dt* de(e1+0.5*k1e1)
        k2e1 = dt*de1(e2+0.5*k1e2)
        k2e2 = dt*de2(e3+0.5*k1e3)
        k2e3 = dt*de3(e4+0.5*k1e4)
        k2e4 = dt*de4(e5+0.5*k1e5)
        k2e5 = dt*de5(kd,t+0.5*dt,e1+0.5*k1e1,e2+0.5*k1e2,e3+0.5*k1e3,e4+0.5*k1e4,e5+0.5*k1e5)

        k3e  = dt* de(e1+0.5*k2e1)
        k3e1 = dt*de1(e2+0.5*k2e2)
        k3e2 = dt*de2(e3+0.5*k2e3)
        k3e3 = dt*de3(e4+0.5*k2e4)
        k3e4 = dt*de4(e5+0.5*k2e5)
        k3e5 = dt*de5(kd,t+0.5*dt,e1+0.5*k2e1,e2+0.5*k2e2,e3+0.5*k2e3,e4+0.5*k2e4,e5+0.5*k2e5)

        k4e  = dt* de(e1+k3e1)
        k4e1 = dt*de1(e2+k3e2)
        k4e2 = dt*de2(e3+k3e3)
        k4e3 = dt*de3(e4+k3e4)
        k4e4 = dt*de4(e5+k3e5)
        k4e5 = dt*de5(kd,t+dt,e1+k3e1,e2+k3e2,e3+k3e3,e4+k3e4,e5+k3e5)

        e  = e +      (k1e+2.0*k2e+2.0*k3e+k4e)/6.0
        e1 = e1 + (k1e1+2.0*k2e1+2.0*k3e1+k4e1)/6.0
        e2 = e2 + (k1e2+2.0*k2e2+2.0*k3e2+k4e2)/6.0
        e3 = e3 + (k1e3+2.0*k2e3+2.0*k3e3+k4e3)/6.0
        e4 = e4 + (k1e4+2.0*k2e4+2.0*k3e4+k4e4)/6.0
        e5 = e5 + (k1e5+2.0*k2e5+2.0*k3e5+k4e5)/6.0
        sigma = napetost(t)
        t = t + dt


        write(15,*) t, sigma, e, e1, e2, e3, e4, e5

    end do

    close(15)


    end program MaxwellModel


    ! napetost sigma
    function napetost(t)
        implicit none
        real :: napetost, t

        if (t < 40.0) then
            napetost = 10.0
        else
            napetost = 0.0
        end if

    end function

    ! 1. odvod
    function de(e1)
        implicit none
        real :: de,e1
        de = e1
    end function

    ! 2. odvod
    function de1(e2)
        implicit none
        real :: de1,e2
        de1 = e2
    end function

    ! 3. odvod
    function de2(e3)
        implicit none
        real :: de2,e3
        de2 = e3
    end function

    ! 4. odvod
    function de3(e4)
        implicit none
        real :: de3,e4
        de3 = e4
    end function

    ! 5. odvod
    function de4(e5)
        implicit none
        real :: de4,e5
        de4 = e5
    end function

    ! 6. odvod
    function de5(kd,t,e1,e2,e3,e4,e5)
        implicit none
        real :: kd(6)
        real :: de5,napetost,e1,e2,e3,e4,e5, t
        de5 = (napetost(t)-kd(5)*e5-kd(4)*e4-kd(3)*e3-kd(2)*e2-kd(1)*e1)/kd(6)
    end function

