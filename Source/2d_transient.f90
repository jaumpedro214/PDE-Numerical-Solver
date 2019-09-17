program upwind2d
	implicit none

	! ============================================================================================== !
	!        Finite Differences Algorithm to solve Bidimensional Partial Differencial Equations      !
	!                      		Author: João Pedro da Silva Lima                                     !
	!							joaopedrodasivalima@gmail.com										 !
	! ============================================================================================== !

	!Archives ID
	integer output_mash, output_solution, plot_info

	!Mesh info
	integer nodes;
	integer nodes_x, nodes_y
	integer time_steps

	real(8) intr_x, intr_y
	real(8) delta_t
	real(8) delta_x, delta_y

	!Constants 
	real(8) k1, k2, k3 !Diffusion
	real(8) vx, vy !Velocity
	real(8) sig	!Reaction

	!Auxiliar matrices
	integer, dimension(9,2) :: directions
	integer, dimension(:,:), allocatable :: conectivity
	
	!Solution variables
	real(8), dimension(:), allocatable :: uold, unew

	!Auxiliar variables
	integer i, j, k, t
	integer pivot_x, pivot_y
	integer i_aux, j_aux
	real(8) A, B, C, D, E, F

	!Neuman conditon
	logical neuman_on_y

	!Direction vectors
	directions(1,1) = 0
	directions(1,2) = 0
	directions(2,1) = -1
	directions(2,2) = -1
	directions(3,1) = -1
	directions(3,2) = 0
	directions(4,1) = -1
	directions(4,2) = 1
	directions(5,1) = 0
	directions(5,2) = 1
	directions(6,1) = 1
	directions(6,2) = 1
	directions(7,1) = 1
	directions(7,2) = 0
	directions(8,1) = 1
	directions(8,2) = -1
	directions(9,1) = 0
	directions(9,2) = -1

	!Archives
	output_mash = 66
	output_solution = 67
	plot_info = 70
	open(unit = output_mash, file = "output_mash.out" )
	open(unit = output_solution, file = "output_solution.out" )
	open(unit = plot_info, file = "plot_info.out" )

	!Here is where you change the parameters equations
	!Equation Parameters =====================================================!
	intr_x = 1.0000
	intr_y = 1.0000
	
	nodes_x = 50

	neuman_on_y = .true.
	
	k1 = 1.000 ! k_{1,1}
	k2 = 0.000 ! k_{1,2} + k_{2,1}
	k3 = 1.000 ! k_{2,2}
	vx = 0.000 
	vy = -5.000
	sig = 0.000 !sigma 
	
	time_steps = 300
	delta_t = 0.20*delta_x**2/max( k1 , max(k2,k3) )
	!========================================================================!

	nodes_y = nodes_x
	delta_x = intr_x/(nodes_x-1)
	delta_y = intr_y/(nodes_y-1)
	
	A = -2*k1/delta_x**2 - 2*k3/delta_y**2 + vx/delta_x + vy/delta_y + sig
	B = k1/delta_x**2 - vx/delta_x
	C = k3/delta_y**2 - vy/delta_y
	D = k1/delta_x**2
	E = k3/delta_y**2
	F = k2/(4*delta_x*delta_y)

	allocate( conectivity( nodes_x, nodes_y ) )
	allocate( uold( nodes_x*nodes_y ) )
	allocate( unew( nodes_x*nodes_y ) )

	!Conectivity
	do i = 1, nodes_x
		do j = 1, nodes_y
			conectivity(i,j) = (i-1)*nodes_x + j
		enddo
	enddo

	!Initial condition
	do i = 1, nodes_x*nodes_y
			pivot_x = i/nodes_x+1
			pivot_y = mod(i,nodes_x)
			if(pivot_y.eq.0) then
				pivot_y = nodes_y
			endif
			uold(i) = initial_c( (pivot_x-1)*delta_x, (pivot_y-1)*delta_y )
	enddo
	write(output_solution,'(F8.4)')( uold(j), j = 1, nodes_x*nodes_y )

	!Time Loop
	do t = 1, time_steps

		!Mesh Loop
		do i = 1, nodes_x*nodes_y

			pivot_x = i/nodes_x+1
			pivot_y = mod(i,nodes_x)
			if(pivot_y.eq.0) then
				pivot_y = nodes_y
			endif

			unew(i) = uold(i) - flux( (pivot_x-1)*delta_x, (pivot_y-1)*delta_y )*delta_t
			do k = 1, size( directions )

				i_aux = pivot_x + directions(k,1)
				j_aux = pivot_y + directions(k,2)
				
				!Finite difference discretizations
				if( in_mash( i_aux, j_aux ).eqv..true. ) then
					unew(i) = unew(i) + uold( conectivity(i_aux,j_aux) )*coefficient( directions(k,1), directions(k,2) )*delta_t
				endif
				
			enddo

			!Border conditions
			if( pivot_x == 1 ) then

				unew(i) = ul( (pivot_y-1)*delta_y )

			elseif( (i-1)/nodes_x == nodes_y-1 ) then

				unew(i) = ur( (pivot_y-1)*delta_y )

			elseif( pivot_y == 1 ) then

				if( neuman_on_y.eqv..true. ) then
					unew(i) = unew( conectivity(pivot_x,pivot_y+1) ) - udn( (pivot_x-1)*delta_x )*delta_y
				else
					unew(i) = udn( (pivot_x-1)*delta_x )
				endif

			elseif( pivot_y == nodes_y ) then

				if( neuman_on_y.eqv..true. ) then
					unew(i) = unew( i-1 ) + uup( (pivot_x-1)*delta_x )*delta_y
				else
					unew(i) = uup( (pivot_x-1)*delta_x )
				endif

			endif

		enddo

		write(output_solution,'(F8.4)')( unew(j), j = 1, nodes_x*nodes_y )
		uold = unew

	enddo
	
	do i = 1, nodes_x
		write(output_mash,'(F8.4)') (i-1)*delta_x
	enddo
	do i = 1, nodes_y
		write(output_mash,'(F8.4)') (i-1)*delta_y
	enddo

	write( plot_info, '(I5)' ) nodes_x
	write( plot_info, '(I5)' ) nodes_y
	write( plot_info, '(I5)' ) time_steps

	deallocate(conectivity)
	deallocate(uold)
	deallocate(unew)

	contains

	logical function in_mash( x, y )
		
		integer, intent(in) :: x, y
		if ( (x.gt.0).and.(y.gt.0).and.(x.le.nodes_x).and.(y.le.nodes_y) ) then
			in_mash = .true.
		else
			in_mash = .false.
		endif

	end function in_mash

	!Initial problem condition
	real(8) function initial_c(x,y)
		real(8) x, y

		initial_c = 0
		if( y > 0.40 .and. y < 0.55 .and. x > 0.001 ) then
			initial_c = 20
		endif

	end function initial_c

	!Flux function
	real(8) function flux(x,y)
		real(8) x, y
		flux = 0
		
	end function flux

	!UR (U Right) is the value of u(Lx,x)
	real(8) function ur(xx)
		real(8) xx
		ur = real(10,8)
	end function ur

	!UL (U Left) is the value of u(0,x)
	real(8) function ul(xx)
		real(8) xx
		ul = real(10,8)
	end function ul

	!Uup (U Left) is the value of u(x,Ly) or ux(x,Ly)
	real(8) function uup(xx)
		real(8) xx
		uup = 0
	end function uup

	!Udn (U Left) is the value of u(x,0) or ux(x,0)
	real(8) function udn(xx)
		real(8) xx
		udn = 0
	end function udn

	!Coefficients according each direction
	real(8) function coefficient( d1, d2 )
		integer d1, d2

		if( d1.eq.(0).and.d2.eq.(0) ) then
			coefficient = A
		elseif( d1.eq.(-1).and.d2.eq.(-1) ) then
			coefficient = F
		elseif( d1.eq.(-1).and.d2.eq.(0) ) then
			coefficient = B
		elseif( d1.eq.(-1).and.d2.eq.(1) ) then
			coefficient = -F
		elseif( d1.eq.(0).and.d2.eq.(1) ) then
			coefficient = E
		elseif( d1.eq.(1).and.d2.eq.(1) ) then
			coefficient = F
		elseif( d1.eq.(1).and.d2.eq.(0) ) then
			coefficient = D
		elseif( d1.eq.(1).and.d2.eq.(-1) ) then
			coefficient = -F
		elseif( d1.eq.(0).and.d2.eq.(-1) ) then
			coefficient = C
		endif

	end function coefficient
	
	
	
end program upwind2d