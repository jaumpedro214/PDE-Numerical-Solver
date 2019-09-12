module linear_system_mod

contains
	
	function gauss_seidel( matrix , indp_coef ) result(xnew)
		real(8), dimension(:), intent(in) :: indp_coef
		real(8), dimension( size(indp_coef) ) :: xold, xnew, error_vector
		real(8), dimension(:,:), intent(in) :: matrix
		real(8) error, max_error, summ
		integer :: dim
		integer :: it,i,j,k
		integer :: tries

		max_error = 0.1
		tries = 10
		dim = size(indp_coef)
		
		do it = 1, tries
			error = dble(0);
			do i = 1,dim
				summ = dble(0)
				do j = 1, dim
					if(i.ne.j) then
						summ = summ + xold(j)*matrix(i,j)
					endif
				enddo

				xnew(i) = (indp_coef(i)-summ)/matrix(i,i)
				error = max(error, abs(xnew(i)-xold(i)) )
				xold(i) = xnew(i)
			enddo

			if(error < max_error) then
				exit
			endif

		enddo

	end function gauss_seidel

	subroutine triangularize(matrix, vector)
		real(8), dimension(:,:), intent(inout) :: matrix
		real(8), dimension(:), intent(inout) :: vector
			
		real(8), dimension(2,2):: mat
		real(8) a
		integer i, j, k
		integer, dimension(2) :: dim

		dim = shape(matrix)
		
		do i = 1,dim(1)
			do j = i+1,dim(2)
				a = matrix(j,i)/matrix(i,i)							
				do k = 1, dim(2)
					matrix(j,k) = matrix(j,k) - a*matrix(i,k)
				enddo
				vector(j) = vector(j) - a*vector(i)
			enddo
		enddo

		!mat = matrix
	end subroutine triangularize

	function superior_triangular_solver(matrix, vector) result(solution)
		real(8), dimension(:,:), intent(inout) :: matrix
		real(8), dimension(:), intent(inout) :: vector
		integer i, j , k, n
		real :: summ
		real(8), dimension( size(vector) ) :: solution
		n = size(vector)

		solution(n) = vector(n)/matrix(n,n)
		do k =1,n-1
			i = n - k
			summ = 0
			do j =i+1,n
				summ = summ + solution(j)*matrix(i,j)
			end do
			solution(i) = (vector(i) - summ)/matrix(i,i)
		end do

	end function superior_triangular_solver

end module linear_system_mod