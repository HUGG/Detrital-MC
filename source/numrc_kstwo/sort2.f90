	SUBROUTINE sort2(arr,slave)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
        EXTERNAL indexx_sp
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
	INTEGER(I4B) :: ndum
	INTEGER(I4B), DIMENSION(size(arr)) :: index
	ndum=assert_eq(size(arr),size(slave),'sort2')
	call indexx_sp(arr,index)
	arr=arr(index)
	slave=slave(index)
	END SUBROUTINE sort2
