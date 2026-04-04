 module output

     use parkind, only : im=>kind_im, rb=>kind_rb

     implicit none
     public:: write_out_3d

     contains

     subroutine write_out_3d(var,filename,ntime,nlev,nlat,nlon)
 
       integer(kind=im),intent(in):: ntime,nlev,nlat,nlon
       real(kind=rb),intent(in)::    var(ntime,nlev+1,nlat,nlon)
       character(len=16), intent(in)::    filename
           
       real(kind=rb) :: var_ncl(nlon,nlat,nlev+1,ntime)
       integer(kind=im) :: itime,ilev,ilat,ilon,recl_num
           
        do itime = 1,ntime
            do ilev = 1, nlev+1
                do ilat = 1, nlat
                    do ilon = 1, nlon
                        var_ncl(ilon,ilat,ilev,itime)=var(itime,ilev,ilat,ilon)
                     end do
                 end do
             end do
        end do

        open(99,file='data_output/'//trim(filename),access='sequential',form='unformatted',status='replace',recl=ntime*(nlev+1)*nlat*nlon*2)
        write(99) var_ncl
        close(99)
        return
     end

 end module
