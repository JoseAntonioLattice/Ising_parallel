module mod_parallel

  implicit none

contains

  pure function tile_indices(i)    
    integer, intent(in)  :: i
    integer, dimension(2):: tile_indices

    tile_indices(1) = (i * (this_image() - 1) ) / num_images() + 1 
    tile_indices(2) = (i * this_image() ) / num_images()
  end function tile_indices

  function neighbors()
    integer :: neighbors(2)

    neighbors(1) = mod(this_image()-1,num_images())
    neighbors(2) = mod(this_image()+1,num_images())

    if (neighbors(1) == 0) neighbors(1) = num_images()
    if (neighbors(2) == 0) neighbors(2) = num_images()
  end function neighbors

end module mod_parallel
