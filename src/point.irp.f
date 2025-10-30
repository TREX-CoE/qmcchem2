BEGIN_PROVIDER [ real, point, (3) ]
  implicit none
  BEGIN_DOC
  ! Coordinates of the current point
  END_DOC
  point(1) = 0.d0
  point(2) = 0.d0
  point(3) = 0.d0
END_PROVIDER


BEGIN_PROVIDER [ real, point_nucl_dist_vec, (nucl_num,3) ]
  implicit none
  BEGIN_DOC
  ! Distance vector between the current point and the nuclei
  END_DOC
  integer                        :: k
  do k=1,nucl_num
    point_nucl_dist_vec(k,1) = point(1)-nucl_coord(k,1)
    point_nucl_dist_vec(k,2) = point(2)-nucl_coord(k,2)
    point_nucl_dist_vec(k,3) = point(3)-nucl_coord(k,3)
  enddo
END_PROVIDER


BEGIN_PROVIDER [ real, point_nucl_dist, (nucl_num) ]
  implicit none
  BEGIN_DOC
  ! Distance between the current point and the nuclei
  END_DOC
  integer                        :: k,l
  do k=1,nucl_num
    point_nucl_dist(k) = point_nucl_dist_vec(k,1)*point_nucl_dist_vec(k,1) +&
        point_nucl_dist_vec(k,2)*point_nucl_dist_vec(k,2) +          &
        point_nucl_dist_vec(k,3)*point_nucl_dist_vec(k,3)
    point_nucl_dist(k) = sqrt(point_nucl_dist(k))
  enddo
END_PROVIDER

