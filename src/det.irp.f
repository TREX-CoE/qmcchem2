 BEGIN_PROVIDER  [ integer, det_i ]
&BEGIN_PROVIDER  [ integer, det_i_prev ]

  BEGIN_DOC
  ! Current running alpha determinant
  END_DOC
  det_i=det_alpha_order(1)
  det_i_prev=det_alpha_order(1)

END_PROVIDER

 BEGIN_PROVIDER  [ integer, det_j ]
&BEGIN_PROVIDER  [ integer, det_j_prev ]

  BEGIN_DOC
  ! Current running beta determinant
  END_DOC
  det_j=det_beta_order(1)
  det_j_prev=det_beta_order(1)

END_PROVIDER

subroutine det_update(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m(LDS) ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,n)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,n) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  if (d == 0.d0) then
    return
  endif
  select case (n)
      case default
      call det_update_general(n,LDS,m,l,S,S_inv,d)
    BEGIN_TEMPLATE
      case ($n)
        call det_update$n(n,LDS,m,l,S,S_inv,d)
    SUBST [n]
    1;;
    2;;
    3;;
    4;;
    5;;
    6;;
    7;;
    8;;
    9;;
    10;;
    11;;
    12;;
    13;;
    14;;
    15;;
    16;;
    17;;
    18;;
    19;;
    20;;
    21;;
    22;;
    23;;
    24;;
    25;;
    26;;
    27;;
    28;;
    29;;
    30;;
    31;;
    32;;
    33;;
    34;;
    35;;
    36;;
    37;;
    38;;
    39;;
    40;;
    41;;
    42;;
    43;;
    44;;
    45;;
    46;;
    47;;
    48;;
    49;;
    50;;
    51;;
    52;;
    53;;
    54;;
    55;;
    56;;
    57;;
    58;;
    59;;
    60;;
    61;;
    62;;
    63;;
    64;;
    65;;
    66;;
    67;;
    68;;
    69;;
    70;;
    71;;
    72;;
    73;;
    74;;
    75;;
    76;;
    77;;
    78;;
    79;;
    80;;
    81;;
    82;;
    83;;
    84;;
    85;;
    86;;
    87;;
    88;;
    89;;
    90;;
    91;;
    92;;
    93;;
    94;;
    95;;
    96;;
    97;;
    98;;
    99;;
    100;;
    101;;
    102;;
    103;;
    104;;
    105;;
    106;;
    107;;
    108;;
    109;;
    110;;
    111;;
    112;;
    113;;
    114;;
    115;;
    116;;
    117;;
    118;;
    119;;
    120;;
    121;;
    122;;
    123;;
    124;;
    125;;
    126;;
    127;;
    128;;
    129;;
    130;;
    131;;
    132;;
    133;;
    134;;
    135;;
    136;;
    137;;
    138;;
    139;;
    140;;
    141;;
    142;;
    143;;
    144;;
    145;;
    146;;
    147;;
    148;;
    149;;
    150;;
    END_TEMPLATE
  end select
end

subroutine det_update2(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m(2)   ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,2)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,2) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  S(1,l) = m(1)
  S(2,l) = m(2)
  S_inv(1,1) = S(1,1)
  S_inv(1,2) = S(2,1)
  S_inv(2,1) = S(1,2)
  S_inv(2,2) = S(2,2)
  call invert2(S_inv,LDS,n,d)

end

subroutine det_update1(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m(1)   ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,1)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,1) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  S(1,l) = m(1)
  S_inv(1,1) = S(1,1)
  call invert1(S_inv,LDS,n,d)

end

subroutine det_update3(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m(3)   ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,3)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,3) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  integer                        :: i
  do i=1,3
    S(i,l) = m(i)
  enddo
  do i=1,3
    S_inv(1,i) = S(i,1)
    S_inv(2,i) = S(i,2)
    S_inv(3,i) = S(i,3)
  enddo

  call invert3(S_inv,LDS,n,d)

end

subroutine det_update4(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m(4)   ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,4)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,4) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  double precision               :: u(4), z(4), w(4), lambda, d_inv
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: z, w, u
  integer                        :: i,j
  u(1) = m(1) - S(1,l)
  u(2) = m(2) - S(2,l)
  u(3) = m(3) - S(3,l)
  u(4) = m(4) - S(4,l)
  z(1) = S_inv(1,1)*u(1) + S_inv(2,1)*u(2) + S_inv(3,1)*u(3) + S_inv(4,1)*u(4)
  z(2) = S_inv(1,2)*u(1) + S_inv(2,2)*u(2) + S_inv(3,2)*u(3) + S_inv(4,2)*u(4)
  z(3) = S_inv(1,3)*u(1) + S_inv(2,3)*u(2) + S_inv(3,3)*u(3) + S_inv(4,3)*u(4)
  z(4) = S_inv(1,4)*u(1) + S_inv(2,4)*u(2) + S_inv(3,4)*u(3) + S_inv(4,4)*u(4)

  d_inv = 1.d0/d
  d = d+z(l)
  lambda = d_inv*d
  if (dabs(lambda) < 1.d-3) then
    d = 0.d0
    return
  endif

  !DIR$ VECTOR ALIGNED
  do i=1,4
    w(i) = S_inv(i,l)*d_inv
    S(i,l) = m(i)
  enddo

  do i=1,4
   !DIR$ VECTOR ALIGNED
   do j=1,4
    S_inv(j,i) = S_inv(j,i)*lambda -z(i)*w(j)
   enddo
  enddo

end

BEGIN_TEMPLATE
! Version for mod(n,4) = 0
subroutine det_update$n(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m($n)  ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,$n)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,$n) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  double precision               :: u($n), z($n), w($n), lambda, d_inv
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: z, w, u
  !DIR$ ASSUME_ALIGNED S : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED S_inv : $IRP_ALIGN
  !DIR$ ASSUME (mod(LDS,$IRP_ALIGN/8) == 0)
  !DIR$ ASSUME (LDS >= $n)
  integer                        :: i,j
  double precision :: zj, zj1, zj2, zj3

  !DIR$ NOPREFETCH
  do i=1,$n
    u(i) = m(i) - S(i,l)
  enddo

  zj = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do i=1,$n-1,4
    zj = zj + S_inv(i,l)*u(i) + S_inv(i+1,l)*u(i+1)  &
            + S_inv(i+2,l)*u(i+2) + S_inv(i+3,l)*u(i+3)
  enddo

  d_inv = 1.d0/d
  d = d+zj
  lambda = d*d_inv
  if (dabs(lambda) < 1.d-3) then
    d = 0.d0
    return
  endif

  !DIR$ VECTOR ALIGNED
  do j=1,$n,4
   zj  = 0.d0
   zj1 = 0.d0
   zj2 = 0.d0
   zj3 = 0.d0
   !DIR$ VECTOR ALIGNED
   !DIR$ NOPREFETCH
   do i=1,$n
    zj  = zj  + S_inv(i,j  )*u(i)
    zj1 = zj1 + S_inv(i,j+1)*u(i)
    zj2 = zj2 + S_inv(i,j+2)*u(i)
    zj3 = zj3 + S_inv(i,j+3)*u(i)
   enddo
   z(j  ) = zj
   z(j+1) = zj1
   z(j+2) = zj2
   z(j+3) = zj3
  enddo

  !DIR$ NOPREFETCH
  do i=1,$n
    w(i) = S_inv(i,l)*d_inv
    S(i,l) = m(i)
  enddo

  do i=1,$n,4
   zj  = z(i  )
   zj1 = z(i+1)
   zj2 = z(i+2)
   zj3 = z(i+3)
   !DIR$ VECTOR ALIGNED
   !DIR$ NOPREFETCH
   do j=1,$n
    S_inv(j,i  ) = S_inv(j,i  )*lambda - w(j)*zj
    S_inv(j,i+1) = S_inv(j,i+1)*lambda - w(j)*zj1
    S_inv(j,i+2) = S_inv(j,i+2)*lambda - w(j)*zj2
    S_inv(j,i+3) = S_inv(j,i+3)*lambda - w(j)*zj3
   enddo
  enddo

end

SUBST [ n ]
8 ;;
12 ;;
16 ;;
20 ;;
24 ;;
28 ;;
32 ;;
36 ;;
40 ;;
44 ;;
48 ;;
52 ;;
56 ;;
60 ;;
64 ;;
68 ;;
72 ;;
76 ;;
80 ;;
84 ;;
88 ;;
92 ;;
96 ;;
100 ;;
104 ;;
108 ;;
112 ;;
116 ;;
120 ;;
124 ;;
128 ;;
132 ;;
136 ;;
140 ;;
144 ;;
148 ;;

END_TEMPLATE

BEGIN_TEMPLATE
! Version for mod(n,4) = 1
subroutine det_update$n(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m($n)  ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,$n)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,$n) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  double precision               :: u($n), z($n), w($n), lambda, d_inv
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: z, w, u
  !DIR$ ASSUME_ALIGNED S : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED S_inv : $IRP_ALIGN
  !DIR$ ASSUME (mod(LDS,$IRP_ALIGN/8) == 0)
  !DIR$ ASSUME (LDS >= $n)
  integer                        :: i,j
  double precision :: zj, zj1, zj2, zj3

  do i=1,$n
    u(i) = m(i) - S(i,l)
  enddo

  zj = 0.d0
  !DIR$ NOPREFETCH
  do i=1,$n-1,4
    zj = zj + S_inv(i,l)*u(i) + S_inv(i+1,l)*u(i+1)  &
            + S_inv(i+2,l)*u(i+2) + S_inv(i+3,l)*u(i+3)
  enddo
  zj = zj + S_inv($n,l)*u($n)

  d_inv = 1.d0/d
  d = d+zj
  lambda = d*d_inv
  if (dabs(lambda) < 1.d-3) then
    d = 0.d0
    return
  endif

  !DIR$ VECTOR ALIGNED
  do j=1,$n-1,4
   zj  = 0.d0
   zj1 = 0.d0
   zj2 = 0.d0
   zj3 = 0.d0
   !DIR$ VECTOR ALIGNED
   !DIR$ NOPREFETCH
   do i=1,$n-1
    zj  = zj  + S_inv(i,j  )*u(i)
    zj1 = zj1 + S_inv(i,j+1)*u(i)
    zj2 = zj2 + S_inv(i,j+2)*u(i)
    zj3 = zj3 + S_inv(i,j+3)*u(i)
   enddo
   z(j  ) = zj  + S_inv($n,j  )*u($n)
   z(j+1) = zj1 + S_inv($n,j+1)*u($n)
   z(j+2) = zj2 + S_inv($n,j+2)*u($n)
   z(j+3) = zj3 + S_inv($n,j+3)*u($n)
  enddo

  zj  = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do i=1,$n-1
   zj = zj + S_inv(i,$n)*u(i)
  enddo
  z($n) = zj + S_inv($n,$n)*u($n)

  !DIR$ NOPREFETCH
  do i=1,$n
    w(i) = S_inv(i,l)*d_inv
    S(i,l) = m(i)
  enddo

  do i=1,$n-1,4
   zj  = z(i  )
   zj1 = z(i+1)
   zj2 = z(i+2)
   zj3 = z(i+3)
   !DIR$ VECTOR ALIGNED
   !DIR$ NOPREFETCH
   do j=1,$n-1
    S_inv(j,i  ) = S_inv(j,i  )*lambda - w(j)*zj
    S_inv(j,i+1) = S_inv(j,i+1)*lambda - w(j)*zj1
    S_inv(j,i+2) = S_inv(j,i+2)*lambda - w(j)*zj2
    S_inv(j,i+3) = S_inv(j,i+3)*lambda - w(j)*zj3
   enddo
   S_inv($n,i  ) = S_inv($n,i  )*lambda - w($n)*zj
   S_inv($n,i+1) = S_inv($n,i+1)*lambda - w($n)*zj1
   S_inv($n,i+2) = S_inv($n,i+2)*lambda - w($n)*zj2
   S_inv($n,i+3) = S_inv($n,i+3)*lambda - w($n)*zj3
  enddo

  zj = z($n)
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do i=1,$n
   S_inv(i,$n) = S_inv(i,$n)*lambda -w(i)*zj
  enddo


end

SUBST [ n ]
5 ;;
9 ;;
13 ;;
17 ;;
21 ;;
25 ;;
29 ;;
33 ;;
37 ;;
41 ;;
45 ;;
49 ;;
53 ;;
57 ;;
61 ;;
65 ;;
69 ;;
73 ;;
77 ;;
81 ;;
85 ;;
89 ;;
93 ;;
97 ;;
101 ;;
105 ;;
109 ;;
113 ;;
117 ;;
121 ;;
125 ;;
129 ;;
133 ;;
137 ;;
141 ;;
145 ;;
149 ;;

END_TEMPLATE


BEGIN_TEMPLATE
! Version for mod(n,4) = 2
subroutine det_update$n(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m($n)  ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,$n)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,$n) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  double precision               :: u($n), z($n), w($n), lambda, d_inv
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: z, w, u
  !DIR$ ASSUME_ALIGNED S : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED S_inv : $IRP_ALIGN
  !DIR$ ASSUME (mod(LDS,$IRP_ALIGN/8) == 0)
  !DIR$ ASSUME (LDS >= $n)
  integer                        :: i,j

  double precision :: zj, zj1, zj2, zj3
  !DIR$ NOPREFETCH
  do i=1,$n
    u(i) = m(i) - S(i,l)
  enddo

  zj = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do i=1,$n-2,4
    zj = zj + S_inv(i,l)*u(i) + S_inv(i+1,l)*u(i+1)  &
            + S_inv(i+2,l)*u(i+2) + S_inv(i+3,l)*u(i+3)
  enddo
  i=$n-1
  zj = zj + S_inv(i,l)*u(i) + S_inv(i+1,l)*u(i+1)

  d_inv = 1.d0/d
  d = d+zj
  lambda = d*d_inv
  if (dabs(lambda) < 1.d-3) then
    d = 0.d0
    return
  endif

  !DIR$ VECTOR ALIGNED
  do j=1,$n-2,4
   zj  = 0.d0
   zj1 = 0.d0
   zj2 = 0.d0
   zj3 = 0.d0
   !DIR$ VECTOR ALIGNED
   do i=1,$n-2
    zj  = zj  + S_inv(i,j  )*u(i)
    zj1 = zj1 + S_inv(i,j+1)*u(i)
    zj2 = zj2 + S_inv(i,j+2)*u(i)
    zj3 = zj3 + S_inv(i,j+3)*u(i)
   enddo
   z(j  ) = zj     + S_inv($n-1,j  )*u($n-1)
   z(j  ) = z(j  ) + S_inv($n,j  )*u($n)
   z(j+1) = zj1    + S_inv($n-1,j+1)*u($n-1)
   z(j+1) = z(j+1) + S_inv($n,j+1)*u($n)
   z(j+2) = zj2    + S_inv($n-1,j+2)*u($n-1)
   z(j+2) = z(j+2) + S_inv($n,j+2)*u($n)
   z(j+3) = zj3    + S_inv($n-1,j+3)*u($n-1)
   z(j+3) = z(j+3) + S_inv($n,j+3)*u($n)
  enddo

  j=$n-1
  zj  = 0.d0
  zj1 = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do i=1,$n-2
   zj  = zj  + S_inv(i,j  )*u(i)
   zj1 = zj1 + S_inv(i,j+1)*u(i)
  enddo
  z(j  ) = zj     + S_inv($n-1,j  )*u($n-1)
  z(j  ) = z(j  ) + S_inv($n,j  )*u($n)
  z(j+1) = zj1    + S_inv($n-1,j+1)*u($n-1)
  z(j+1) = z(j+1) + S_inv($n,j+1)*u($n)

  !DIR$ NOPREFETCH
  do i=1,$n
    w(i) = S_inv(i,l)*d_inv
    S(i,l) = m(i)
  enddo

  do i=1,$n-2,4
   zj  = z(i)
   zj1 = z(i+1)
   zj2 = z(i+2)
   zj3 = z(i+3)
   !DIR$ VECTOR ALIGNED
   do j=1,$n-2
    S_inv(j,i  ) = S_inv(j,i  )*lambda -zj *w(j)
    S_inv(j,i+1) = S_inv(j,i+1)*lambda -zj1*w(j)
    S_inv(j,i+2) = S_inv(j,i+2)*lambda -zj2*w(j)
    S_inv(j,i+3) = S_inv(j,i+3)*lambda -zj3*w(j)
   enddo
   S_inv($n-1,i  ) = S_inv($n-1,i  )*lambda -zj *w($n-1)
   S_inv($n  ,i  ) = S_inv($n  ,i  )*lambda -zj *w($n  )
   S_inv($n-1,i+1) = S_inv($n-1,i+1)*lambda -zj1*w($n-1)
   S_inv($n  ,i+1) = S_inv($n  ,i+1)*lambda -zj1*w($n  )
   S_inv($n-1,i+2) = S_inv($n-1,i+2)*lambda -zj2*w($n-1)
   S_inv($n  ,i+2) = S_inv($n  ,i+2)*lambda -zj2*w($n  )
   S_inv($n-1,i+3) = S_inv($n-1,i+3)*lambda -zj3*w($n-1)
   S_inv($n  ,i+3) = S_inv($n  ,i+3)*lambda -zj3*w($n  )
  enddo

  i=$n-1
  zj = z(i)
  zj1= z(i+1)
  !DIR$ VECTOR ALIGNED
  do j=1,$n-2
   S_inv(j,i  ) = S_inv(j,i  )*lambda -zj*w(j)
   S_inv(j,i+1) = S_inv(j,i+1)*lambda -zj1*w(j)
  enddo
  S_inv($n-1,i  ) = S_inv($n-1,i  )*lambda -zj*w($n-1)
  S_inv($n-1,i+1) = S_inv($n-1,i+1)*lambda -zj1*w($n-1)
  S_inv($n  ,i  ) = S_inv($n  ,i  )*lambda -zj*w($n  )
  S_inv($n  ,i+1) = S_inv($n  ,i+1)*lambda -zj1*w($n  )

end

SUBST [ n ]
6 ;;
10 ;;
14 ;;
18 ;;
22 ;;
26 ;;
30 ;;
34 ;;
38 ;;
42 ;;
46 ;;
50 ;;
54 ;;
58 ;;
62 ;;
66 ;;
70 ;;
74 ;;
78 ;;
82 ;;
86 ;;
90 ;;
94 ;;
98 ;;
102 ;;
106 ;;
110 ;;
114 ;;
118 ;;
122 ;;
126 ;;
130 ;;
134 ;;
138 ;;
142 ;;
146 ;;
150 ;;

END_TEMPLATE

BEGIN_TEMPLATE
! Version for mod(n,4) = 3
subroutine det_update$n(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS  ! Dimension of the vector
  real, intent(in)               :: m($n)  ! New vector
  integer, intent(in)            :: l      ! New position in S

  real,intent(inout)             :: S(LDS,$n)     ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,$n) ! Inverse Slater matrix
  double precision,intent(inout) :: d            ! Det(S)

  double precision               :: u($n), z($n), w($n), lambda, d_inv
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: z, w, u
  !DIR$ ASSUME_ALIGNED S : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED S_inv : $IRP_ALIGN
  !DIR$ ASSUME (mod(LDS,$IRP_ALIGN/8) == 0)
  !DIR$ ASSUME (LDS >= $n)
  integer                        :: i,j

  double precision :: zj, zj1, zj2, zj3

  do i=1,$n
    u(i) = m(i) - S(i,l)
  enddo

  zj = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do i=1,$n-3,4
    zj = zj + S_inv(i,l)*u(i) + S_inv(i+1,l)*u(i+1)  &
            + S_inv(i+2,l)*u(i+2) + S_inv(i+3,l)*u(i+3)
  enddo
  i=$n-2
  zj = zj + S_inv(i,l)*u(i) + S_inv(i+1,l)*u(i+1) + S_inv(i+2,l)*u(i+2)


  d_inv = 1.d0/d
  d = d+zj
  lambda = d*d_inv
  if (dabs(lambda) < 1.d-3) then
    d = 0.d0
    return
  endif

  !DIR$ VECTOR ALIGNED
  do j=1,$n-3,4
   zj  = 0.d0
   zj1 = 0.d0
   zj2 = 0.d0
   zj3 = 0.d0
   !DIR$ VECTOR ALIGNED
   do i=1,$n-3
    zj  = zj  + S_inv(i,j  )*u(i)
    zj1 = zj1 + S_inv(i,j+1)*u(i)
    zj2 = zj2 + S_inv(i,j+2)*u(i)
    zj3 = zj3 + S_inv(i,j+3)*u(i)
   enddo
   z(j  ) = zj     + S_inv($n-2,j  )*u($n-2)
   z(j  ) = z(j  ) + S_inv($n-1,j  )*u($n-1)
   z(j  ) = z(j  ) + S_inv($n,j  )*u($n)
   z(j+1) = zj1    + S_inv($n-2,j+1)*u($n-2)
   z(j+1) = z(j+1) + S_inv($n-1,j+1)*u($n-1)
   z(j+1) = z(j+1) + S_inv($n,j+1)*u($n)
   z(j+2) = zj2    + S_inv($n-2,j+2)*u($n-2)
   z(j+2) = z(j+2) + S_inv($n-1,j+2)*u($n-1)
   z(j+2) = z(j+2) + S_inv($n,j+2)*u($n)
   z(j+3) = zj3    + S_inv($n-2,j+3)*u($n-2)
   z(j+3) = z(j+3) + S_inv($n-1,j+3)*u($n-1)
   z(j+3) = z(j+3) + S_inv($n,j+3)*u($n)
  enddo

  j=$n-2
  zj  = 0.d0
  zj1 = 0.d0
  zj2 = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do i=1,$n-3
   zj  = zj  + S_inv(i,j  )*u(i)
   zj1 = zj1 + S_inv(i,j+1)*u(i)
   zj2 = zj2 + S_inv(i,j+2)*u(i)
  enddo
  z(j  ) = zj     + S_inv($n-2,j  )*u($n-2)
  z(j  ) = z(j  ) + S_inv($n-1,j  )*u($n-1)
  z(j  ) = z(j  ) + S_inv($n,j  )*u($n)
  z(j+1) = zj1    + S_inv($n-2,j+1)*u($n-2)
  z(j+1) = z(j+1) + S_inv($n-1,j+1)*u($n-1)
  z(j+1) = z(j+1) + S_inv($n,j+1)*u($n)
  z(j+2) = zj2    + S_inv($n-2,j+2)*u($n-2)
  z(j+2) = z(j+2) + S_inv($n-1,j+2)*u($n-1)
  z(j+2) = z(j+2) + S_inv($n,j+2)*u($n)

  !DIR$ NOPREFETCH
  do i=1,$n
    w(i) = S_inv(i,l)*d_inv
    S(i,l) = m(i)
  enddo

  do i=1,$n-3,4
   zj  = z(i)
   zj1 = z(i+1)
   zj2 = z(i+2)
   zj3 = z(i+3)
   !DIR$ VECTOR ALIGNED
   !DIR$ NOPREFETCH
   do j=1,$n-3
    S_inv(j,i  ) = S_inv(j,i  )*lambda - w(j)*zj
    S_inv(j,i+1) = S_inv(j,i+1)*lambda - w(j)*zj1
    S_inv(j,i+2) = S_inv(j,i+2)*lambda - w(j)*zj2
    S_inv(j,i+3) = S_inv(j,i+3)*lambda - w(j)*zj3
   enddo
   S_inv($n-2,i  ) = S_inv($n-2,i  )*lambda -zj *w($n-2)
   S_inv($n-1,i  ) = S_inv($n-1,i  )*lambda -zj *w($n-1)
   S_inv($n  ,i  ) = S_inv($n  ,i  )*lambda -zj *w($n  )
   S_inv($n-2,i+1) = S_inv($n-2,i+1)*lambda -zj1*w($n-2)
   S_inv($n-1,i+1) = S_inv($n-1,i+1)*lambda -zj1*w($n-1)
   S_inv($n  ,i+1) = S_inv($n  ,i+1)*lambda -zj1*w($n  )
   S_inv($n-2,i+2) = S_inv($n-2,i+2)*lambda -zj2*w($n-2)
   S_inv($n-1,i+2) = S_inv($n-1,i+2)*lambda -zj2*w($n-1)
   S_inv($n  ,i+2) = S_inv($n  ,i+2)*lambda -zj2*w($n  )
   S_inv($n-2,i+3) = S_inv($n-2,i+3)*lambda -zj3*w($n-2)
   S_inv($n-1,i+3) = S_inv($n-1,i+3)*lambda -zj3*w($n-1)
   S_inv($n  ,i+3) = S_inv($n  ,i+3)*lambda -zj3*w($n  )
  enddo

  i=$n-2
  zj  = z(i)
  zj1 = z(i+1)
  zj2 = z(i+2)
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do j=1,$n
   S_inv(j,i  ) = S_inv(j,i  )*lambda - w(j)*zj
   S_inv(j,i+1) = S_inv(j,i+1)*lambda - w(j)*zj1
   S_inv(j,i+2) = S_inv(j,i+2)*lambda - w(j)*zj2
  enddo


end

SUBST [ n ]
7 ;;
11 ;;
15 ;;
19 ;;
23 ;;
27 ;;
31 ;;
35 ;;
39 ;;
43 ;;
47 ;;
51 ;;
55 ;;
59 ;;
63 ;;
67 ;;
71 ;;
75 ;;
79 ;;
83 ;;
87 ;;
91 ;;
95 ;;
99 ;;
103 ;;
107 ;;
111 ;;
115 ;;
119 ;;
123 ;;
127 ;;
131 ;;
135 ;;
139 ;;
143 ;;
147 ;;

END_TEMPLATE



subroutine det_update_general(n,LDS,m,l,S,S_inv,d)
  implicit none

  integer, intent(in)            :: n,LDS      ! Dimension of the vector
  real, intent(in)               :: m(LDS)     ! New vector
  integer, intent(in)            :: l          ! New position in S

  real,intent(inout)             :: S(LDS,n)       ! Slater matrix
  double precision,intent(inout) :: S_inv(LDS,n)   ! Inverse Slater matrix
  double precision,intent(inout) :: d              ! Det(S)

  double precision               :: lambda, d_inv
  double precision               :: u(3840), z(3840), w(3840)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: z, w, u
  !DIR$ ASSUME_ALIGNED S : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED S_inv : $IRP_ALIGN
  !DIR$ ASSUME (LDS >= n)
  !DIR$ ASSUME (LDS <= 3840)
  !DIR$ ASSUME (MOD(LDS,$IRP_ALIGN/8) == 0)
  !DIR$ ASSUME (n>150)

  integer                        :: i,j,n4
  double precision :: zl

  !DIR$ NOPREFETCH
  do i=1,n
    u(i) = m(i) - S(i,l)
  enddo

  zl = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do i=1,n
   zl = zl + S_inv(i,l)*u(i)
  enddo

  d_inv = 1.d0/d
  d = d+zl
  lambda = d*d_inv
!
  if ( dabs(lambda) < 1.d-3 ) then
    d = 0.d0
    return
  endif

  double precision :: zj, zj1, zj2, zj3

  n4 = iand(n,not(3))
  !DIR$ VECTOR ALIGNED
  !DIR$ NOPREFETCH
  do j=1,n4,4
   zj  = 0.d0
   zj1 = 0.d0
   zj2 = 0.d0
   zj3 = 0.d0
   !DIR$ VECTOR ALIGNED
   !DIR$ NOPREFETCH
   do i=1,n
    zj  = zj  + S_inv(i,j  )*u(i)
    zj1 = zj1 + S_inv(i,j+1)*u(i)
    zj2 = zj2 + S_inv(i,j+2)*u(i)
    zj3 = zj3 + S_inv(i,j+3)*u(i)
   enddo
   z(j  ) = zj
   z(j+1) = zj1
   z(j+2) = zj2
   z(j+3) = zj3
  enddo

  do j=n4+1,n
   zj = 0.d0
   !DIR$ VECTOR ALIGNED
   !DIR$ NOPREFETCH
   do i=1,n
    zj = zj + S_inv(i,j)*u(i)
   enddo
   z(j  ) = zj
  enddo

  !DIR$ NOPREFETCH
  do i=1,n
    w(i) = S_inv(i,l)*d_inv
    S(i,l) = m(i)
  enddo

  !DIR$ NOPREFETCH
  do i=1,n
    w(i) = S_inv(i,l)*d_inv
    S(i,l) = m(i)
  enddo

  do i=1,n4,4
    zj  = z(i)
    zj1 = z(i+1)
    zj2 = z(i+2)
    zj3 = z(i+3)
    !DIR$ VECTOR ALIGNED
    !DIR$ NOPREFETCH
    do j=1,n
     S_inv(j,i  ) = S_inv(j,i  )*lambda -zj *w(j)
     S_inv(j,i+1) = S_inv(j,i+1)*lambda -zj1*w(j)
     S_inv(j,i+2) = S_inv(j,i+2)*lambda -zj2*w(j)
     S_inv(j,i+3) = S_inv(j,i+3)*lambda -zj3*w(j)
    enddo
  enddo

  do i=n4+1,n
    zj = z(i)
    !DIR$ VECTOR ALIGNED
    !DIR$ NOPREFETCH
    do j=1,n
      S_inv(j,i) = S_inv(j,i)*lambda -zj*w(j)
    enddo
  enddo

end



subroutine bitstring_to_list( string, list, n_elements, Nint)
  implicit none
  BEGIN_DOC
  ! Gives the inidices(+1) of the bits set to 1 in the bit string
  END_DOC
  integer, intent(in)            :: Nint
  integer*8, intent(in)          :: string(Nint)
  integer, intent(out)           :: list(Nint*64)
  integer, intent(out)           :: n_elements

  integer                        :: i, ishift
  integer*8                      :: l

  n_elements = 0
  ishift = 2
  do i=1,Nint
    l = string(i)
    do while (l /= 0_8)
      n_elements = n_elements+1
      list(n_elements) = ishift+popcnt(l-1_8) - popcnt(l)
      l = iand(l,l-1_8)
    enddo
    ishift = ishift + 64
  enddo
end

subroutine get_excitation_degree_spin(key1,key2,degree,Nint)
  implicit none
  BEGIN_DOC
  ! Returns the excitation degree between two determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(8), intent(in)  :: key1(Nint)
  integer(8), intent(in)  :: key2(Nint)
  integer, intent(out)           :: degree

  integer(8)              :: xorvec(32)
  integer                        :: l

  ASSERT (Nint > 0)
  ASSERT (Nint <= 32)

  select case (Nint)

    case (1)
      xorvec(1) = xor( key1(1), key2(1))
      degree = popcnt(xorvec(1))

    case (2)
      xorvec(1) = xor( key1(1), key2(1))
      xorvec(2) = xor( key1(2), key2(2))
      degree = popcnt(xorvec(1))+popcnt(xorvec(2))

    case (3)
      xorvec(1) = xor( key1(1), key2(1))
      xorvec(2) = xor( key1(2), key2(2))
      xorvec(3) = xor( key1(3), key2(3))
      degree = sum(popcnt(xorvec(1:3)))

    case (4)
      xorvec(1) = xor( key1(1), key2(1))
      xorvec(2) = xor( key1(2), key2(2))
      xorvec(3) = xor( key1(3), key2(3))
      xorvec(4) = xor( key1(4), key2(4))
      degree = sum(popcnt(xorvec(1:4)))

    case default
      do l=1,Nint
        xorvec(l) = xor( key1(l), key2(l))
      enddo
      degree = sum(popcnt(xorvec(1:Nint)))

  end select

  degree = shiftr(degree,1)

end

 BEGIN_PROVIDER [ integer, mo_list_alpha_curr, (elec_alpha_num) ]
&BEGIN_PROVIDER [ integer, mo_list_alpha_prev, (elec_alpha_num) ]
&BEGIN_PROVIDER [ integer, mo_exc_alpha_curr ]
 implicit none
 BEGIN_DOC
 ! List of MOs in the current alpha determinant
 END_DOC
 integer                        :: l
 if (det_i /= det_alpha_order(1)) then
   mo_list_alpha_prev = mo_list_alpha_curr
 else
   mo_list_alpha_prev = 0
 endif
 !DIR$ FORCEINLINE
 call bitstring_to_list ( psi_det_alpha(1,det_i), mo_list_alpha_curr, l, N_int )
 if (l /= elec_alpha_num) then
   stop 'error in number of alpha electrons'
 endif
 call get_excitation_degree_spin(psi_det_alpha(1,det_i),    &
                                 psi_det_alpha(1,det_i_prev), &
                                 mo_exc_alpha_curr, N_int)

END_PROVIDER


 BEGIN_PROVIDER [ integer, mo_list_beta_curr, (elec_beta_num) ]
&BEGIN_PROVIDER [ integer, mo_list_beta_prev, (elec_beta_num) ]
&BEGIN_PROVIDER [ integer, mo_exc_beta_curr ]
 implicit none
 BEGIN_DOC
 ! List of MOs in the current beta determinant
 END_DOC
 integer                        :: l
 if (elec_beta_num == 0) then
   return
 endif
 if (det_j /= det_beta_order(1)) then
   mo_list_beta_prev = mo_list_beta_curr
 else
   mo_list_beta_prev = 0
 endif

 !DIR$ FORCEINLINE
 call bitstring_to_list ( psi_det_beta(1,det_j), mo_list_beta_curr, l, N_int )
 if (l /= elec_beta_num) then
   stop 'error in number of beta electrons'
 endif
 call get_excitation_degree_spin(psi_det_beta(1,det_j),    &
                                 psi_det_beta(1,det_j_prev), &
                                 mo_exc_beta_curr, N_int)

END_PROVIDER

 BEGIN_PROVIDER [ double precision, det_alpha_value_curr_qmckl ]
&BEGIN_PROVIDER [ real, slater_matrix_alpha_qmckl, (elec_alpha_num_8,elec_alpha_num) ]
&BEGIN_PROVIDER [ double precision, slater_matrix_alpha_inv_det_qmckl, (elec_alpha_num_8,elec_alpha_num) ]

   implicit none

   BEGIN_DOC
   ! det_alpha_value_curr_qmckl : Value of the current alpha determinant
   !
   ! slater_matrix_alpha_qmckl : Slater matrix for the current alpha determinant.
   !  1st index runs over electrons and
   !  2nd index runs over MOs.
   !  Built with 1st determinant
   !
   ! slater_matrix_alpha_inv_det_qmckl: Inverse of the alpha Slater matrix x determinant
   END_DOC

   use qmckl

   double precision               :: ddet
   integer                        :: i,j,k,imo,l,nupds_selected,series
   integer                        :: to_do(elec_alpha_num), n_to_do_old, n_to_do
   double precision               :: tmp_inv(elec_alpha_num_8)
   real                           :: tmp_det(elec_alpha_num_8)
   integer, save                  :: ifirst, series_counter = 1

   double precision               :: cpu0, cpu1
   double precision, save         :: cpu_accumulator = 0.0d0
   integer (qmckl_exit_code)      :: rc
   integer(kind=8)                :: nupdates
   real(c_double)                 :: breakdown
   real(c_double)                 :: updates(elec_alpha_num_8, elec_alpha_num)
   !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: tmp_inv, tmp_det

   if (ifirst == 0) then
     ifirst = 1
     !DIR$ VECTOR ALIGNED
     slater_matrix_alpha_qmckl = 0.
     !DIR$ VECTOR ALIGNED
     slater_matrix_alpha_inv_det_qmckl = 0.d0
   endif ! (1)
   PROVIDE mo_value

   ! allocate(updates(elec_alpha_num_8, elec_alpha_num))
   updates = 0.0d0 !! Needed for zero-padding / correct local energy

   call cpu_time(cpu0)

   if (det_i /= det_alpha_order(1) ) then ! alpha determinant order has changed

     n_to_do = 0
     do k=1,elec_alpha_num ! run over all alpha electrons
       imo = mo_list_alpha_curr(k) !  pick MO-number from list of alpha-MOs
       if ( imo /= mo_list_alpha_prev(k) ) then ! check if MO-number of kth alpha electron has changed
         ! write(*,*) "(2a.1) imo DIFF FROM mp_list_alpha_prev(k) : imo = ", imo
         n_to_do += 1        ! update MO of kth alpha electron
         to_do(n_to_do) = k  ! remember which electrons to update
       endif
     enddo

     ! make swaps and keep 1 update
     if (n_to_do > 1 .and. mo_exc_alpha_curr == 1) then

       if (iand(n_to_do+1,1)==1) then !! Test if n_to_do is even
         ! write(*,*) "(2a.2.1) n_to_do is EVEN"
         det_alpha_value_curr_qmckl = -det_alpha_value_curr_qmckl
         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         slater_matrix_alpha_inv_det_qmckl = - slater_matrix_alpha_inv_det_qmckl
       endif

       if (mo_list_alpha_curr(to_do(1)) == mo_list_alpha_prev(to_do(1)+1)) then

         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         do l=1,elec_alpha_num_8
           tmp_det(l) = slater_matrix_alpha_qmckl(l,to_do(1))
           tmp_inv(l) = slater_matrix_alpha_inv_det_qmckl(l,to_do(1))
         enddo

         do k=to_do(1),to_do(n_to_do-1)
           !DIR$ VECTOR ALWAYS
           !DIR$ VECTOR ALIGNED
           do l=1,elec_alpha_num_8
             slater_matrix_alpha_qmckl(l,k) = slater_matrix_alpha_qmckl(l,k+1)
             slater_matrix_alpha_inv_det_qmckl(l,k) = slater_matrix_alpha_inv_det_qmckl(l,k+1)
           enddo
         enddo
         k = to_do(n_to_do)
         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         do l=1,elec_alpha_num_8
           slater_matrix_alpha_qmckl(l,k) = tmp_det(l)
           slater_matrix_alpha_inv_det_qmckl(l,k) = tmp_inv(l)
         enddo
         to_do(1) = to_do(n_to_do)

       else if (mo_list_alpha_curr(to_do(n_to_do)) == mo_list_alpha_prev(to_do(n_to_do)-1)) then
         ! write(*,*) "(2a.2.2b) mo_list_alpha_curr(to_do(n_to_do)) == mo_list_alpha_prev(to_do(n_to_do)-1)"
         k = to_do(n_to_do)
         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         do l=1,elec_alpha_num_8
           tmp_det(l) = slater_matrix_alpha_qmckl(l,k)
           tmp_inv(l) = slater_matrix_alpha_inv_det_qmckl(l,k)
         enddo
         do k=to_do(n_to_do),to_do(2),-1
           !DIR$ VECTOR ALWAYS
           !DIR$ VECTOR ALIGNED
           do l=1,elec_alpha_num_8
             slater_matrix_alpha_qmckl(l,k) = slater_matrix_alpha_qmckl(l,k-1)
             slater_matrix_alpha_inv_det_qmckl(l,k) = slater_matrix_alpha_inv_det_qmckl(l,k-1)
           enddo
         enddo
         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         do l=1,elec_alpha_num_8
           slater_matrix_alpha_qmckl(l,to_do(1)) = tmp_det(l)
           slater_matrix_alpha_inv_det_qmckl(l,to_do(1)) = tmp_inv(l)
         enddo

       endif
       n_to_do = 1
     endif

     ddet = 0.d0

     if (n_to_do < shiftl(elec_alpha_num,1)) then
       ! write(*,*) "(2a.3) n_to_do < 2 * elec_alpha_num"

       ddet = det_alpha_value_curr_qmckl ! set ddet to the current value
       slater_matrix_alpha_inv_det_qmckl = slater_matrix_alpha_inv_det_qmckl / ddet


       do j = 1, n_to_do   ! for all updates do ...
         k = to_do(j)            ! select the electron that needs an update
         imo = mo_list_alpha_curr(k) ! select the MO for electron k
         do i = 1, elec_alpha_num  ! run over all electrons
           updates(i, j) = mo_value(i, imo) - slater_matrix_alpha_qmckl(i, k)
         end do
       end do

       do j = 1, n_to_do   ! for all updates do ...
         k = to_do(j)            ! select the electron that needs an update
         imo = mo_list_alpha_curr(k) ! select the MO for electron k
         do i = 1, elec_alpha_num  ! run over all electrons
           slater_matrix_alpha_qmckl(i, k) = mo_value(i, imo) ! replacement upds applied to cols of S
         end do
       end do

       breakdown = 1d-3
! integer*8 :: context
! context = qmckl_context_create()
! rc = qmckl_sherman_morrison_smw32s(context,                 &
       rc = qmckl_sherman_morrison_smw32s(qmckl_ctx,                 &
           int(elec_alpha_num_8, kind=8),                            &
           int(elec_alpha_num, kind=8),                              &
           int(n_to_do, kind=8),                                     &
           updates,                                                  &
           int(to_do, kind=8),                                       &
           breakdown,                                                &
           slater_matrix_alpha_inv_det_qmckl,                              &
           ddet)
       rc = qmckl_check(qmckl_ctx, rc)
!  rc = qmckl_context_destroy(context)

       slater_matrix_alpha_inv_det_qmckl = slater_matrix_alpha_inv_det_qmckl * ddet
       det_alpha_value_curr_qmckl = ddet
     endif

   else
     ! write(*,*) "(2b) det_i SAME AS det_alpha_order(1). Set ddet to 0 : det_i, ddet (just before) = ", det_i, ", ", ddet
     ddet = 0.d0

   endif

   ! Avoid NaN
   if (ddet /= 0.d0) then
     ! write(*,*) "(3a) ddet DIFF FROM zero. continue : ddet = ", ddet
     continue
   else
     ! write(*,*) "(3b) ddet EQUAL TO zero. lapack time"
     do j=1,mo_closed_num
       !DIR$ VECTOR ALIGNED
       !DIR$ LOOP COUNT(100)
       do i=1,elec_alpha_num
         slater_matrix_alpha_qmckl(i,j) = mo_value(i,j)
         slater_matrix_alpha_inv_det_qmckl(j,i) = mo_value(i,j)
       enddo
     enddo
     do k=mo_closed_num+1,elec_alpha_num
       !DIR$ VECTOR ALIGNED
       !DIR$ LOOP COUNT(100)
       do i=1,elec_alpha_num
         slater_matrix_alpha_qmckl(i,k) = mo_value(i,mo_list_alpha_curr(k))
         slater_matrix_alpha_inv_det_qmckl(k,i) = mo_value(i,mo_list_alpha_curr(k))
       enddo
     enddo

     ! write(*,*) "CALLING LAPACK..."
     call invert(slater_matrix_alpha_inv_det_qmckl,elec_alpha_num_8,elec_alpha_num,ddet)

   endif

   call cpu_time(cpu1)
   cpu_accumulator += cpu1 - cpu0

   ASSERT (ddet /= 0.d0)

   det_alpha_value_curr_qmckl = ddet

END_PROVIDER


 BEGIN_PROVIDER [ double precision, det_alpha_value_curr ]
&BEGIN_PROVIDER [ real, slater_matrix_alpha, (elec_alpha_num_8,elec_alpha_num) ]
&BEGIN_PROVIDER [ double precision, slater_matrix_alpha_inv_det, (elec_alpha_num_8,elec_alpha_num) ]

  implicit none

  BEGIN_DOC
  ! det_alpha_value_curr : Value of the current alpha determinant
  !
  ! det_alpha_value_curr : Slater matrix for the current alpha determinant.
  !  1st index runs over electrons and
  !  2nd index runs over MOs.
  !  Built with 1st determinant
  !
  ! slater_matrix_alpha_inv_det: Inverse of the alpha Slater matrix x determinant
  END_DOC

  double precision               :: ddet
  integer                        :: i,j,k,imo,l
  integer                        :: to_do(elec_alpha_num), n_to_do_old, n_to_do
  double precision               :: tmp_inv(elec_alpha_num_8)
  real                           :: tmp_det(elec_alpha_num_8)
  integer, save                  :: ifirst
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: tmp_inv, tmp_det

  if (ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    slater_matrix_alpha = 0.
    !DIR$ VECTOR ALIGNED
    slater_matrix_alpha_inv_det = 0.d0
  endif

  PROVIDE mo_value
  if (det_i /= det_alpha_order(1) ) then

    n_to_do = 0
    do k=1,elec_alpha_num
      imo = mo_list_alpha_curr(k)
      if ( imo /= mo_list_alpha_prev(k) ) then
          n_to_do += 1
          to_do(n_to_do) = k
      endif
    enddo

    ! make swaps and keep 1 update
    if (n_to_do > 1 .and. mo_exc_alpha_curr == 1) then

      if (iand(n_to_do+1,1)==1) then
        det_alpha_value_curr = -det_alpha_value_curr
        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        slater_matrix_alpha_inv_det = - slater_matrix_alpha_inv_det
      endif

      if (mo_list_alpha_curr(to_do(1)) == mo_list_alpha_prev(to_do(1)+1)) then

        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        do l=1,elec_alpha_num_8
          tmp_det(l) = slater_matrix_alpha(l,to_do(1))
          tmp_inv(l) = slater_matrix_alpha_inv_det(l,to_do(1))
        enddo

        do k=to_do(1),to_do(n_to_do-1)
          !DIR$ VECTOR ALWAYS
          !DIR$ VECTOR ALIGNED
          do l=1,elec_alpha_num_8
            slater_matrix_alpha(l,k) = slater_matrix_alpha(l,k+1)
            slater_matrix_alpha_inv_det(l,k) = slater_matrix_alpha_inv_det(l,k+1)
          enddo
        enddo
        k = to_do(n_to_do)
        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        do l=1,elec_alpha_num_8
          slater_matrix_alpha(l,k) = tmp_det(l)
          slater_matrix_alpha_inv_det(l,k) = tmp_inv(l)
        enddo
        to_do(1) = to_do(n_to_do)

      else if (mo_list_alpha_curr(to_do(n_to_do)) == mo_list_alpha_prev(to_do(n_to_do)-1)) then
        k = to_do(n_to_do)
        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        do l=1,elec_alpha_num_8
          tmp_det(l) = slater_matrix_alpha(l,k)
          tmp_inv(l) = slater_matrix_alpha_inv_det(l,k)
        enddo
        do k=to_do(n_to_do),to_do(2),-1
          !DIR$ VECTOR ALWAYS
          !DIR$ VECTOR ALIGNED
          do l=1,elec_alpha_num_8
            slater_matrix_alpha(l,k) = slater_matrix_alpha(l,k-1)
            slater_matrix_alpha_inv_det(l,k) = slater_matrix_alpha_inv_det(l,k-1)
          enddo
        enddo
        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        do l=1,elec_alpha_num_8
          slater_matrix_alpha(l,to_do(1)) = tmp_det(l)
          slater_matrix_alpha_inv_det(l,to_do(1)) = tmp_inv(l)
        enddo

      endif
      n_to_do = 1
    endif

    ddet = 0.d0

    if (n_to_do < shiftl(elec_alpha_num,1)) then

      do while ( n_to_do > 0 )
        ddet = det_alpha_value_curr
        n_to_do_old = n_to_do
        n_to_do = 0
        do l=1,n_to_do_old
          k = to_do(l)
          imo = mo_list_alpha_curr(k)
          call det_update(elec_alpha_num, elec_alpha_num_8,            &
              mo_value(1,imo),                                         &
              k,                                                       &
              slater_matrix_alpha,                                     &
              slater_matrix_alpha_inv_det,                             &
              ddet)
          if (ddet /= 0.d0) then
            det_alpha_value_curr = ddet
          else
            n_to_do += 1
            to_do(n_to_do) = k
            ddet = det_alpha_value_curr
          endif
        enddo
        if (n_to_do == n_to_do_old) then
          ddet = 0.d0
          exit
        endif
      enddo

    endif

  else

    ddet = 0.d0

  endif

  ! Avoid NaN
  if (ddet /= 0.d0) then
    continue
  else
    do j=1,mo_closed_num
      !DIR$ VECTOR ALIGNED
      !DIR$ LOOP COUNT(100)
      do i=1,elec_alpha_num
        slater_matrix_alpha(i,j) = mo_value(i,j)
        slater_matrix_alpha_inv_det(j,i) = mo_value(i,j)
      enddo
    enddo
    do k=mo_closed_num+1,elec_alpha_num
      !DIR$ VECTOR ALIGNED
      !DIR$ LOOP COUNT(100)
      do i=1,elec_alpha_num
        slater_matrix_alpha(i,k) = mo_value(i,mo_list_alpha_curr(k))
        slater_matrix_alpha_inv_det(k,i) = mo_value(i,mo_list_alpha_curr(k))
      enddo
    enddo
    call invert(slater_matrix_alpha_inv_det,elec_alpha_num_8,elec_alpha_num,ddet)

  endif
  ASSERT (ddet /= 0.d0)

  det_alpha_value_curr = ddet
END_PROVIDER


 BEGIN_PROVIDER [ double precision, det_beta_value_curr_qmckl ]
&BEGIN_PROVIDER [ real, slater_matrix_beta_qmckl, (elec_beta_num_8,elec_beta_num) ]
&BEGIN_PROVIDER [ double precision, slater_matrix_beta_inv_det_qmckl, (elec_beta_num_8,elec_beta_num) ]
   BEGIN_DOC
   !  det_beta_value_curr_qmckl : Value of the current beta determinant
   !
   !  slater_matrix_beta_qmckl : Slater matrix for the current beta determinant.
   !  1st index runs over electrons and
   !  2nd index runs over MOs.
   !  Built with 1st determinant
   !
   !  slater_matrix_beta_inv_det_qmckl : Inverse of the beta Slater matrix x determinant
   END_DOC

   use qmckl

   double precision               :: ddet
   integer                        :: i,j,k,imo,l
   integer                        :: to_do(elec_alpha_num - mo_closed_num), n_to_do_old, n_to_do
   double precision               :: tmp_inv(elec_alpha_num_8)
   real                           :: tmp_det(elec_alpha_num_8)
   integer, save                  :: ifirst

   integer (qmckl_exit_code)      :: rc
   integer(c_int64_t)             :: nupdates
   real(c_double)                 :: breakdown
   real(c_double)                 :: updates(elec_beta_num_8, elec_beta_num)

   if (elec_beta_num == 0) then
     det_beta_value_curr_qmckl = 0.d0
     return
   endif

   if (ifirst == 0) then
     ifirst = 1
     slater_matrix_beta_qmckl = 0.
     slater_matrix_beta_inv_det_qmckl = 0.d0
   endif

   PROVIDE mo_value

   ! allocate(updates(elec_alpha_num_8, elec_alpha_num))
   updates = 0.0d0 !! Needed for zero-padding / correct local energy

   if (det_j /= det_beta_order(1)) then
     n_to_do = 0 !! Number of updates to apply
     do k=mo_closed_num+1,elec_beta_num
       imo = mo_list_beta_curr(k) !! current MO
       if ( imo /= mo_list_beta_prev(k) ) then
         n_to_do += 1 !! One more update to apply
         to_do(n_to_do) = k !! MO-list number to change
       endif
     enddo

     ! make swaps and keep 1 update
     if (n_to_do > 1 .and. mo_exc_beta_curr == 1) then

       if (iand(n_to_do+1,1)==1) then
         det_beta_value_curr_qmckl = -det_beta_value_curr_qmckl
         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         slater_matrix_beta_inv_det_qmckl = - slater_matrix_beta_inv_det_qmckl
       endif

       if (mo_list_beta_curr(to_do(1)) == mo_list_beta_prev(to_do(1)+1)) then

         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         do l=1,elec_beta_num_8
           tmp_det(l) = slater_matrix_beta_qmckl(l,to_do(1))
           tmp_inv(l) = slater_matrix_beta_inv_det_qmckl(l,to_do(1))
         enddo

         do k=to_do(1),to_do(n_to_do-1)
           !DIR$ VECTOR ALWAYS
           !DIR$ VECTOR ALIGNED
           do l=1,elec_beta_num_8
             slater_matrix_beta_qmckl(l,k) = slater_matrix_beta_qmckl(l,k+1)
             slater_matrix_beta_inv_det_qmckl(l,k) = slater_matrix_beta_inv_det_qmckl(l,k+1)
           enddo
         enddo
         k = to_do(n_to_do)
         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         do l=1,elec_beta_num_8
           slater_matrix_beta_qmckl(l,k) = tmp_det(l)
           slater_matrix_beta_inv_det_qmckl(l,k) = tmp_inv(l)
         enddo
         to_do(1) = to_do(n_to_do)

       else if (mo_list_beta_curr(to_do(n_to_do)) == mo_list_beta_prev(to_do(n_to_do)-1)) then
         k = to_do(n_to_do)
         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         do l=1,elec_beta_num_8
           tmp_det(l) = slater_matrix_beta_qmckl(l,k)
           tmp_inv(l) = slater_matrix_beta_inv_det_qmckl(l,k)
         enddo
         do k=to_do(n_to_do),to_do(2),-1
           !DIR$ VECTOR ALWAYS
           !DIR$ VECTOR ALIGNED
           do l=1,elec_beta_num_8
             slater_matrix_beta_qmckl(l,k) = slater_matrix_beta_qmckl(l,k-1)
             slater_matrix_beta_inv_det_qmckl(l,k) = slater_matrix_beta_inv_det_qmckl(l,k-1)
           enddo
         enddo
         !DIR$ VECTOR ALWAYS
         !DIR$ VECTOR ALIGNED
         do l=1,elec_beta_num_8
           slater_matrix_beta_qmckl(l,to_do(1)) = tmp_det(l)
           slater_matrix_beta_inv_det_qmckl(l,to_do(1)) = tmp_inv(l)
         enddo

       endif
       n_to_do = 1
     endif

     ddet = 0.d0

     if (n_to_do < shiftl(elec_beta_num,1)) then !! Why compare to double the number of electrons?
       ddet = det_beta_value_curr_qmckl ! set ddet to the current value
       slater_matrix_beta_inv_det_qmckl = slater_matrix_beta_inv_det_qmckl / ddet

       do j = 1, n_to_do   ! do for all updates
         k = to_do(j)            ! find back the electron we need to update
         imo = mo_list_beta_curr(k) ! find column/row/MO for electron k
         do i = 1, elec_beta_num  ! run over all electrons
           updates(i, j) = mo_value(elec_alpha_num + i, imo) - slater_matrix_beta_qmckl(i, k)
           slater_matrix_beta_qmckl(i, k) = mo_value(elec_alpha_num + i, imo)
         end do
       end do

       breakdown = 1d-3
! integer*8 :: context
! context = qmckl_context_create()
! rc = qmckl_sherman_morrison_smw32s(context,                 &
       rc = qmckl_sherman_morrison_smw32s(qmckl_ctx,                 &
           int(elec_beta_num_8, kind=8),                             &
           int(elec_beta_num, kind=8),                               &
           int(n_to_do, kind=8),                                     &
           updates,                                                  &
           int(to_do, kind=8),                                       &
           breakdown,                                                &
           slater_matrix_beta_inv_det_qmckl,                               &
           ddet)
       rc = qmckl_check(qmckl_ctx, rc)
!  rc = qmckl_context_destroy(context)

       slater_matrix_beta_inv_det_qmckl = slater_matrix_beta_inv_det_qmckl * ddet
       det_beta_value_curr_qmckl = ddet
     endif
   else
     ddet = 0.d0
   endif

   ! Avoid NaN
   if (ddet /= 0.d0) then
     continue
   else
     do j=1,mo_closed_num
       !DIR$ VECTOR UNALIGNED
       !DIR$ LOOP COUNT (100)
       do i=1,elec_beta_num
         slater_matrix_beta_qmckl(i,j) = mo_value(i+elec_alpha_num,j)
         slater_matrix_beta_inv_det_qmckl(j,i) = mo_value(i+elec_alpha_num,j)
       enddo
     enddo
     do k=mo_closed_num+1,elec_beta_num
       !DIR$ VECTOR UNALIGNED
       !DIR$ LOOP COUNT (100)
       do i=1,elec_beta_num
         slater_matrix_beta_qmckl(i,k) = mo_value(i+elec_alpha_num,mo_list_beta_curr(k))
         slater_matrix_beta_inv_det_qmckl(k,i) = mo_value(i+elec_alpha_num,mo_list_beta_curr(k))
       enddo
     enddo
     call invert(slater_matrix_beta_inv_det_qmckl,elec_beta_num_8,elec_beta_num,ddet)
   endif
   ASSERT (ddet /= 0.d0)

   det_beta_value_curr_qmckl = ddet

END_PROVIDER


 BEGIN_PROVIDER [ double precision, det_beta_value_curr ]
&BEGIN_PROVIDER [ real, slater_matrix_beta, (elec_beta_num_8,elec_beta_num) ]
&BEGIN_PROVIDER [ double precision, slater_matrix_beta_inv_det, (elec_beta_num_8,elec_beta_num) ]
  BEGIN_DOC
  !  det_beta_value_curr : Value of the current beta determinant
  !
  !  slater_matrix_beta : Slater matrix for the current beta determinant.
  !  1st index runs over electrons and
  !  2nd index runs over MOs.
  !  Built with 1st determinant
  !
  !  slater_matrix_beta_inv_det : Inverse of the beta Slater matrix x determinant
  END_DOC

  double precision               :: ddet
  integer                        :: i,j,k,imo,l
  integer                        :: to_do(elec_alpha_num-mo_closed_num), n_to_do_old, n_to_do
  double precision               :: tmp_inv(elec_alpha_num_8)
  real                           :: tmp_det(elec_alpha_num_8)
  integer, save                  :: ifirst

  if (elec_beta_num == 0) then
    det_beta_value_curr = 0.d0
    return
  endif

  if (ifirst == 0) then
    ifirst = 1
    slater_matrix_beta = 0.
    slater_matrix_beta_inv_det = 0.d0
  endif
  PROVIDE mo_value

  if (det_j /= det_beta_order(1)) then

    n_to_do = 0
    do k=mo_closed_num+1,elec_beta_num
      imo = mo_list_beta_curr(k)
      if ( imo /= mo_list_beta_prev(k) ) then
          n_to_do += 1
          to_do(n_to_do) = k
      endif
    enddo

    ! make swaps and keep 1 update
    if (n_to_do > 1 .and. mo_exc_beta_curr == 1) then

      if (iand(n_to_do+1,1)==1) then
        det_beta_value_curr = -det_beta_value_curr
        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        slater_matrix_beta_inv_det = - slater_matrix_beta_inv_det
      endif

      if (mo_list_beta_curr(to_do(1)) == mo_list_beta_prev(to_do(1)+1)) then

        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        do l=1,elec_beta_num_8
          tmp_det(l) = slater_matrix_beta(l,to_do(1))
          tmp_inv(l) = slater_matrix_beta_inv_det(l,to_do(1))
        enddo

        do k=to_do(1),to_do(n_to_do-1)
          !DIR$ VECTOR ALWAYS
          !DIR$ VECTOR ALIGNED
          do l=1,elec_beta_num_8
            slater_matrix_beta(l,k) = slater_matrix_beta(l,k+1)
            slater_matrix_beta_inv_det(l,k) = slater_matrix_beta_inv_det(l,k+1)
          enddo
        enddo
        k = to_do(n_to_do)
        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        do l=1,elec_beta_num_8
          slater_matrix_beta(l,k) = tmp_det(l)
          slater_matrix_beta_inv_det(l,k) = tmp_inv(l)
        enddo
        to_do(1) = to_do(n_to_do)

      else if (mo_list_beta_curr(to_do(n_to_do)) == mo_list_beta_prev(to_do(n_to_do)-1)) then
        k = to_do(n_to_do)
        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        do l=1,elec_beta_num_8
          tmp_det(l) = slater_matrix_beta(l,k)
          tmp_inv(l) = slater_matrix_beta_inv_det(l,k)
        enddo
        do k=to_do(n_to_do),to_do(2),-1
          !DIR$ VECTOR ALWAYS
          !DIR$ VECTOR ALIGNED
          do l=1,elec_beta_num_8
            slater_matrix_beta(l,k) = slater_matrix_beta(l,k-1)
            slater_matrix_beta_inv_det(l,k) = slater_matrix_beta_inv_det(l,k-1)
          enddo
        enddo
        !DIR$ VECTOR ALWAYS
        !DIR$ VECTOR ALIGNED
        do l=1,elec_beta_num_8
          slater_matrix_beta(l,to_do(1)) = tmp_det(l)
          slater_matrix_beta_inv_det(l,to_do(1)) = tmp_inv(l)
        enddo

      endif
      n_to_do = 1
    endif

    ddet = 0.d0

    if (n_to_do < shiftl(elec_beta_num,1)) then

      do while ( n_to_do > 0 )
        ddet = det_beta_value_curr
        n_to_do_old = n_to_do
        n_to_do = 0
        do l=1,n_to_do_old
          k = to_do(l)
          imo = mo_list_beta_curr(k)
          call det_update(elec_beta_num, elec_beta_num_8,                &
              mo_value(elec_alpha_num+1,imo),                            &
              k,                                                         &
              slater_matrix_beta,                                        &
              slater_matrix_beta_inv_det,                                &
              ddet)
          if (ddet /= 0.d0) then
            det_beta_value_curr = ddet
          else
            n_to_do += 1
            to_do(n_to_do) = k
            ddet = det_beta_value_curr
          endif
        enddo
        if (n_to_do == n_to_do_old) then
          ddet = 0.d0
          exit
        endif
      enddo

    endif

  else

    ddet = 0.d0

  endif

  ! Avoid NaN
  if (ddet /= 0.d0) then
    continue
  else
    do j=1,mo_closed_num
      !DIR$ VECTOR UNALIGNED
      !DIR$ LOOP COUNT (100)
      do i=1,elec_beta_num
        slater_matrix_beta(i,j) = mo_value(i+elec_alpha_num,j)
        slater_matrix_beta_inv_det(j,i) = mo_value(i+elec_alpha_num,j)
      enddo
    enddo
    do k=mo_closed_num+1,elec_beta_num
      !DIR$ VECTOR UNALIGNED
      !DIR$ LOOP COUNT (100)
      do i=1,elec_beta_num
        slater_matrix_beta(i,k) = mo_value(i+elec_alpha_num,mo_list_beta_curr(k))
        slater_matrix_beta_inv_det(k,i) = mo_value(i+elec_alpha_num,mo_list_beta_curr(k))
      enddo
    enddo
    call invert(slater_matrix_beta_inv_det,elec_beta_num_8,elec_beta_num,ddet)
  endif
  ASSERT (ddet /= 0.d0)

  det_beta_value_curr = ddet

END_PROVIDER



 BEGIN_PROVIDER [ integer, det_alpha_num_pseudo ]
&BEGIN_PROVIDER [ integer, det_beta_num_pseudo ]
 implicit none
 BEGIN_DOC
 ! Dimensioning of large arrays made smaller without pseudo
 END_DOC
 if (do_pseudo) then
    det_alpha_num_pseudo = det_alpha_num
    det_beta_num_pseudo = det_beta_num
 else
    det_alpha_num_pseudo = 1
    det_beta_num_pseudo = 1
 endif
END_PROVIDER


 BEGIN_PROVIDER [ double precision , det_alpha_value,  (det_alpha_num_8) ]
&BEGIN_PROVIDER [ double precision , det_alpha_grad_lapl, (4,elec_alpha_num,det_alpha_num) ]
&BEGIN_PROVIDER [ double precision , det_alpha_pseudo, (elec_alpha_num_8,det_alpha_num_pseudo) ]

  implicit none

  BEGIN_DOC
  ! Values of the alpha determinants
  ! Gradients of the alpha determinants
  ! Laplacians of the alpha determinants
  END_DOC

  integer                        :: j,i,k
  integer, save                  :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    det_alpha_value  = 0.d0
    det_alpha_grad_lapl = 0.d0
    det_alpha_pseudo = 0.d0
  endif


  do j=1,det_alpha_num

    det_i_prev = det_i
    det_i = det_alpha_order(j)
    if (j > 1) then
      TOUCH det_i
    endif

    det_alpha_value(det_i) = det_alpha_value_curr
    det_alpha_grad_lapl(1:4,1:elec_alpha_num,det_i) = det_alpha_grad_lapl_curr(1:4,1:elec_alpha_num)
    if (do_pseudo) then
      det_alpha_pseudo(1:elec_alpha_num,det_i) = det_alpha_pseudo_curr(1:elec_alpha_num)
    endif

  enddo

  det_i = det_alpha_order(1)
  det_i_prev = det_alpha_order(1)
  SOFT_TOUCH det_i det_i_prev

END_PROVIDER

 BEGIN_PROVIDER [ double precision, det_beta_value,  (det_beta_num_8) ]
&BEGIN_PROVIDER [ double precision, det_beta_grad_lapl, (4,elec_beta_num,det_beta_num) ]
&BEGIN_PROVIDER [ double precision, det_beta_pseudo, (elec_beta_num_8,det_beta_num_pseudo) ]


  implicit none

  BEGIN_DOC
  ! Values of the beta determinants
  ! Gradients of the beta determinants
  ! Laplacians of the beta determinants
  END_DOC

  integer                        :: j,i,k
  integer, save                  :: ifirst = 0
  if (elec_beta_num == 0) then
    det_beta_value = 1.d0
    return
  endif

  if (ifirst == 0) then
    ifirst = 1
    det_beta_value  = 0.d0
    det_beta_grad_lapl = 0.d0
    det_beta_pseudo = 0.d0
  endif

  do j=1,det_beta_num

    det_j_prev = det_j
    det_j = det_beta_order(j)
    if (j > 1) then
      TOUCH det_j
    endif

    det_beta_value(det_j)  = det_beta_value_curr
    det_beta_grad_lapl(1:4,1:elec_beta_num,det_j) = &
       det_beta_grad_lapl_curr(1:4,elec_alpha_num+1:elec_num)
    if (do_pseudo) then
      det_beta_pseudo(1:elec_beta_num,det_j) = &
        det_beta_pseudo_curr(elec_alpha_num+1:elec_num)
    endif

  enddo

  det_j = det_beta_order(1)
  det_j_prev = det_beta_order(1)
  SOFT_TOUCH det_j det_j_prev

END_PROVIDER

 BEGIN_PROVIDER [ double precision, det_alpha_lapl_sum, (det_alpha_num_8) ]
&BEGIN_PROVIDER [ double precision, det_beta_lapl_sum,  (det_beta_num_8) ]
  implicit none
  BEGIN_DOC
  ! Sum of Laplacian_i per spin-determinant
  END_DOC
  integer :: i, k

  do k=1,det_alpha_num
    det_alpha_lapl_sum(k) = sum(det_alpha_grad_lapl(4,1:elec_alpha_num,k))
  enddo
  do k=1,det_beta_num
    det_beta_lapl_sum(k) = sum(det_beta_grad_lapl(4,1:elec_beta_num,k))
  enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision, psidet_value ]
&BEGIN_PROVIDER [ double precision, psidet_inv ]
&BEGIN_PROVIDER [ double precision, psidet_grad_lapl, (4,elec_num_8) ]
&BEGIN_PROVIDER [ double precision, pseudo_non_local, (elec_num) ]
&BEGIN_PROVIDER [ double precision, CDb, (det_alpha_num_8) ]
&BEGIN_PROVIDER [ double precision, DaC, (det_beta_num_8) ]

  implicit none
  BEGIN_DOC
  ! Value of the determinantal part of the wave function

  ! Gradient of the determinantal part of the wave function

  ! Laplacian of determinantal part of the wave function

  ! Non-local component of the pseudopotentials

  ! Regularized 1/psi = 1/(psi + 1/psi)

  ! C x D_beta

  ! D_alpha^t x C

  ! D_alpha^t x (C x D_beta)
  END_DOC

  integer, save                  :: ifirst=0
  if (ifirst == 0) then
    ifirst = 1
    psidet_grad_lapl = 0.d0
    pseudo_non_local = 0.d0
  endif

  integer :: i,j,k, l
  integer :: i1,i2,i3,i4,det_num4
  integer :: j1,j2,j3,j4
  double precision :: f

  DaC = 0.d0
  CDb = 0.d0

  if (det_num < shiftr(det_alpha_num*det_beta_num,2)) then

    det_num4 = iand(det_num,not(3))
    !DIR$ VECTOR ALIGNED
    do k=1,det_num4,4
      i1 = det_coef_matrix_rows(k  )
      i2 = det_coef_matrix_rows(k+1)
      i3 = det_coef_matrix_rows(k+2)
      i4 = det_coef_matrix_rows(k+3)
      j1 = det_coef_matrix_columns(k  )
      j2 = det_coef_matrix_columns(k+1)
      j3 = det_coef_matrix_columns(k+2)
      j4 = det_coef_matrix_columns(k+3)
      if ( (j1 == j2).and.(j1 == j3).and.(j1 == j4) ) then
        f = det_beta_value (j1)
        CDb(i1) = CDb(i1) + det_coef_matrix_values(k  )*f
        CDb(i2) = CDb(i2) + det_coef_matrix_values(k+1)*f
        CDb(i3) = CDb(i3) + det_coef_matrix_values(k+2)*f
        CDb(i4) = CDb(i4) + det_coef_matrix_values(k+3)*f

        if ( ((i2-i1) == 1).and.((i3-i1) == 2).and.((i4-i1) == 3) ) then
          DaC(j1) = DaC(j1) + det_coef_matrix_values(k)*det_alpha_value(i1) &
          + det_coef_matrix_values(k+1)*det_alpha_value(i1+1) &
          + det_coef_matrix_values(k+2)*det_alpha_value(i1+2) &
          + det_coef_matrix_values(k+3)*det_alpha_value(i1+3)
        else
          DaC(j1) = DaC(j1) + det_coef_matrix_values(k)*det_alpha_value(i1) &
          + det_coef_matrix_values(k+1)*det_alpha_value(i2) &
          + det_coef_matrix_values(k+2)*det_alpha_value(i3) &
          + det_coef_matrix_values(k+3)*det_alpha_value(i4)
        endif
      else
        DaC(j1) = DaC(j1) + det_coef_matrix_values(k  )*det_alpha_value(i1)
        DaC(j2) = DaC(j2) + det_coef_matrix_values(k+1)*det_alpha_value(i2)
        DaC(j3) = DaC(j3) + det_coef_matrix_values(k+2)*det_alpha_value(i3)
        DaC(j4) = DaC(j4) + det_coef_matrix_values(k+3)*det_alpha_value(i4)
        CDb(i1) = CDb(i1) + det_coef_matrix_values(k  )*det_beta_value (j1)
        CDb(i2) = CDb(i2) + det_coef_matrix_values(k+1)*det_beta_value (j2)
        CDb(i3) = CDb(i3) + det_coef_matrix_values(k+2)*det_beta_value (j3)
        CDb(i4) = CDb(i4) + det_coef_matrix_values(k+3)*det_beta_value (j4)
      endif
    enddo

    do k=det_num4+1,det_num
      i = det_coef_matrix_rows(k)
      j = det_coef_matrix_columns(k)
      DaC(j) = DaC(j) + det_coef_matrix_values(k)*det_alpha_value(i)
      CDb(i) = CDb(i) + det_coef_matrix_values(k)*det_beta_value (j)
    enddo

  else

    if (det_num == 1) then

      DaC(1) = det_alpha_value_curr
      CDb(1) = det_beta_value_curr

    else if (use_svd) then

      DaC = 0.d0
      CDb = 0.d0
      double precision :: DaU(4), VtDb(4)
      integer :: n_svd4


      n_svd4 = iand(n_svd_coefs,not(3))

      do i=1,n_svd4,4

        DaU = 0.d0
        do j=1,det_alpha_num
          DaU(1) = DaU(1) + psi_svd_alpha(j,i+0) * det_alpha_value(j)
          DaU(2) = DaU(2) + psi_svd_alpha(j,i+1) * det_alpha_value(j)
          DaU(3) = DaU(3) + psi_svd_alpha(j,i+2) * det_alpha_value(j)
          DaU(4) = DaU(4) + psi_svd_alpha(j,i+3) * det_alpha_value(j)
        end do
        DaU(1) = DaU(1)*psi_svd_coefs(i+0)
        DaU(2) = DaU(2)*psi_svd_coefs(i+1)
        DaU(3) = DaU(3)*psi_svd_coefs(i+2)
        DaU(4) = DaU(4)*psi_svd_coefs(i+3)


        VtDb = 0.d0
        do j=1,det_beta_num
          VtDb(1) = VtDb(1) + psi_svd_beta(j,i+0) * det_beta_value(j)
          VtDb(2) = VtDb(2) + psi_svd_beta(j,i+1) * det_beta_value(j)
          VtDb(3) = VtDb(3) + psi_svd_beta(j,i+2) * det_beta_value(j)
          VtDb(4) = VtDb(4) + psi_svd_beta(j,i+3) * det_beta_value(j)

          DaC(j) = DaC(j) + DaU(1) * psi_svd_beta(j,i+0) + &
                            DaU(2) * psi_svd_beta(j,i+1) + &
                            DaU(3) * psi_svd_beta(j,i+2) + &
                            DaU(4) * psi_svd_beta(j,i+3)
        end do
        VtDb(1) = VtDb(1)*psi_svd_coefs(i+0)
        VtDb(2) = VtDb(2)*psi_svd_coefs(i+1)
        VtDb(3) = VtDb(3)*psi_svd_coefs(i+2)
        VtDb(4) = VtDb(4)*psi_svd_coefs(i+3)

        do j=1,det_alpha_num
          CDb(j) = CDb(j) + VtDb(1) * psi_svd_alpha(j,i+0) + &
                            VtDb(2) * psi_svd_alpha(j,i+1) + &
                            VtDb(3) * psi_svd_alpha(j,i+2) + &
                            VtDb(4) * psi_svd_alpha(j,i+3)
        end do

      end do


      do i=n_svd4+1,n_svd_coefs

        DaU(1) = 0.d0
        do j=1,det_alpha_num
          DaU(1) = DaU(1) + psi_svd_alpha(j,i) * det_alpha_value(j)
        end do
        DaU(1) = DaU(1)*psi_svd_coefs(i)

        VtDb(1) = 0.d0
        do j=1,det_beta_num
          VtDb(1) = VtDb(1) + psi_svd_beta(j,i) * det_beta_value(j)
          DaC(j) = DaC(j) + DaU(1) * psi_svd_beta(j,i)
        end do
        VtDb(1) = VtDb(1)*psi_svd_coefs(i)

        do j=1,det_alpha_num
          CDb(j) = CDb(j) + VtDb(1) * psi_svd_alpha(j,i)
        end do

      end do

    else

      call dgemv('T',det_alpha_num,det_beta_num,1.d0,det_coef_matrix_dense, &
        size(det_coef_matrix_dense,1), det_alpha_value, 1, 0.d0, DaC, 1)

      call dgemv('N',det_alpha_num,det_beta_num,1.d0,det_coef_matrix_dense, &
        size(det_coef_matrix_dense,1), det_beta_value, 1, 0.d0, CDb, 1)

    endif

  endif

  ! Value
  ! -----

  psidet_value = 0.d0
  do j=1,det_beta_num
    psidet_value = psidet_value + det_beta_value(j) * DaC(j)
  enddo


  if (psidet_value == 0.d0) then
    call abrt(irp_here,'Determinantal component of the wave function is zero.')
  endif
  psidet_inv = 1.d0/psidet_value


  ! Gradients
  ! ---------
  if(det_num .eq. 1) then
    do i = 1, elec_alpha_num
      psidet_grad_lapl(1,i) = det_alpha_grad_lapl_curr(1,i) * det_beta_value_curr
      psidet_grad_lapl(2,i) = det_alpha_grad_lapl_curr(2,i) * det_beta_value_curr
      psidet_grad_lapl(3,i) = det_alpha_grad_lapl_curr(3,i) * det_beta_value_curr
      psidet_grad_lapl(4,i) = det_alpha_grad_lapl_curr(4,i) * det_beta_value_curr
    enddo
    do i = elec_alpha_num+1, elec_num
      psidet_grad_lapl(1,i) = det_beta_grad_lapl_curr(1,i) * det_alpha_value_curr
      psidet_grad_lapl(2,i) = det_beta_grad_lapl_curr(2,i) * det_alpha_value_curr
      psidet_grad_lapl(3,i) = det_beta_grad_lapl_curr(3,i) * det_alpha_value_curr
      psidet_grad_lapl(4,i) = det_beta_grad_lapl_curr(4,i) * det_alpha_value_curr
    enddo
  else
    ! psidet_grad_lapl(4,elec_alpha_num) =
    ! det_alpha_grad_lapl(4,elec_alpha_num,det_alpha_num) @ CDb(det_alpha_num)
    call dgemv('N',elec_alpha_num*4,det_alpha_num,1.d0,                &
        det_alpha_grad_lapl,                                           &
        size(det_alpha_grad_lapl,1)*size(det_alpha_grad_lapl,2),       &
        CDb, 1, 0.d0, psidet_grad_lapl, 1)
    if (elec_beta_num /= 0) then
      call dgemv('N',elec_beta_num*4,det_beta_num,1.d0,                  &
          det_beta_grad_lapl,                                            &
          size(det_beta_grad_lapl,1)*size(det_beta_grad_lapl,2),         &
          DaC, 1, 0.d0, psidet_grad_lapl(1,elec_alpha_num+1), 1)
    endif
  endif


  if (do_pseudo) then
    call dgemv('N',elec_alpha_num,det_alpha_num,psidet_inv,          &
        det_alpha_pseudo, size(det_alpha_pseudo,1),                  &
        CDb, 1, 0.d0, pseudo_non_local, 1)
    if (elec_beta_num /= 0) then
      call dgemv('N',elec_beta_num,det_beta_num,psidet_inv,            &
          det_beta_pseudo, size(det_beta_pseudo,1),                    &
          DaC, 1, 0.d0, pseudo_non_local(elec_alpha_num+1), 1)
    endif
  endif

END_PROVIDER

BEGIN_PROVIDER  [ double precision, det_alpha_pseudo_curr, (elec_alpha_num) ]
  implicit none
  BEGIN_DOC
! Pseudopotential non-local contribution
  END_DOC
  integer                        :: i,j,l,m,k,n
  integer                        :: imo,kk
  double precision               :: c
  integer, save                  :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    det_alpha_pseudo_curr = 0.d0
  endif
  if (do_pseudo) then
    do i=1,elec_alpha_num
      det_alpha_pseudo_curr(i) = 0.d0
      do n=1,elec_alpha_num
        imo = mo_list_alpha_curr(n)
        c = slater_matrix_alpha_inv_det(i,n)
        det_alpha_pseudo_curr(i) =                                &
              det_alpha_pseudo_curr(i) + c*pseudo_mo_term(imo,i)
      enddo
    enddo
  endif
END_PROVIDER

BEGIN_PROVIDER  [ double precision, det_beta_pseudo_curr, (elec_alpha_num+1:elec_num) ]
  implicit none
  BEGIN_DOC
! Pseudopotential non-local contribution
  END_DOC
  integer                        :: i,j,l,m,k,n
  integer                        :: imo,kk
  double precision               :: c
  integer, save                  :: ifirst = 0
  if (elec_beta_num == 0) then
    return
  endif
  if (ifirst == 0) then
    ifirst = 1
    det_beta_pseudo_curr = 0.d0
  endif
  if (do_pseudo) then
    do i=elec_alpha_num+1,elec_num
      det_beta_pseudo_curr(i) = 0.d0
      do n=1,elec_beta_num
        imo = mo_list_beta_curr(n)
        c = slater_matrix_beta_inv_det(i-elec_alpha_num,n)
        det_beta_pseudo_curr(i) =                                    &
            det_beta_pseudo_curr(i) + c*pseudo_mo_term(imo,i)
      enddo
    enddo
  endif
END_PROVIDER

BEGIN_PROVIDER  [ double precision, det_alpha_grad_lapl_curr, (4,elec_alpha_num) ]
  implicit none
  BEGIN_DOC
  ! Gradient of the current alpha determinant
  END_DOC

  integer                        :: i, j, k
  !DIR$ VECTOR ALIGNED
  do i=1,elec_alpha_num
    det_alpha_grad_lapl_curr(1,i) = 0.d0
    det_alpha_grad_lapl_curr(2,i) = 0.d0
    det_alpha_grad_lapl_curr(3,i) = 0.d0
    det_alpha_grad_lapl_curr(4,i) = 0.d0
  enddo

  integer :: imo, imo2

! -------
! The following code does the same as this:
!
!    do j=1,elec_alpha_num
!      imo  = mo_list_alpha_curr(j)
!      do i=1,elec_alpha_num
!        do k=1,4
!          det_alpha_grad_lapl_curr(k,i) = &
!             det_alpha_grad_lapl_curr(k,i) + mo_grad_lapl_alpha(k,i,imo)*slater_matrix_alpha_inv_det(i,j)
!        enddo
!      enddo
!    enddo
!
! -------

  if (iand(elec_alpha_num,1) == 0) then

    do j=1,elec_alpha_num,2
      imo  = mo_list_alpha_curr(j  )
      imo2 = mo_list_alpha_curr(j+1)
      do i=1,elec_alpha_num,2
        !DIR$ VECTOR ALIGNED
        do k=1,4
          det_alpha_grad_lapl_curr(k,i  ) = det_alpha_grad_lapl_curr(k,i  )       &
            + mo_grad_lapl_alpha(k,i  ,imo )*slater_matrix_alpha_inv_det(i  ,j  ) &
            + mo_grad_lapl_alpha(k,i  ,imo2)*slater_matrix_alpha_inv_det(i  ,j+1)
          det_alpha_grad_lapl_curr(k,i+1) = det_alpha_grad_lapl_curr(k,i+1)       &
            + mo_grad_lapl_alpha(k,i+1,imo )*slater_matrix_alpha_inv_det(i+1,j  ) &
            + mo_grad_lapl_alpha(k,i+1,imo2)*slater_matrix_alpha_inv_det(i+1,j+1)
        enddo
      enddo
    enddo

  else

    do j=1,elec_alpha_num-1,2
      imo  = mo_list_alpha_curr(j  )
      imo2 = mo_list_alpha_curr(j+1)
      do i=1,elec_alpha_num-1,2
        !DIR$ VECTOR ALIGNED
        do k=1,4
          det_alpha_grad_lapl_curr(k,i  ) = det_alpha_grad_lapl_curr(k,i  )       &
            + mo_grad_lapl_alpha(k,i  ,imo )*slater_matrix_alpha_inv_det(i  ,j  ) &
            + mo_grad_lapl_alpha(k,i  ,imo2)*slater_matrix_alpha_inv_det(i  ,j+1)
          det_alpha_grad_lapl_curr(k,i+1) = det_alpha_grad_lapl_curr(k,i+1)       &
            + mo_grad_lapl_alpha(k,i+1,imo )*slater_matrix_alpha_inv_det(i+1,j  ) &
            + mo_grad_lapl_alpha(k,i+1,imo2)*slater_matrix_alpha_inv_det(i+1,j+1)
        enddo
      enddo
      i=elec_alpha_num
        !DIR$ VECTOR ALIGNED
      do k=1,4
        det_alpha_grad_lapl_curr(k,i) = det_alpha_grad_lapl_curr(k,i)       &
          + mo_grad_lapl_alpha(k,i,imo )*slater_matrix_alpha_inv_det(i,j  ) &
          + mo_grad_lapl_alpha(k,i,imo2)*slater_matrix_alpha_inv_det(i,j+1)
      enddo
    enddo

    j=elec_alpha_num
    imo  = mo_list_alpha_curr(j)
    do i=1,elec_alpha_num
        !DIR$ VECTOR ALIGNED
      do k=1,4
        det_alpha_grad_lapl_curr(k,i  ) = det_alpha_grad_lapl_curr(k,i  ) + &
           mo_grad_lapl_alpha(k,i  ,imo)*slater_matrix_alpha_inv_det(i  ,j)
      enddo
    enddo

  endif


END_PROVIDER


BEGIN_PROVIDER  [ double precision, det_beta_grad_lapl_curr, (4,elec_alpha_num+1:elec_num) ]
  implicit none
  BEGIN_DOC
  ! Gradient and Laplacian of the current beta determinant
  END_DOC

  integer                        :: i, j, k, l

  !DIR$ VECTOR ALIGNED
  do i=elec_alpha_num+1,elec_num
    det_beta_grad_lapl_curr(1,i) = 0.d0
    det_beta_grad_lapl_curr(2,i) = 0.d0
    det_beta_grad_lapl_curr(3,i) = 0.d0
    det_beta_grad_lapl_curr(4,i) = 0.d0
  enddo

  integer                          :: imo, imo2

! -------
! The following code does the same as this:
!
!  do j=1,elec_beta_num
!    imo = mo_list_beta_curr(j)
!    do i=elec_alpha_num+1,elec_num
!      do k=1,4
!        det_beta_grad_lapl_curr(k,i) = det_beta_grad_lapl_curr(k,i) +&
!            mo_grad_lapl_alpha(k,i,imo)*slater_matrix_beta_inv_det(i-elec_alpha_num,j)
!      enddo
!    enddo
!  enddo
!
! --------

  if (iand(elec_beta_num,1) == 0) then

    do j=1,elec_beta_num,2
      imo  = mo_list_beta_curr(j  )
      imo2 = mo_list_beta_curr(j+1)
      !DIR$ LOOP COUNT (16)
      do i=elec_alpha_num+1,elec_num,2
        l = i-elec_alpha_num
        !DIR$ VECTOR ALIGNED
        do k=1,4
          det_beta_grad_lapl_curr(k,i) = det_beta_grad_lapl_curr(k,i) +&
              mo_grad_lapl_beta(k,i,imo )*slater_matrix_beta_inv_det(l,j  ) + &
              mo_grad_lapl_beta(k,i,imo2)*slater_matrix_beta_inv_det(l,j+1)
          det_beta_grad_lapl_curr(k,i+1) = det_beta_grad_lapl_curr(k,i+1) +&
              mo_grad_lapl_beta(k,i+1,imo )*slater_matrix_beta_inv_det(l+1,j  ) + &
              mo_grad_lapl_beta(k,i+1,imo2)*slater_matrix_beta_inv_det(l+1,j+1)
        enddo
      enddo
    enddo

  else

    do j=1,elec_beta_num-1,2
      imo  = mo_list_beta_curr(j  )
      imo2 = mo_list_beta_curr(j+1)
      !DIR$ LOOP COUNT (16)
      do i=elec_alpha_num+1,elec_num-1,2
        l = i-elec_alpha_num
        !DIR$ VECTOR ALIGNED
        do k=1,4
          det_beta_grad_lapl_curr(k,i) = det_beta_grad_lapl_curr(k,i) +&
              mo_grad_lapl_beta(k,i,imo )*slater_matrix_beta_inv_det(l,j  ) + &
              mo_grad_lapl_beta(k,i,imo2)*slater_matrix_beta_inv_det(l,j+1)
          det_beta_grad_lapl_curr(k,i+1) = det_beta_grad_lapl_curr(k,i+1) +&
              mo_grad_lapl_beta(k,i+1,imo )*slater_matrix_beta_inv_det(l+1,j  ) + &
              mo_grad_lapl_beta(k,i+1,imo2)*slater_matrix_beta_inv_det(l+1,j+1)
        enddo
      enddo
      i = elec_num
      l = elec_num-elec_alpha_num
      !DIR$ VECTOR ALIGNED
      do k=1,4
        det_beta_grad_lapl_curr(k,i) = det_beta_grad_lapl_curr(k,i) +&
            mo_grad_lapl_beta(k,i,imo )*slater_matrix_beta_inv_det(l,j  ) + &
            mo_grad_lapl_beta(k,i,imo2)*slater_matrix_beta_inv_det(l,j+1)
      enddo
    enddo

    j = elec_beta_num
    imo = mo_list_beta_curr(j)
    do i=elec_alpha_num+1,elec_num
      l = i-elec_alpha_num
      !DIR$ VECTOR ALIGNED
      do k=1,4
        det_beta_grad_lapl_curr(k,i) = det_beta_grad_lapl_curr(k,i) +&
            mo_grad_lapl_beta(k,i,imo)*slater_matrix_beta_inv_det(l,j)
      enddo
    enddo

  endif

END_PROVIDER


 BEGIN_PROVIDER [ double precision, single_det_value ]
&BEGIN_PROVIDER [ double precision, single_det_grad, (elec_num_8,3) ]
&BEGIN_PROVIDER [ double precision, single_det_lapl, (elec_num) ]
  BEGIN_DOC
  ! Value of a single determinant wave function from the 1st determinant
  END_DOC
  det_i = 1
  det_j = 1
  integer                        :: i
  single_det_value = det_alpha_value_curr * det_beta_value_curr
  do i=1,elec_alpha_num
    single_det_grad(i,1)  = det_alpha_grad_lapl_curr(1,i) * det_beta_value_curr
    single_det_grad(i,2)  = det_alpha_grad_lapl_curr(2,i) * det_beta_value_curr
    single_det_grad(i,3)  = det_alpha_grad_lapl_curr(3,i) * det_beta_value_curr
    single_det_lapl(i)    = det_alpha_grad_lapl_curr(4,i) * det_beta_value_curr
  enddo
  do i=elec_alpha_num+1,elec_num
    single_det_grad(i,1)  = det_alpha_value_curr * det_beta_grad_lapl_curr(1,i)
    single_det_grad(i,2)  = det_alpha_value_curr * det_beta_grad_lapl_curr(2,i)
    single_det_grad(i,3)  = det_alpha_value_curr * det_beta_grad_lapl_curr(3,i)
    single_det_lapl(i)    = det_alpha_value_curr * det_beta_grad_lapl_curr(4,i)
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, psidet_lapl ]
  implicit none
  BEGIN_DOC
  ! Laplacian of the wave functionwithout Jastrow
  END_DOC

  integer                        :: i, j
  psidet_lapl = 0.d0
  do j=1,elec_num
    psidet_lapl = psidet_lapl + psidet_grad_lapl(4,j)
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, E_kin_elec_psidet, (elec_num) ]
  implicit none
  BEGIN_DOC
  ! Electronic Kinetic energy of the determinantal part only: -1/2 (Lapl.Psidet)/Psidet
  END_DOC
  integer                        :: i
  do i=1,elec_num
    E_kin_elec_psidet(i) = -0.5d0*psidet_grad_lapl(4,i) * psidet_inv
  enddo
END_PROVIDER
