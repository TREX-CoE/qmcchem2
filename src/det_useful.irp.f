
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

  integer(8)              :: xorvec(256)
  integer                        :: l

  if (Nint > 256) then
    print *, 'in det_useful.irp.f, get_excitation_degree_spin, Nint too lage'
    stop 1
  endif
  ASSERT (Nint > 0)

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

! ---

 BEGIN_PROVIDER [ integer, det_alpha_num_pseudo ]
&BEGIN_PROVIDER [ integer, det_beta_num_pseudo ]

  BEGIN_DOC
  ! Dimensioning of large arrays made smaller without pseudo
  END_DOC

  implicit none

  if(do_pseudo) then
    det_alpha_num_pseudo = det_alpha_num
    det_beta_num_pseudo  = det_beta_num
  else
    det_alpha_num_pseudo = 1
    det_beta_num_pseudo  = 1
  endif

END_PROVIDER

! ---



