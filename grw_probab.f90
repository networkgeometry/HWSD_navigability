!   ---------------------------------------------------------------------------
!   Stochastic greedy routing with time-out
!   12/06/2024                                  
!   Laia Barjuan Ballabriga
!   Optimal navigability of weighted human brain connectomes in physical space
!   ---------------------------------------------------------------------------

PROGRAM GRW 
use mtmod

implicit double precision (a-h,o-z)
integer shortestpath              ! Function computes shortest path
real*8  wvsdprob                  ! Decision function
real*8  deucl                     ! Function computes distance in 3D

integer, parameter :: nodesmax=3000             ! Max number of nodes (list dimension)
integer, parameter :: nedgesmax=50000           ! Max number of edges (list dimension)
integer, parameter :: ndegreemax=500            ! Max degree of a node
integer, parameter :: diameter=5000             ! Diameter of network (longest shortest path allowed)
integer, parameter :: tcontrol=1000             ! Max number of hops in navigation
integer, parameter :: numlambda=25
integer, parameter :: seed=300

real*8  lambda                    ! Parameter for combination weights and distances
integer a                         
real*8  lambdalist(numlambda)     ! List of lambda values
integer nmax                      ! Number of nodes
integer edgemax                   ! Number of edges
integer tar                       ! Target of GR
integer nactual                   ! Current position during GR
integer next                      ! Next position in GR
integer nhops                     ! Number of hops in GR
real*8  probsum                   ! Cumulative value of the probability of transition for each neighbor
real*8  prob                      ! Value of the decision function between actual node and evaluated neighbor
real*8  ctransd                   ! Transmission cost in distance
real*8  ctransw                   ! Transmission cost in weight
real*8  ctransdit                 ! Transmission cost in distance between node i and t
real*8  ctranswit                 ! Transmission cost in weight between node i and t
real*8  totalstretch              ! Total stretch of all paths

real*8  xpos(nodesmax,3)          ! list of x y and z positions
integer edges(nedgesmax,2)        ! list of edges
real*8  weights(nedgesmax)        ! List of weights of edges
real*8  weidist(nedgesmax)        ! List of weight-distances
integer nk(nodesmax)              ! list of degrees
integer npos(nodesmax,2)          ! 1. initial position in neihbors list 2.(after the loop=)final position
integer neighbors(nedgesmax)      ! list of neighbors n11 n12 n13 n14 |n21 n22 n23 n24 n25| n31 n32 ... 
real*8  neighweight(nedgesmax)    ! list of weights of the edges connecting neighbors. Correspondence with the list neighbors(nedgesmax) 
real*8  neighdist(nedgesmax)      ! list of distances of the edges connecting neighbors. Correspondence with the list neighbors(nedgesmax) 
real*8, allocatable :: eucldist(:,:)        ! matrix of euclidean distances between nodes
real*8, allocatable :: probtransition(:,:)  ! matrix of transition probabilities
real*8  probsumlist(ndegreemax)             ! cumulative probability of transition to each neighbor from 0 to 1 (x,x+y,x+y+...+z=1) 
character(50)   filenet,filepos,fileout

filenet="Human0_1015_edges.txt"
filepos="Human0_1015_coord_eucl.txt"
fileout="r_t1000_hum0.txt"

call sgrnd(seed)

!List of lambda values
lambdalist(1)=0.0d0
do i=2,17
  a=-20+i
  lambdalist(i)=dexp(dble(a)/2.0d0)
enddo
do i=18,numlambda
  a=i-15
  lambdalist(i)=dble(a)/10.0d0
enddo

! Initialize lists
xpos=0.0d0
edges=0
weights=0.0d0
weidist=0.0d0
nk=0
npos=0
neighbors=0
neighweight=0.0d0
neighdist=0.0d0

! ====================================================================
! Read positions and edges from file
! ====================================================================
open(13,file=filepos,status='old')
nmax=0
do 
  read(13,*,END=10) node,x,y,z    ! Nodes ordered from 1 to nmax
  xpos(node,1)=x
  xpos(node,2)=y
  xpos(node,3)=z
  if (node>nmax) nmax=node                
enddo
10 continue
close(13)

open(14,file=filenet,status='old')  ! Edges cannot be repeated
edgemax=0                           
do 
  read(14,*,END=11) i,j,weight 
  edgemax=edgemax+1
  edges(edgemax,1)=i
  edges(edgemax,2)=j
  weights(edgemax)=weight
  weidist(edgemax)=-dlog(weight)
  nk(i)=nk(i)+1  
  nk(j)=nk(j)+1
enddo
11 continue
close(14)

! Compute euclidean distances and store in distance matrix 
allocate(eucldist(nmax,nmax))
eucldist=0.0d0                         
do i=1,nmax-1
  do j=i+1,nmax
    dist=deucl(xpos(i,1),xpos(i,2),xpos(i,3),xpos(j,1),xpos(j,2),xpos(j,3))
    eucldist(i,j)=dist
    eucldist(j,i)=dist
  enddo
enddo

! ====================================================================
! Neighbors list
! ====================================================================
! Construct list of initial pointers for neighbors list
ncount=1
do i=1,nmax 
  if (nk(i).ne.0) then
    npos(i,1)=ncount   !initial position in neighbors list
    npos(i,2)=npos(i,1)-1  
    ncount=ncount+nk(i) 
  endif 
enddo

! Neighbors list and final pointer for neighbors list
do i=1,edgemax
  node1=edges(i,1)
  node2=edges(i,2)
  npos(node1,2)=npos(node1,2)+1     !actualize final position in neighbors list for node 1
  npos(node2,2)=npos(node2,2)+1     
  neighbors(npos(node1,2))=node2    !write n2 as neighbor of node 1
  neighbors(npos(node2,2))=node1   
  neighweight(npos(node1,2))=weights(i)
  neighweight(npos(node2,2))=weights(i)
  neighdist(npos(node1,2))=weidist(i)
  neighdist(npos(node2,2))=weidist(i)
enddo


! ====================================================================
! Greedy routing weighted
! ====================================================================
open(15,file=fileout)
write(15,*) "#lambda,  Success rate,  Average stretch, Ctransd,  Ctransw"

allocate(probtransition(nmax,nedgesmax))
do a=1,numlambda
  lambda=lambdalist(a)

  ! Compute the transition probability from each node to its neighbors when going to a given target
  probtransition=0.0d0
  do i=1,nmax
    do tar=1,nmax
      do k=npos(i,1),npos(i,2)
        probtransition(tar,k)=wvsdprob(k,i,tar,npos,neighbors,neighdist,eucldist,nodesmax,nedgesmax,nmax,lambda)
      enddo
    enddo
  enddo

  ctransd=0.0d0       !Information cost Euclidean distances
  ctransw=0.0d0       !Information cost weight distances
  totalstretch=0.0d0  !Stretch
  nfail=0             !Failed paths

  do i=1,nmax         !Source
  if (nk(i).ne.0) then 
    do tar=1,nmax     !Target
    if (nk(tar).ne.0) then 
    if (i.ne.tar) then
      nhops=0
      ctranswit=0.0d0
      ctransdit=0.0d0
      nactual=i       !Current position
      do while (nactual.ne.tar)
        probsum=0.0d0
        probsumlist=0.0d0
        nneighbor=0
        ! Create list of cumulative jump probability and actualize information cost
        do k=npos(nactual,1),npos(nactual,2)    !For each neighbor of the current node
          nneighbor=nneighbor+1
          prob=probtransition(tar,k)
          ctranswit=ctranswit+prob*(neighdist(k))
          ctransdit=ctransdit+prob*(eucldist(neighbors(k),nactual))
          probsum=probsum+prob
          probsumlist(nneighbor)=probsum
        enddo    
        ! Select to which neighbor jump next
        randomnum=grnd()
        nneighbor=0
        do k=npos(nactual,1),npos(nactual,2) 
          nneighbor=nneighbor+1
          if (randomnum<probsumlist(nneighbor)) then
            next=neighbors(k)
            go to 50
          endif
        enddo
  50    continue
        nhops=nhops+1
        nactual=next !Jump
        if (nhops.gt.tcontrol) then !Time-out control
          nfail=nfail+1
          go to 60
        endif
      enddo
      totalstretch=totalstretch+dble(nhops)/dble(shortestpath(i,tar,npos,neighbors,nodesmax,nedgesmax,diameter))
      ctransw=ctransw+ctranswit
      ctransd=ctransd+ctransdit
    endif
    endif
60  continue
    enddo
  endif
  enddo
  nsuccess=nmax*(nmax-1)-nfail
  write(15,*) lambda,nsuccess/dble(nmax*(nmax-1)),totalstretch/dble(nsuccess),ctransd/dble(nsuccess),ctransw/dble(nsuccess)
enddo

close(15)
deallocate(eucldist)
deallocate(probtransition)

END PROGRAM GRW


! ====================================================================
! Functions
! ====================================================================
! DECISION FUNCTION --------------------------------------------------
! Computes function that balances weights and Euclidean distance between connected nodes
real*8 function wvsdprob(knode,nactual,tar,npos,neighbors,ndist,eucldist,nodesmax,nedgesmax,nmax,lam)
  implicit none
  !INPUT
  integer knode                  ! Place in list of neighbor currently evaluated
  integer nactual                ! Current node in greedy routing
  integer tar                    ! Target node
  integer nmax                   ! Number of nodes
  integer nodesmax               ! Max number of nodes (list dimension)
  integer nedgesmax              ! Max number of edges (list dimension)
  integer npos(nodesmax,2)       ! 1. initial position in neihbors list 2.final position
  integer neighbors(nedgesmax)   ! list of neighbors n11 n12 n13 n14 |n21 n22 n23 n24 n25| n31 n32 ... 
  real*8  ndist(nedgesmax)       ! list of wedistances of the edges connecting neighbors. Correspondence with the list neighbors(nedgesmax) 
  real*8  eucldist(nmax,nmax)      
  !OTHER VARIABLES
  integer i 
  real*8  sum
  real*8  lam
  real*8  deuclktar
  real*8  dk, di

  dk=ndist(knode)
  deuclktar=eucldist(neighbors(knode),tar)
  sum=0.0d0
  do i=npos(nactual,1),npos(nactual,2)
    di=ndist(i)
    sum=sum+dexp(-(lam*(eucldist(neighbors(i),tar)-deuclktar)+(1.0d0-lam)*(di-dk)))
  enddo
  wvsdprob=1.0d0/sum
  return
end function wvsdprob

! EUCLIDEAN DISTANCE 3D-------------------------------------------------
real*8 function deucl(x1,y1,z1,x2,y2,z2)
  implicit double precision(x,y,z)
  deucl=dsqrt((x1-x2)**2.0d0+(y1-y2)**2.0d0+(z1-z2)**2.0d0)
  return
end function deucl

! --------------------------------------------------------------------
! SHORTEST PATH 
! --------------------------------------------------------------------
integer function shortestpath(source,tar,npos,neighbors,nodesmax,nedgesmax,diameter)
  implicit double precision(a-h,o-z)
  ! INPUT
  integer  nodesmax               ! Max number of nodes (list dimension)
  integer  nedgesmax              ! Max number of edges (list dimension)
  integer  diameter               ! Max number of steps in a path
  integer  source, tar            ! source and target 
  integer  npos(nodesmax,2)       ! 1. initial position in neihbors list 2.final position
  integer  neighbors(nedgesmax)   ! list of neighbors n11 n12 n13 n14 |n21 n22 n23 n24 n25| n31 n32 ... 
  ! OTHER VARIABLES
  integer  nvisited(0:diameter)          ! List of nodes already taken into account
  integer  nallowed(nodesmax)            ! nallowed(i)=1 node "i" has already been visited. If =0 still not visited
  integer  nodes_at_distance(0:diameter) ! nodes_at_distance(i)=x "x" is the number of nodes at topological distance "i" from node source
  integer  ncount                        ! counter, # of visited node in the sequence
  integer  l                             ! counter, length of path
  integer  nactual                       ! current node
  integer  neigh                         ! neighbor currently evaluated

  nvisited=0
  nallowed=0
  nodes_at_distance=0

  !Start from source
  nvisited(0)=source ! Include the source position as 0
  nallowed(source)=1 ! exclude own node
  nodes_at_distance(0)=1 !Include the source as l=0

  ncount=0
  do l=1,diameter
    ncountold=ncount
    do i=1+ncountold-nodes_at_distance(l-1),ncountold  !position in sequence of nodes at path length (l-1)
      nactual=nvisited(i)
      do k=npos(nactual,1),npos(nactual,2)
        neigh=neighbors(k)
        if (neigh.eq.tar) then
          shortestpath=l
          go to 70
        endif
        if (nallowed(neigh).eq.0) then
          ncount=ncount+1
          nvisited(ncount)=neigh  ! at step "ncount" I visited neighbor neigh
          nallowed(neigh)=1       ! not allowed bc already visited
          nodes_at_distance(l)=nodes_at_distance(l)+1
        else
        endif
      enddo
    enddo
    if (nodes_at_distance(l).eq.0) then !If not more nodes can be reached
      shortestpath=0
      go to 70  
    endif
    if (l.eq.diameter) write(*,*) "Max diameter too small" !If not reached node after #diam hops
  enddo
70 continue
  return
end function shortestpath
