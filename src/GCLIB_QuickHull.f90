    !  ͨ�ñ����Fortran���ݽṹ��������  (GCLIB)
    !   General and Convenient Fortran Data Structure and Function Library
    !---------------------------------------------------------------
    !
    ! GCLIB_QuickHull 
    ! 2023.6.15 - coded - XIE JIHONG 
    ! 2023.8.22 - modified - XIE JIHONG 
    ! 2023.9. 7 - openMP optimizated - XIE Jihong
    ! ʹ�÷����� ֱ�ӵ���QuickHull()���ɻ��͹��������Ƭ
    !
    !---------------------------------------------------------------
    MODULE GCLIB_QuickHull
        USE  OMP_LIB
    IMPLICIT NONE  ! Used to encapsulate variables
        
        ! ��
        TYPE, PRIVATE :: Point
            REAL*8 :: xyz(3) = 0.D0
            LOGICAL*1 :: isOn
        END TYPE Point
            
            TYPE, PRIVATE :: Point_ptr
                TYPE(Point), POINTER :: p => NULL() ! 
            END TYPE Point_ptr
        
            TYPE, PRIVATE :: Face_ptr
                TYPE(Face), POINTER :: p => NULL()! 
            END TYPE Face_ptr
        
        ! ��
        TYPE, PRIVATE :: Face 
            REAL*8 :: vertex(3,3) = 0.D0  ! �����
            REAL*8 :: nml(3) = 0.D0 ! ������������ķ�ʸ
            TYPE(Face_ptr) :: N(3) ! ������
            TYPE(Point_ptr), ALLOCATABLE :: OS(:) ! �ⲿ�㼯
            LOGICAL*1 :: flag ! ��־λ
        END TYPE Face

        ! �ٽ��
        TYPE, PRIVATE :: CriticalEdge
            REAL*8 :: ep(2,3)
            TYPE(Face_ptr) :: f(2) ! ��ʱδ�ڳ������õ��������Ż���
        END TYPE CriticalEdge

        PUBLIC :: QuickHull
        
        PRIVATE :: FILT_DUPL
        PRIVATE :: initial_p, Get_Start_Convex, Point_into_Type, Calc_NormalVector_DerivForm_MassPoint
        PRIVATE :: Distr_OS_for_Face, Store_Nonempty_Faces, Update_Adjacent_Faces, Initial_Flag_FaceSet
        PRIVATE :: Search_Farthest_Point_from_FaceSet, Initial_Visible_FaceSet, Store_Face_Into_Visible_FaceSet
        PRIVATE :: find_finished_visible_faceSet, together_OS_to_L, find_criticalEdge_to_H
        PRIVATE :: generate_newFace_from_P_and_criticalEdge, remove_V_add_NS
        
        PRIVATE :: FOOT_PL, FOOT_LL, UNINML
        PRIVATE :: VEC_PL, OVERLAP
        PRIVATE :: DIST_PF_SIGN, UTZVEC, IS_INSIDE_PF, isZero8
        PRIVATE :: isPointInSimplex, CROSS_PRODUCT_3D, SET_DIFF
        PRIVATE :: getHullMeshesVertex_post_debug
        
        !---------------------------------------------------------------
        ! Random data
        PRIVATE :: GET_RANDOM_UNIT_VECTOR
        
        
    CONTAINS  ! Used to encapsulate functions or subroutines
    
    
        SUBROUTINE QuickHull(p_, hullMesh_, info_)
        !DEC$ ATTRIBUTES DLLEXPORT :: QuickHull
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(:,:)
            REAL*8, INTENT(OUT), ALLOCATABLE :: hullMesh_(:,:,:)
            INTEGER*4, INTENT(OUT) :: info_ ! 0 - ͹�����㲢���ֳɹ���  1 - ����ʧ��
            
            TYPE(Face), SAVE :: startHull(4) !��ʼ͹��
            !$OMP THREADPRIVATE(startHull) ! 
            INTEGER*4, SAVE :: istat, i, j
            !$OMP THREADPRIVATE(istat, i, j) ! 
            REAL*8, SAVE :: p_far(3)
            !$OMP THREADPRIVATE(p_far) ! 
            REAL*8, SAVE, ALLOCATABLE :: p_filted(:,:), p_temp(:,:)
            !$OMP THREADPRIVATE(p_filted, p_temp) ! 
            
            TYPE(Face), SAVE, ALLOCATABLE, TARGET :: Hull(:) !͹��
            TYPE(Face_ptr), SAVE, ALLOCATABLE :: Q(:) !�����漯
            TYPE(Face_ptr), SAVE, ALLOCATABLE :: V(:) !�ɼ��漯
            TYPE(Face), SAVE, ALLOCATABLE, TARGET :: NS(:) ! �����ٽ�ߵ����漯��
            TYPE(Point), SAVE, ALLOCATABLE, TARGET :: P(:) ! �㼯
            TYPE(Point_ptr), SAVE, ALLOCATABLE :: L(:) ! V�е��ⲿ�㼯
            TYPE(CriticalEdge), SAVE, ALLOCATABLE :: H(:) ! �ٽ��
            !$OMP THREADPRIVATE(Hull, Q, V, NS, P, L, H) ! 
            
            INTEGER*4, SAVE :: num_thread_now, iterC
            !$OMP THREADPRIVATE(num_thread_now, iterC) ! 
            
            ! ��ʼ��
            iterC = 0 ! ����������
            info_ = 0
            
            ! �����Ҫ������͹�������㷨�������Ƿ�淶
            IF ( (SIZE(p_, 1) < 4) .OR. (SIZE(p_, 2) /= 3) ) THEN
                info_ = 1 ! 1 - ����͹��ʧ��
                RETURN
            END IF
            
            ! ȥ���ظ��ĵ�
            CALL FILT_DUPL(p_, p_filted, info_)
            IF ( info_ /= 0 ) RETURN
            
            ! ȥ�غ��ٴμ����Ҫ������͹�������㷨�������Ƿ�淶
            ! 1���ж��Ƿ�㲻��
            IF ( (SIZE(p_filted, 1) < 4) .OR. (SIZE(p_filted, 2) /= 3) ) THEN
                info_ = 1 ! 1 - ����͹��ʧ��
                RETURN
            END IF
            
            ! 2���ж��Ƿ����е㶼��ͬһƽ���ϣ����Ϊ1.D-12������������������ǣ���ô�˻�Ϊ��ά��,��ǰ�汾���벻����
            info_ = INFO_COP(p_filted)
            IF ( info_ /= 0 ) THEN
                info_ = 1 ! 1 - ����͹��ʧ��
                RETURN
            END IF
            
            ! ��ʼ�����
            CALL initial_p(P)
            
            ! Ϊ�㼯P����ռ�
            ALLOCATE( P( SIZE(p_filted, 1) ), STAT = istat )
            
            ! P��isOn��������
            CALL Initial_P_isOn(P)
            
            ! ��ȡ��ʼ͹��
            CALL Get_Start_Convex(p_filted, startHull, info_)
            IF ( info_ /= 0 ) RETURN
            
            IF(ALLOCATED(Hull)) DEALLOCATE(Hull)
            Hull = startHull
            
            ! �����ṹ�� �ܵ㼯P
            CALL Point_into_Type(p_filted, P)
            
            ! ���´����convex��ÿ����ĳ���ķ�ʸ��������ṹ���ڵ�nml����
            CALL Calc_NormalVector_DerivForm_MassPoint(Hull)
            
            ! ��p���䵽��ʼ͹����ÿ���� F ���ⲿ�㼯OS��
            CALL Distr_OS_for_Face(P, Hull)
            
            ! �ѵ�ǰHull�ⲿ�㼯�ǿյ��汣��������漯Q�У�
            CALL Store_Nonempty_Faces(Hull, Q)

            ! ˢ�µ�ǰ͹�����Լ���Ӧ�Ĵ����漯Q�е� ������
            CALL Update_Adjacent_Faces(Hull)
                    
            DO WHILE(.TRUE.)
                ! ������������ ��������������࣬����Ԥ��
                iterC = iterC + 1 
                IF ( iterC > 500 ) THEN
                      WRITE(UNIT=6, FMT="(A)") "Error warning: too many iterations for GCLIB_QuickHull subroutine." ! 
                      PAUSE
                END IF
                
                ! û���µĴ����漯�ˣ���˵���Ѿ��ҵ�͹��
                IF ( .NOT. ALLOCATED(Q) ) EXIT
                
                ! ���漯������flag ��Ϊ .false.
                CALL Initial_Flag_FaceSet(Hull)
                
                ! �� F ���ⲿ�㼯���ҵ�������F��Զ�ĵ� p
                CALL Search_Farthest_Point_from_FaceSet( Q(1), p_far )
                
                ! ��ʼ���ɼ��漯 V������ F(Q���ĵ�һ����) ����� V ��
                CALL Initial_Visible_FaceSet(V)
                CALL Store_Face_Into_Visible_FaceSet( Q(1), V )
                
                ! Ѱ�� p_far ��Ե���ɵĿɼ��漯 V
                CALL find_finished_visible_faceSet(p_far, V)
                
                ! �Ѽ��� V ��ÿ������ⲿ�㼯���ܵ�һ���㼯 L ��
                !CALL together_OS_to_L(V, P, L)
                
                ! �������� V �����ٽ�ߣ�����һ������ H;
                CALL find_criticalEdge_to_H(p_far, V, H, info_)
                IF ( info_ /= 0 ) RETURN
                
                ! ���ӵ� p ���ٽ�� H�����е�ÿ��R������һ���µ��棬���뵽�������漯NS������δ����NS�����棩
                CALL generate_newFace_from_P_and_criticalEdge(p_far, H, NS)
                
                ! ��͹��hull���Ƴ��ɼ��漯V����������NS��ӽ�ȥ������͹��hull
                CALL remove_V_add_NS( V, NS, Hull)
                !CALL getHullMeshesVertex_post_debug(Hull, info_)
                ! ������͹�����ⷨʸ
                CALL Calc_NormalVector_DerivForm_MassPoint(Hull)
                
                ! P��isOn��������
                CALL Initial_P_isOn(P)
                
                ! ��p���䵽��͹����ÿ���� F ���ⲿ�㼯OS��
                CALL Distr_OS_for_Face(P, Hull)

                ! �ѵ�ǰHull�ⲿ�㼯�ǿյ��汣��������漯Q�У�
                CALL Store_Nonempty_Faces(Hull, Q)
                
                ! ������͹���ĵ�������
                CALL Update_Adjacent_Faces(Hull)
                CALL Initial_Flag_FaceSet(Hull)
                
                ! û���µĴ����漯�ˣ���˵���Ѿ��ҵ�͹��
                IF ( .NOT. ALLOCATED(Q) ) EXIT
                
                ! ���û�У���ôˢ�µ�ǰ͹�����Լ���Ӧ�Ĵ����漯Q�е� ������
                CALL Update_Adjacent_Faces(Hull)
                !pause
            END DO

            ! ���Hull�����Hull�е���������⣬��ô�Ƴ�
            !CALL Hull_list%traverse(removeProblemFace)

            ! �������
            IF(ALLOCATED(hullMesh_)) DEALLOCATE(hullMesh_)
            ALLOCATE( hullMesh_( SIZE(Hull), 3, 3 ), STAT = istat )
            hullMesh_ = 0.D0
            
            DO i = 1, SIZE(Hull)
                hullMesh_(i, :, :) = Hull(i)%vertex(:,:)
            END DO
            
            RETURN      
        END SUBROUTINE QuickHull
        
        !---------------------------------------------------------------
        ! ���Hull�����Hull�е���������⣬��ô�Ƴ�
        !SUBROUTINE removeProblemFace(Hull_)
        !    IMPLICIT NONE
        !    TYPE(Face), ALLOCATABLE, INTENT(INOUT) :: Hull_(:) !͹��
        !    
        !    IF( ALL( DABS(p(1,:) - p(2,:)) < 1.D-12 )  .OR. &
        !        ALL( DABS(p(2,:) - p(3,:)) < 1.D-12 )  .OR. &
        !        ALL( DABS(p(3,:) - p(1,:)) < 1.D-12 ) ) THEN
        !        CALL Hull_list%remove( Hull_list%search(p) )
        !        Alert = 1
        !    END IF
        !        
        !    RETURN
        !END SUBROUTINE removeProblemFace
                          
        !---------------------------------------------------------------
        !
        ! ����P�е�isOn����
        !
        !---------------------------------------------------------------
        SUBROUTINE Initial_P_isOn(P_)
            IMPLICIT NONE
            TYPE(Point), INTENT(INOUT) :: P_(:) ! 
            INTEGER*4 :: i ! 
            DO i = 1, SIZE(P_), 1
                P_(i)%isOn = .False.
            END DO
            RETURN
        END SUBROUTINE Initial_P_isOn
        
        !---------------------------------------------------------------
        !
        ! �ж�����ĵ㼯�Ƿ�ȫ������,��������4��������
        ! info_ : 0 - ������   1 - �����ҹ���   2 - ���治����  3 - ����   4 - ���棬���в��ֵ��غ�
        !---------------------------------------------------------------
        FUNCTION INFO_COP(p_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(:,:)
            INTEGER*4 :: res_ !  
            
            INTEGER*4 :: i ! 
            REAL*8, ALLOCATABLE :: p_filted(:,:)
            REAL*8 :: nml(3), dot
            REAL*8 :: refeFace(3,3)
            REAL*8 :: AB(3), BC(3), cross(3)
            
            ! check data 
            IF ( SIZE(p_,1) < 4 .OR. SIZE(p_, 2) /= 3 ) THEN
                WRITE(UNIT=6, FMT="(A)") "FUNCTION res_COP(p_, res_)" ! 
                READ(UNIT=5, FMT=*)   
                STOP
            END IF
            
            ! ȥ���ظ��ĵ�
            CALL FILT_DUPL(p_, p_filted, res_)
            IF ( res_ /= 0 ) RETURN ! ���ù�
            
            ! �ٴ�check data
            IF ( (SIZE(p_filted, 1) < 4) .OR. (SIZE(p_filted, 2) /= 3) ) THEN
                IF ( SIZE(p_filted, 1) == 1 ) THEN
                    res_ = 3 ! ����
                    RETURN
                ELSE
                    res_ = 4 ! ���棬���в��ֵ��غ�
                    RETURN
                END IF
            END IF
            
            ! ��ѡ3������Ϊ��׼��
            refeFace = 0.D0
                ! ��ѡ������
            refeFace(1, :) = p_filted(1,:)
            refeFace(2, :) = p_filted(2,:)
            AB(:) = refeFace(2, :) - refeFace(1, :)
                ! ��ѡ��3���㣬Ҫ������㲻�ܸ�ǰ�����㹲��
            DO i = 3, SIZE(p_filted, 1), 1
                BC(:) = p_filted(i,:) - refeFace(2, :)
                
                ! ֱ�� AB X BC /= 0 Ϊֹ                
                cross(1) = AB(2) * BC(3) - AB(3) * BC(2)
                cross(2) = AB(3) * BC(1) - AB(1) * BC(3)
                cross(3) = AB(1) * BC(2) - AB(2) * BC(1)
                !cross = CROSS_PRODUCT_3D(AB,BC)
                IF( ALL( DABS( cross ) < 1.D-12 ) ) THEN
                    ! ֱ�����Ҳû���ҵ������ߵĵ㣬��ô���Է��ؽ����
                    IF ( i == SIZE(p_filted, 1) ) THEN
                        res_ = 1 ! �����ҹ���
                        RETURN
                    END IF
                    CYCLE
                ELSE
                    refeFace(3, :) = p_filted(i, :)
                    EXIT
                END IF
            END DO

            ! �ж������ĵ��Ƿ�������湲��
            nml = cross / NORM2(cross)
            DO i = 4, SIZE(p_filted, 1), 1
                dot = DOT_PRODUCT(nml, p_filted(i,:) - p_filted(1,:))
                IF( DABS( dot ) < 1.D-12 ) THEN
                    ! ֱ�����Ҳû���ҵ�������ĵ㣬��ô���Է��ؽ����
                    IF ( i == SIZE(p_filted, 1) ) THEN
                        res_ = 2 ! ���浫������
                        RETURN
                    END IF
                    CYCLE
                ELSE
                    ! ��������һ�������׼�治���棬��ô���ز�����
                    res_ = 0
                    RETURN
                END IF
            END DO

            RETURN
        END FUNCTION INFO_COP
        
        !---------------------------------------------------------------
        !
        ! �˵��ظ��ĵ�
        !
        !---------------------------------------------------------------
        SUBROUTINE FILT_DUPL(p_in_, p_out_, info_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_in_(:,:)
            REAL*8, INTENT(OUT), ALLOCATABLE :: p_out_(:,:)
            INTEGER*4, INTENT(OUT) :: info_! 
            
            INTEGER*4 :: i, j, k, C, istat
            REAL*8, ALLOCATABLE :: v(:,:), v_temp(:,:)
            REAL*8 :: thisMesh(3,3), subs(3)
            LOGICAL*1 :: isOverloap ! 
            
            ! �����p_in_���淶���򷵻ش�����Ϣ
            IF ( SIZE(p_in_, 2) /= 3 ) THEN
                info_ = 1 ! 1 - �������鲻�淶
                RETURN
            END IF
            
            ! ��������
            C = 1
            ALLOCATE( v(1,3), STAT = istat )
            v(1,:) = p_in_(1,:)
                
            ! �����㣬�ж�ÿ�����Ƿ���Ҫ����ӽ�v
            DO j = 1, SIZE(p_in_, 1), 1
                    
                ! �ǵ�ȥ��v��ÿ��������ж�
                isOverloap = .FALSE.
                k = 1
                DO WHILE( k <= SIZE(v, 1) )
                    subs(:) = p_in_(j,:) - v(k,:)
                    ! ��������ص�����ǣ�������WHILEѭ��
                    IF ( ALL( DABS(subs(:) ) < 1.D-12) ) THEN
                        isOverloap = .TRUE.           
                        EXIT
                    END IF
                    k = k + 1
                END DO
                    
                ! ������ص�������Ϊ���µĵ㣬�����ȥv
                IF ( isOverloap == .FALSE. ) THEN
                    IF( ALLOCATED(v_temp) ) DEALLOCATE(v_temp)
                    v_temp = v
                    DEALLOCATE(v)
                    ALLOCATE( v( 1 + SIZE(v_temp, 1), 3), STAT = istat )
                    v( 1 : SIZE(v_temp, 1), : ) = v_temp(:,:)
                    v( 1 + SIZE(v_temp, 1) , : ) = p_in_(j,:)
                END IF
                    
            END DO
            
            p_out_ = v
            info_ = 0
            
            RETURN
        END SUBROUTINE FILT_DUPL
        
        !---------------------------------------------------------------
        !
        ! ��͹��hull���Ƴ��ɼ��漯V����������NS��ӽ�ȥ������͹��hull
        !
        !---------------------------------------------------------------
        SUBROUTINE remove_V_add_NS( V_, NS_, Hull_)
            IMPLICIT NONE
            TYPE(Face_ptr), INTENT(IN) :: V_(:) !�ɼ��漯
            TYPE(Face), INTENT(IN) :: NS_(:) ! �����ٽ�ߵ����漯��
            
            TYPE(Face), ALLOCATABLE, INTENT(INOUT) :: Hull_(:) !͹��

            INTEGER*4 :: i, j, istat, index
            INTEGER*4, ALLOCATABLE :: beingRemvIdx(:), keptIdx(:), iota(:)
            TYPE(Face), ALLOCATABLE :: Hull_new(:) !͹��
            
            ! �ڴ洦��
            IF(ALLOCATED(beingRemvIdx)) DEALLOCATE(beingRemvIdx)
            ALLOCATE( beingRemvIdx( SIZE(V_) ) )
            IF(ALLOCATED(iota)) DEALLOCATE(iota)
            ALLOCATE( iota( SIZE(Hull_) ) )
            
            ! ��¼
            DO i = 1, SIZE(V_) 
                beingRemvIdx(i) = locateBeingRemvIndex(Hull_, V_, i)
            END DO
            
            iota = [(i, i = 1, SIZE(Hull_))]
            keptIdx = SET_DIFF(iota, beingRemvIdx)
            
            ! Ϊ��͹�������ڴ�
            IF(ALLOCATED(Hull_new)) DEALLOCATE(Hull_new)
            ALLOCATE( Hull_new( SIZE(keptIdx) + SIZE(NS_) ) )
            
            ! �õ�һ����͹��
            DO i = 1, SIZE(keptIdx) 
                Hull_new(i) = Hull_( keptIdx(i) )
            END DO
            
            DO i = 1 , SIZE(NS_) 
                 Hull_new( SIZE(keptIdx) + i) = NS_(i)
            END DO

            IF(ALLOCATED(Hull_)) DEALLOCATE(Hull_)
            !ALLOCATE( Hull_( SIZE(Hull_new) ) )
            
            Hull_ = Hull_new
            
            RETURN
            CONTAINS
                ! ����hull����Ҫ���Ƴ�����ı�ţ���Ѱ��ӦV__�е�i__����
                PURE FUNCTION locateBeingRemvIndex(Hull__, V__, idx__) RESULT(res__)
                    TYPE(Face), ALLOCATABLE, INTENT(IN) :: Hull__(:) !͹��
                    TYPE(Face_ptr), INTENT(IN) :: V__(:) !�ɼ��漯
                    INTEGER*4, INTENT(IN) :: idx__
                    INTEGER*4 :: res__, i
                    REAL*8 :: debug(3,3)
                    
                    DO i = 1, SIZE(Hull__) 
                        IF( ALL( DABS(Hull__(i)%vertex(:,:) - V__(idx__)%p%vertex(:,:)) < 1.D-8) ) THEN
                            res__ = i
                            EXIT
                        END IF
                    END DO
                    
                    RETURN
                END FUNCTION locateBeingRemvIndex
        END SUBROUTINE remove_V_add_NS
  
        !---------------------------------------------------------------
        !
        ! ���ӵ� p ���ٽ�� R������һ���µ��棬���뵽�������漯NS
        !
        !---------------------------------------------------------------
        SUBROUTINE generate_newFace_from_P_and_criticalEdge(p_, H_, NS_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(3)
            TYPE(CriticalEdge), INTENT(IN) :: H_(:) ! �ٽ��
            TYPE(Face), ALLOCATABLE, TARGET, INTENT(INOUT) :: NS_(:) ! �����ٽ�ߵ����漯��
            
            INTEGER*4 :: i , j, istat
            LOGICAL*1 :: L ! 
            ! ��NS����ʼ��
            IF ( ALLOCATED(NS_) ) THEN
                DO j = 1, SIZE( NS_ ), 1
                    NULLIFY( NS_(j)%N(1)%p, NS_(j)%N(2)%p, NS_(j)%N(3)%p) ! ����ָ�룬��������ָ��
                    IF ( ALLOCATED(NS_(j)%OS) ) THEN
                        DO i = 1, SIZE( NS_(j)%OS ), 1
                            NS_(j)%OS(i)%p => NULL()  ! ����ָ�룬��������ָ��
                        END DO 
                    END IF
                END DO
                DEALLOCATE( NS_ ) ! ���ٶ���
            END IF
            
            ! �ж������ٽ�ߣ����ж��ٸ�����
            ALLOCATE( NS_( SIZE( H_ ) ), STAT = istat ) 
            
             ! ����p���ٽ�ߣ�
            DO i = 1, SIZE( H_ ), 1
                NS_(i)%vertex(1,:) = p_ ! ����ĵ�һ������Ϊ�����p
                NS_(i)%vertex(2,:) = H_(i)%ep(1,:) ! ������������Ϊ�ٽ�ߵ������˵�
                NS_(i)%vertex(3,:) = H_(i)%ep(2,:) ! ������������Ϊ�ٽ�ߵ������˵�
            END DO
            
            RETURN
            CONTAINS
                ! �жϵ��Ƿ����߶���
                PURE FUNCTION isPointOnEdge(p_, edge_) RESULT(res_)
                    REAL*8, INTENT(IN) :: p_(3), edge_(2,3)
                    LOGICAL*1 :: res_ ! 
                    REAL*8 :: PA(3), PB(3)
                    PA = edge_(1,:) - p_
                    PB = edge_(2,:) - p_
                    IF ( (ABS(NORM2(CROSS_PRODUCT_3D(PA, PB) ) ) < 1.D-12) .AND. (DOT_PRODUCT(PA, PB) <= 1.D-12 ) ) THEN
                        res_ = .TRUE.
                    ELSE    
                        res_ = .FALSE.
                    END IF
                    RETURN
                END FUNCTION

        END SUBROUTINE generate_newFace_from_P_and_criticalEdge
        
        
        !---------------------------------------------------------------
        !
        ! �������� V �����ٽ�ߣ�����һ������ H;
        !
        !---------------------------------------------------------------
        SUBROUTINE find_criticalEdge_to_H(p_, V_, H_, info_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(3)
            TYPE(Face_ptr), INTENT(IN) :: V_(:) !�ɼ��漯
            TYPE(CriticalEdge), ALLOCATABLE, INTENT(INOUT) :: H_(:) ! �ٽ��
            INTEGER*4, INTENT(OUT) :: info_ !
            
            INTEGER*4 :: i, j, k, C, istat
            REAL*8 :: common_side(2,3)
            TYPE(CriticalEdge), ALLOCATABLE :: H_temp(:) ! �ٽ��temp
            
            ! ������Ϣ
            info_ = 0
            
            ! ��ʼ��H
            IF ( ALLOCATED( H_ )) THEN
                FORALL (k = 1 : SIZE( H_ ) ) H_(k)%f(1)%p => NULL() ! ��ָ��
                FORALL (k = 1 : SIZE( H_ ) ) H_(k)%f(2)%p => NULL() ! ��ָ��
                DEALLOCATE(H_) ! ��ɾ��
            END IF
            
            C = 0
            DO i = 1, SIZE(V_), 1
            DO j = 1, 3, 1 ! ����ÿ���ɼ����3������
                ! ������治�ɼ�����ô�ҵ��ٽ��
                IF ( DOT_PRODUCT( V_(i)%p%N(j)%p%nml(:), p_ - V_(i)%p%N(j)%p%vertex(1,:) ) < 1.D-12 ) THEN ! ����<0���򱨴�
                    C = C + 1
                        
                    ! �����ٽ��
                    CALL get_common_TriaSide( V_(i)%p%N(j)%p%vertex(:,:), V_(i)%p%vertex(:,:), common_side, info_  ) 
                        
                    ! ����H���ep
                    IF ( .NOT. ALLOCATED( H_ ) ) THEN
                        ALLOCATE( H_(1), STAT = istat )
                        H_(1)%ep(:,:) = common_side(:,:)
                    ELSE
                        ! ��ʼ��H_temp
                        IF ( ALLOCATED( H_temp ) ) THEN
                            FORALL (k = 1 : SIZE( H_temp ) ) H_temp(k)%f(1)%p => NULL() ! ��ָ��
                            FORALL (k = 1 : SIZE( H_temp ) ) H_temp(k)%f(2)%p => NULL() ! ��ָ��
                            DEALLOCATE(H_temp) ! ��ɾ��
                        END IF
                        
                        ! �ݴ�
                        H_temp = H_
                        
                        ! ��ʼ��H
                        IF ( ALLOCATED( H_ )) THEN
                            FORALL (k = 1 : SIZE( H_ ) ) H_(k)%f(1)%p => NULL() ! ��ָ��
                            FORALL (k = 1 : SIZE( H_ ) ) H_(k)%f(2)%p => NULL() ! ��ָ��
                            DEALLOCATE(H_) ! ��ɾ��
                        END IF
                        ALLOCATE( H_( 1 + SIZE(H_temp) ), STAT = istat ) ! �ռ� + 1
                        
                        ! �µ�H
                        H_(1 : C - 1) = H_temp(:)
                        H_(C)%ep(:,:) = common_side(:,:)
                    END IF

                END IF
            END DO
            END DO
            
            RETURN
            
            CONTAINS
                !��ȡ�����ռ������εĹ�����
                SUBROUTINE get_common_TriaSide(tr1_, tr2_, res_, info_) 
                IMPLICIT NONE
                    REAL*8, INTENT(IN) :: tr1_(3,3), tr2_(3,3)
                    REAL*8, INTENT(OUT) :: res_(2,3)
                    INTEGER*4, INTENT(OUT) :: info_ ! 
                    REAL*8 :: eg1(2,3), eg2(2,3), sub1(3), sub2(3), sub3(3), sub4(3)
                    INTEGER*4 :: i, j, k
                    
                    ! ������Ϣ
                    info_ = 0
                    
                    ! tr1��3����
                    DO i = 1, 3, 1
                        IF ( i == 3 ) THEN
                            k = 1
                        ELSE
                            k = i + 1
                        END IF
                        eg1(1,:) = tr1_(i,:)
                        eg1(2,:) = tr1_(k,:)
                        
                        ! tr2��3����
                        DO j = 1, 3, 1
                            IF ( j == 3 ) THEN
                                k = 1
                            ELSE
                                k = j + 1
                            END IF
                            eg2(1,:) = tr2_(j,:)
                            eg2(2,:) = tr2_(k,:)
                            
                            ! �Ƚ�eg1 eg2
                            sub1 = eg1(1,:) - eg2(1,:)
                            sub2 = eg1(2,:) - eg2(2,:)
                            sub3 = eg1(1,:) - eg2(2,:)
                            sub4 = eg1(2,:) - eg2(1,:)
                            ! ����ж����������ı����غϣ���ô����Ϊ������
                            IF ( (ALL( DABS(sub1) < 1.D-12) .AND. ALL( DABS(sub2) < 1.D-12)) .OR. &
                                 (ALL( DABS(sub3) < 1.D-12) .AND. ALL( DABS(sub4) < 1.D-12))  ) THEN
                                res_(1,:) =  eg1(1,:)
                                res_(2,:) =  eg1(2,:)
                                RETURN
                            END IF
                        END DO
                    END DO
                    
                    res_ = 0.D0
                    
                    ! û���ҵ�������
                    info_ = 1

                    RETURN
                END SUBROUTINE get_common_TriaSide
                
            
        END SUBROUTINE find_criticalEdge_to_H
        
        !---------------------------------------------------------------
        !
        ! �Ѽ��� V ��ÿ������ⲿ�㼯���ܵ�һ���㼯 L ��
        !
        !---------------------------------------------------------------
        SUBROUTINE together_OS_to_L( V_, P_, L_ )
            IMPLICIT NONE
            TYPE(Face_ptr), INTENT(IN) :: V_(:) !�ɼ��漯
            TYPE(Point), INTENT(IN) :: P_(:) ! s
            TYPE(Point_ptr), ALLOCATABLE, INTENT(INOUT) :: L_(:) ! V�е��ⲿ�㼯
            
            INTEGER*4 :: i, j, C, istat
            TYPE(Point_ptr), ALLOCATABLE :: L_temp(:) ! V�е��ⲿ�㼯temp
            ! ��ʼ��L
            IF ( ALLOCATED(L_) ) THEN
                FORALL (i = 1 : SIZE(L_) ) L_(i)%p => NULL() ! ��ָ��
                DEALLOCATE(L_) ! ��ɾ��
            END IF
            
            ! �Ѽ��� V ��ÿ������ⲿ�㼯���ܵ�һ���㼯 L ��
            ALLOCATE( L_( SIZE(P_) ), STAT = istat )    
                
            C = 0
            DO i = 1, SIZE(V_), 1
            DO j = 1, SIZE( V_(i)%p%OS ), 1
                IF ( ASSOCIATED( V_(i)%p%OS(j)%p ) ) THEN
                    C = C + 1
                    L_(C) = V_(i)%p%OS(j)
                ELSE
                    EXIT
                END IF
            END DO
            END DO
            
            ! ȥ������Ҫ�Ŀռ�
            C = 0
            DO i = 1, SIZE(L_), 1
                IF ( ASSOCIATED( L_(i)%p ) ) C = C + 1
            END DO
            
            ALLOCATE( L_temp(C), STAT = istat )

            L_temp(:) = L_(1:C)
            
            IF ( ALLOCATED(L_) ) THEN
                FORALL (i = 1 : SIZE(L_) ) L_(i)%p => NULL() ! ��ָ��
                DEALLOCATE(L_) ! ��ɾ��
            END IF
                
            L_ = L_temp
            
            RETURN
        END SUBROUTINE together_OS_to_L
        
        !---------------------------------------------------------------
        !
        ! Ѱ��p_far��Ե���ɵĿɼ��漯V
        !
        !---------------------------------------------------------------
        SUBROUTINE find_finished_visible_faceSet(p_, V_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(3)
            TYPE(Face_ptr), ALLOCATABLE, INTENT(INOUT) :: V_(:) !�ɼ��漯
            INTEGER*4 :: i, j ! 
            
            ! ������V����ӿɼ��棬һ���µ���ͱ�����������, ֱ�������涼true
            i = 1
            DO WHILE( i <= SIZE(V_) )
                ! �����ǰ�汻���Ϊtrue����ô������һ��ѭ��
                IF ( V_(i)%p%flag == .TRUE. ) THEN
                    i = i + 1
                    CYCLE
                END IF
                
                ! �����ʹ��Ŀɼ�������ҲҪ��ǣ����ұ����ȱ��
                V_(i)%p%flag = .TRUE.
                
                ! ������ÿһ��������
                DO j = 1, 3, 1
                    IF ( V_(i)%p%N(j)%p%flag == .FALSE. ) THEN
                        ! �����ʹ�������Ϊ.TRUE.
                        !V(i)%p%N(j)%p%flag = .TRUE.
                        
                        ! ���p�ڸ��������Ϸ�����ô�Ѹ���Ҳ����V
                        IF( DOT_PRODUCT( V_(i)%p%N(j)%p%nml(:), p_ - V_(i)%p%N(j)%p%vertex(1,:) ) > 1.D-12 ) THEN
                            CALL Store_Face_Into_Visible_FaceSet( V_(i)%p%N(j), V_ )
                        END IF
                    END IF
                END DO

                ! ������
                i = i + 1
            END DO
            
            RETURN
        END SUBROUTINE find_finished_visible_faceSet
        
        
        !---------------------------------------------------------------
        !
        ! ���´����convex��ÿ����ĳ���ķ�ʸ��������ṹ���ڵ�nml����
        !
        !---------------------------------------------------------------
        SUBROUTINE Calc_NormalVector_DerivForm_MassPoint(convex_)
            IMPLICIT NONE
            TYPE(Face), INTENT(INOUT) :: convex_(:)! 
            REAL*8 :: mp(3), fc(3), nml( SIZE( convex_ ), 3), AB(3), BC(3)
            INTEGER*4 :: i, j
            
            mp = 0.D0
            DO j = 1, SIZE( convex_ ), 1
                FORALL(i = 1:3) fc(i) = SUM( convex_(j)%vertex(:,i) ) / 3.D0      
                mp(:) = mp(:) + fc(:) !�������ռӳ����������������
            END DO
            mp = mp / REAL( SIZE( convex_ ), KIND = 8) ! ��Ҫ����������
            
            ! ����ÿ��������ĳ���ķ�ʸ
            DO i = 1, SIZE( convex_ ), 1
                AB(:) = convex_(i)%vertex(2,:) - convex_(i)%vertex(1,:)
                BC(:) = convex_(i)%vertex(3,:) - convex_(i)%vertex(2,:)
                nml(i,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
                IF ( DOT_PRODUCT( nml(i,:), convex_(i)%vertex(1,:) - mp  ) < .0D0 ) nml(i,:) = - nml(i,:)
            END DO
            
            FORALL(i = 1 : SIZE( convex_ )) convex_(i)%nml(:) = nml(i,:)
            
            RETURN
        END SUBROUTINE Calc_NormalVector_DerivForm_MassPoint
        
        !---------------------------------------------------------------
        !
        ! ���漯������flag ��Ϊ .false.
        !
        !---------------------------------------------------------------
        SUBROUTINE Initial_Flag_FaceSet(faceSet_)
            IMPLICIT NONE
            TYPE(Face), INTENT(INOUT) :: faceSet_(:) ! 
            INTEGER*4 :: i ! 
            DO i = 1, SIZE( faceSet_ ), 1
                faceSet_(i)%flag = .FALSE.
            END DO
            RETURN
        END SUBROUTINE Initial_Flag_FaceSet
        
        !---------------------------------------------------------------
        !
        ! ����F��ӽ��ɼ��漯V��
        !
        !---------------------------------------------------------------
        SUBROUTINE Store_Face_Into_Visible_FaceSet(face_ptr_, V_)
            IMPLICIT NONE
            TYPE(Face_ptr), INTENT(IN), TARGET :: face_ptr_ ! 
            TYPE(Face_ptr), ALLOCATABLE, INTENT(INOUT) :: V_(:) !�ɼ��漯
            INTEGER*4 :: i, C, istat, j, k, m
            TYPE(Face_ptr), ALLOCATABLE :: faceSet(:) ! �ݴ��漯
            
            IF ( .NOT. ALLOCATED(V_) ) THEN
                ALLOCATE( V_(1), STAT = istat )
                V_(1)%p => face_ptr_%p
            ELSE
                ! ��ʼ��faceSet
                IF ( ALLOCATED(faceSet) ) THEN
                    DO i = 1, SIZE(faceSet), 1
                        faceSet(i)%p => NULL()
                    END DO
                    DEALLOCATE( faceSet ) ! ���ٶ���
                END IF
                    
                ! �洢��ǰV
                faceSet = V_

                ! ����V��������V�е�ָ��
                IF ( ALLOCATED(V_) ) THEN
                    DO i = 1, SIZE( V_ ), 1
                        V_(i)%p => NULL()
                    END DO
                    DEALLOCATE( V_ ) ! ���ٶ���
                END IF
                    
                ! ����V�Ŀռ�
                ALLOCATE( V_( 1 + SIZE(faceSet) ), STAT = istat )
                    
                ! �õ���V
                V_(1 : SIZE(faceSet)) = faceSet
                V_(1 + SIZE(faceSet)) = face_ptr_
            END IF
            
            RETURN
        END SUBROUTINE Store_Face_Into_Visible_FaceSet
        
        !---------------------------------------------------------------
        !
        ! ��ʼ���ɼ��漯
        !
        !---------------------------------------------------------------
        SUBROUTINE Initial_Visible_FaceSet(V_)
            IMPLICIT NONE
            TYPE(Face_ptr), ALLOCATABLE, INTENT(INOUT) :: V_(:) !�ɼ��漯
            INTEGER*4 ::  j! 
            ! ��V��
            IF ( ALLOCATED(V_) ) THEN
                DO j = 1, SIZE( V_ ), 1
                    V_(j)%p => NULL()
                END DO
                DEALLOCATE( V_ ) ! ���ٶ���
            END IF
            RETURN
        END SUBROUTINE Initial_Visible_FaceSet
        
        !---------------------------------------------------------------
        !
        ! ��F���ⲿ�㼯���ҵ�������F��Զ�ĵ�p
        !
        !---------------------------------------------------------------
        SUBROUTINE Search_Farthest_Point_from_FaceSet(faceSet_ptr_, p_far_)
            IMPLICIT NONE
            TYPE(Face_ptr), INTENT(IN) :: faceSet_ptr_ ! 
            REAL*8, INTENT(OUT) :: p_far_(3)
            INTEGER*4 :: i, j, index_max
            REAL*8 :: dist, max
            
            dist = 0.D0
            max = 0.D0
            DO i = 1, SIZE( faceSet_ptr_%p%OS ), 1
                IF ( ASSOCIATED(faceSet_ptr_%p%OS(i)%p) ) THEN
                    dist = DABS( DIST_PF_SIGN( faceSet_ptr_%p%OS(i)%p%xyz(:), faceSet_ptr_%p%vertex(:,:) ) )
                    IF( dist > max ) THEN
                        max = dist
                        index_max = i
                    END IF
                END IF
            END DO
            
            p_far_ = faceSet_ptr_%p%OS( index_max )%p%xyz(:)
            
            RETURN
        END SUBROUTINE Search_Farthest_Point_from_FaceSet
        
        
        
        !---------------------------------------------------------------
        !
        ! ˢ���漯��ÿ�����������
        !
        !---------------------------------------------------------------
        SUBROUTINE Update_Adjacent_Faces(hull_)
            IMPLICIT NONE
            TYPE(Face), INTENT(INOUT), TARGET :: hull_(:)
            INTEGER*4 :: i, j, k, index, index_hull, C
            REAL*8 :: ep1(3), ep2(3)
            
            ! ��ʼ��͹��������
            DO i = 1, SIZE(hull_), 1
                DO j = 1, 3, 1
                    hull_(i)%N(j)%p => NULL()
                END DO
            END DO
            
            ! ˢ��͹����������
            ! ����ÿ����F
            DO i = 1, SIZE(hull_), 1
                ! ����ÿ�����3���� 1-2   2-3   3-1
                ! ��1-2 ��������
                ep1 = hull_(i)%vertex(1,:)
                ep2 = hull_(i)%vertex(2,:)
                index = get_index_adjacent_face(ep1, ep2, hull_, i)
                hull_(i)%N(1)%p => hull_(index)
                
                ! ��2-3 ��������
                ep1 = hull_(i)%vertex(2,:)
                ep2 = hull_(i)%vertex(3,:)
                index = get_index_adjacent_face(ep1, ep2, hull_, i)
                hull_(i)%N(2)%p => hull_(index)
                
                ! ��3-1 ��������
                ep1 = hull_(i)%vertex(3,:)
                ep2 = hull_(i)%vertex(1,:)
                index = get_index_adjacent_face(ep1, ep2, hull_, i)
                hull_(i)%N(3)%p => hull_(index)
            END DO
            
            RETURN
            
            CONTAINS
                PURE FUNCTION get_index_adjacent_face(ep1_, ep2_, hull_, index_faceNow_) RESULT(res_)
                IMPLICIT NONE
                    REAL*8, INTENT(IN) :: ep1_(3), ep2_(3)
                    TYPE(Face), INTENT(IN) :: hull_(:)
                    INTEGER*4, INTENT(IN) :: index_faceNow_ 
                    INTEGER*4 :: res_
                    INTEGER*4 :: i, j, k
                    REAL*8 :: sub1(3), sub2(3), sub3(3), sub4(3)
                    res_ = 0
                    ! ����hull_
                    DO i = 1, SIZE(hull_), 1
                        ! ��ǰ�治����
                        IF ( i == index_faceNow_ ) CYCLE
                        
                        ! �����ep1 2 �Ƿ���hull_(i)��ĳ�����غ�
                        DO j = 1, 3, 1
                            IF( j == 3) THEN
                                k = 1
                            ELSE
                                k = j + 1
                            END IF
                            sub1 = hull_(i)%vertex(j,:) - ep1_
                            sub2 = hull_(i)%vertex(k,:) - ep2_
                            sub3 = hull_(i)%vertex(j,:) - ep2_
                            sub4 = hull_(i)%vertex(k,:) - ep1_
                            ! ����ж����������ı����غϣ���ô����Ϊ������
                            IF ( (ALL( DABS(sub1) < 1.D-8) .AND. ALL( DABS(sub2) < 1.D-8)) .OR. &
                                 (ALL( DABS(sub3) < 1.D-8) .AND. ALL( DABS(sub4) < 1.D-8))  ) THEN
                                res_ = i
                                RETURN
                            END IF
                        END DO
                    END DO
                    
                    RETURN
                END FUNCTION get_index_adjacent_face
                
            
        END SUBROUTINE Update_Adjacent_Faces
        
        !---------------------------------------------------------------
        !
        ! �ѵ�ǰHull�ⲿ�㼯�ǿյ��汣��������漯Q�У�
        !
        !---------------------------------------------------------------
        SUBROUTINE Store_Nonempty_Faces(Hull_, Q_)
            IMPLICIT NONE
            TYPE(Face), INTENT(IN), TARGET :: Hull_(:) !͹��
            TYPE(Face_ptr), INTENT(INOUT), ALLOCATABLE :: Q_(:) !�����漯
            INTEGER*4 :: i, C, istat, j, k, m
            TYPE(Face_ptr), ALLOCATABLE :: Q_last(:) ! �ݴ�����漯
            
            ! ����ʼ��Q��
            IF ( ALLOCATED(Q_) ) THEN
                DO j = 1, SIZE( Q_ ), 1
                    Q_(j)%p => NULL()
                END DO
                DEALLOCATE( Q_ ) ! ���ٶ���
            END IF
            
            C = 0
            ! Ѱ�ҷǿ��漯��
            DO i = 1, SIZE(Hull_), 1
            ! ����ÿ��Hull��ÿ������ⲿ�㼯
            DO m = 1, SIZE( Hull_(i)%OS ), 1
                ! ������Щ�ⲿ�㼯���Ƿ��з���
                IF( ASSOCIATED( Hull_(i)%OS(m)%p ) ) THEN
                    C = C + 1
                    
                    IF ( C == 1 ) THEN
                        ALLOCATE( Q_(1), STAT = istat )
                        Q_(1)%p => Hull_(i)
                    ELSE
                        ! ����Q_last��������Q_last�е�ָ��
                        IF ( ALLOCATED(Q_last) ) THEN
                            DO j = 1, SIZE( Q_last ), 1
                                Q_last(j)%p => NULL() ! ��ָ��
                            END DO
                            DEALLOCATE( Q_last ) ! ���ٶ���
                        END IF
                    
                        ! �洢��ǰQ
                        Q_last = Q_
                    
                        ! ����Q��������Q�е�ָ��
                        IF ( ALLOCATED(Q_) ) THEN
                            DO j = 1, SIZE( Q_ ), 1
                                Q_(j)%p => NULL() ! ��ָ��
                            END DO
                            DEALLOCATE( Q_ ) ! ���ٶ���
                        END IF
                    
                        ! ����Q�Ŀռ�
                        ALLOCATE( Q_( C ), STAT = istat )
                    
                        ! �õ���Q
                        Q_(1:C-1) = Q_last
                        Q_(C)%p => Hull_(i)
                    END IF
                        
                    ! ����������һ����F
                    EXIT
                        
                END IF
            END DO
            END DO
            
            
            RETURN
        END SUBROUTINE Store_Nonempty_Faces
        
        
        !---------------------------------------------------------------
        !
        ! ���ⲿ�㼯�ǿյ��汣��������漯Q�У�
        !
        !---------------------------------------------------------------
        SUBROUTINE initial_p(P_)
            IMPLICIT NONE
            TYPE(Point), ALLOCATABLE, INTENT(INOUT) :: P_(:) ! s
            ! ��P��
            IF ( ALLOCATED(P_) )  DEALLOCATE( P_ ) ! ���ٶ���
            RETURN
        END SUBROUTINE initial_p
        

        
        !---------------------------------------------------------------
        !
        ! ��p���䵽��ʼ͹����ÿ���� F ���ⲿ�㼯OS��
        !
        !---------------------------------------------------------------
        SUBROUTINE Distr_OS_for_Face(P_, Hull_)
            IMPLICIT NONE
            TYPE(Point), INTENT(INOUT), TARGET :: P_(:) ! s
            TYPE(Face), INTENT(INOUT) :: Hull_(:) !͹��
            REAL*8 :: mp(3), fc(3), nml( SIZE( Hull_ ), 3)
            REAL*8 :: AB(3), BC(3), debug, debug_vec(3)
            INTEGER*4 :: i, j, istat, C
            mp = 0.D0
            DO j = 1, SIZE( Hull_ ), 1
                FORALL(i = 1:3) fc(i) = SUM( Hull_(j)%vertex(:,i) ) / 3.D0      
                mp(:) = mp(:) + fc(:) !�������ռӳ����������������
            END DO
            mp = mp / REAL( SIZE( Hull_ ), KIND = 8) ! ��Ҫ����������
            
            ! ����ÿ��������ĳ���ķ�ʸ
            DO i = 1, SIZE( Hull_ ), 1
                AB(:) = Hull_(i)%vertex(2,:) - Hull_(i)%vertex(1,:)
                BC(:) = Hull_(i)%vertex(3,:) - Hull_(i)%vertex(2,:)
                nml(i,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
                IF ( DOT_PRODUCT( nml(i,:), Hull_(i)%vertex(1,:) - mp  ) < .0D0 ) nml(i,:) = - nml(i,:)
            END DO
             
            ! ��ʼ��F���ⲿ�㼯OS,ʹÿ��F���ⲿ�㼯�ռ����ܵ㼯P�Ĵ�С��ͬ
            DO i = 1, SIZE( Hull_ ), 1
                ALLOCATE( Hull_(i)%OS( SIZE(P_) ), STAT = istat )
                DO j = 1, SIZE(P_), 1
                    Hull_(i)%OS(j)%p => NULL()
                END DO
            END DO
            
            ! �ж�p�Ƿ�����F�⣬�����жϳ�����ĳ����F���ⲿ�㼯��������ѭ�������ظ�����OS
            DO i = 1, SIZE(P_), 1
            ! ����Hull�е���F
            DO j = 1, SIZE( Hull_ ), 1
                    !debug = DOT_PRODUCT( P_(i)%xyz(:) - Hull_(j)%vertex(1,:), nml(j,:) )
                IF( DOT_PRODUCT( P_(i)%xyz(:) - Hull_(j)%vertex(1,:), nml(j,:) ) > 1.D-8) THEN ! ���YES��˵������F��
                    !debug_vec = P_(i)%xyz(:)
                    !debug_vec = Hull_(j)%vertex(1,:)
                    !debug_vec = nml(j,:)
                    C = 0
                    DO WHILE(.TRUE.) ! ֱ����ָ���λ��
                        C = C + 1
                        IF (.NOT. ASSOCIATED( Hull_(j)%OS(C)%p )) THEN
                            Hull_(j)%OS(C)%p => P_(i)
                            P_(i)%isOn = .TRUE.
                            EXIT ! ������������ⲿ�������
                        END IF
                    END DO
                    EXIT ! �õ������ϣ�������һ����p
                END IF
            END DO
            END DO
            
            RETURN
        END SUBROUTINE Distr_OS_for_Face
        
        
        !---------------------------------------------------------------
        !
        ! �����ṹ��
        !
        !---------------------------------------------------------------
        SUBROUTINE Point_into_Type(p_, P_type_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(:,:)
            TYPE(Point), INTENT(INOUT) :: P_type_(:) ! 
            INTEGER*4 :: i ! 
            FORALL (i = 1 : SIZE(p_, 1)) P_type_(i)%xyz(:) =  p_(i,:)
            FORALL (i = 1 : SIZE(p_, 1)) P_type_(i)%isOn = .FALSE.
            RETURN
        END SUBROUTINE Point_into_Type
        
        !---------------------------------------------------------------
        !
        ! ��ȡ��ʼ4��������
        !
        !---------------------------------------------------------------
        SUBROUTINE Get_Start_Convex(p_, convex_, info_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(:,:)
            TYPE(Face), INTENT(OUT) :: convex_(4)
            INTEGER*4, INTENT(OUT) :: info_ ! 
            
            REAL*8 :: vertex(4,3) = .0D0
            REAL*8 :: convex_vertex(4,3)
            TYPE(Face) :: convex_temp(4)
            INTEGER*4 :: max_loc_coord(3), min_loc_coord(3)
            INTEGER*4 :: i, j, istat, C
            REAL*8 :: dir(3), r, dot
            
            convex_vertex = .0D0
            
            ! ������Ϣ
            info_ = 0
            
            !---------------------------------------------------------------
            ! ��һ��������ʼ͹����  
            !---------------------------------------------------------------
            !---------------------------------------------------------------
            ! 0��{�������ֱ�Ӽ��ж�}  �������ĵ����4���㣬����4���㲻���棬��ôֱ�ӷ������4���������ɵĵ�������Ϊ͹�������򱨴�
            IF ( SIZE(p_, 1) == 4 ) THEN 
                IF( .NOT. is4PointsOnCommonPlane(p_) ) THEN ! 4�㲻����
                    ! face 1 
                    convex_(1)%vertex(1,:) =  p_(1,:)
                    convex_(1)%vertex(2,:) =  p_(2,:)
                    convex_(1)%vertex(3,:) =  p_(3,:)
                    ! face 2                      
                    convex_(2)%vertex(1,:) =  p_(1,:)
                    convex_(2)%vertex(2,:) =  p_(3,:)
                    convex_(2)%vertex(3,:) =  p_(4,:)
                    ! face 3                    
                    convex_(3)%vertex(1,:) =  p_(1,:)
                    convex_(3)%vertex(2,:) =  p_(4,:)
                    convex_(3)%vertex(3,:) =  p_(2,:)
                    ! face 4                    
                    convex_(4)%vertex(1,:) =  p_(2,:)
                    convex_(4)%vertex(2,:) =  p_(4,:)
                    convex_(4)%vertex(3,:) =  p_(3,:)
                    
                    RETURN
                ELSE 
                    info_ = 2 ! �����4���㹲��
                    RETURN
                END IF
            END IF
            
            !---------------------------------------------------------------
            ! 1��Ѱ��ǰ������
            C = 0
            DO WHILE(.TRUE.)
                C = C + 1
                IF ( C > 99 ) THEN
                    info_ = 1
                    RETURN
                    !WRITE(UNIT=6, FMT="(A)") "Cannot get start convex. - SUBROUTINE Get_Start_Convex(p_, convex_)" !
                    !READ(UNIT=5, FMT=*)   
                    !STOP
                END IF
                !---------------------------------------------------------------
                ! ��1.1 ��ʼ��Ѱ���� 1  �������
                dir = GET_RANDOM_UNIT_VECTOR(C)
                
                ! �����ʼ��Զ�� 1����������ӵ�convex_vertex��
                convex_vertex(1,:) = support_point(p_, dir)
                !---------------------------------------------------------------
                ! ��1.2 ��ʼ��Ѱ���� 2 �� ����
                dir = - dir
            
                ! �����ʼ��Զ�� 2����������ӵ�convex_vertex��
                convex_vertex(2,:) = support_point(p_, dir)
            
                ! �����ʼ������Զ�㶼�غϣ���ô�ı���Ѱ����ֱ�����غ�Ϊֹ
                IF ( ALL( DABS( convex_vertex(1,:) - convex_vertex(2,:) ) < 1.D-8 ) ) THEN
                    CYCLE
                ELSE
                    EXIT
                END IF
            END DO
            
            !---------------------------------------------------------------
            ! 2��Ѱ�ҵ�3����
            C = 0
            DO WHILE(.TRUE.)
                C = C + 1
                IF ( C > 99 ) THEN
                    info_ = 1
                    RETURN
                    !WRITE(UNIT=6, FMT="(A)") "Cannot get start convex. - SUBROUTINE Get_Start_Convex(p_, convex_)" !
                    !READ(UNIT=5, FMT=*)   
                    !STOP
                END IF
                !---------------------------------------------------------------
                ! ��1.3 ��Ѱ���� 3  ���-1 - 1��
                dir = GET_RANDOM_UNIT_VECTOR(C)
                
                ! ������Զ��3����������ӵ�convex_vertex��
                convex_vertex(3,:) = support_point(p_, dir)
                
                ! ������������ǰ�����㹲�ߣ���ô������Ѱ����
                IF ( NORM2( CROSS_PRODUCT_3D(convex_vertex(3,:) - convex_vertex(2,:), convex_vertex(3,:) - convex_vertex(1,:)) ) < 1.D-12 ) THEN
                    CYCLE
                ELSE
                    EXIT
                END IF
            END DO
            
            !---------------------------------------------------------------
            ! 2��Ѱ�ҵ�4����
            C = 0
            DO WHILE(.TRUE.)
                C = C + 1
                ! ������������10 ʱ����
                IF ( C > 99 ) THEN
                    info_ = 1
                    RETURN
                    !WRITE(UNIT=6, FMT="(A)") "Cannot get start convex. - SUBROUTINE Get_Start_Convex(p_, convex_)" !
                    !READ(UNIT=5, FMT=*)   
                    !STOP
                END IF
                !---------------------------------------------------------------
                ! ��1.3 ��Ѱ���� 4  ���-1 - 1��
                dir = GET_RANDOM_UNIT_VECTOR(C)
                
                ! ������Զ�� 4����������ӵ�convex_vertex��
                convex_vertex(4,:) = support_point(p_, dir)
                
                ! �����4�����ǰ4���㹲�棬��ô������Ѱ���򣬵�����������������ʱ����
                dot = DIST_PF_SIGN(convex_vertex(4,:), convex_vertex(1:3,:))
                IF ( DABS(dot)  < 1.D-8 ) THEN ! 20230706
                    ! �ⲿ��֪�㲻���棬�������������Ȼ����һ����������ô����ǿ��ָ����ʼ͹���ĵ��ĸ��㣬�õ���ǰ3������ɵ��治����
                    IF ( C >= 10 ) THEN
                        CALL forciblySpecifyPoint4(p_, convex_vertex(1:3,:), convex_vertex(4,:), info_)
                        ! �����Ȼ�޷��ҵ�͹����˵�������߼��������⣬����ִ��
                        IF ( info_ /= 0 ) THEN
                            info_ = 1
                            RETURN
                        ELSE
                            info_ = 0
                            EXIT
                        END IF
                    END IF
                    CYCLE
                ELSE
                    EXIT
                END IF
            END DO
            
            ! face 1 
            convex_temp(1)%vertex(1,:) =  convex_vertex(1,:)
            convex_temp(1)%vertex(2,:) =  convex_vertex(2,:)
            convex_temp(1)%vertex(3,:) =  convex_vertex(3,:)
            ! face 2 
            convex_temp(2)%vertex(1,:) =  convex_vertex(1,:)
            convex_temp(2)%vertex(2,:) =  convex_vertex(3,:)
            convex_temp(2)%vertex(3,:) =  convex_vertex(4,:)
            ! face 3 
            convex_temp(3)%vertex(1,:) =  convex_vertex(1,:)
            convex_temp(3)%vertex(2,:) =  convex_vertex(4,:)
            convex_temp(3)%vertex(3,:) =  convex_vertex(2,:)
            ! face 4 
            convex_temp(4)%vertex(1,:) =  convex_vertex(2,:)
            convex_temp(4)%vertex(2,:) =  convex_vertex(4,:)
            convex_temp(4)%vertex(3,:) =  convex_vertex(3,:)
            
            convex_ = convex_temp
            
            RETURN
            
        CONTAINS
                !����ǿ��ָ����ʼ͹���ĵ��ĸ��� info = 0 ������ info = 1 �����߼�����
                SUBROUTINE forciblySpecifyPoint4(pAll_, p1_p2_p3_, p4_, info_)
                    IMPLICIT NONE
                    REAL*8, INTENT(IN) :: pAll_(:,:), p1_p2_p3_(3,3)
                    REAL*8, INTENT(OUT) :: p4_(3)
                    INTEGER*4 :: info_ !
                    
                    INTEGER*4 :: i ! 
                    REAL*8 :: cross(3), AB(3), BC(3)
                    
                    AB(:) =  p1_p2_p3_(2,:) - p1_p2_p3_(1,:)
                    BC(:) =  p1_p2_p3_(3,:) - p1_p2_p3_(2,:)
                    cross(:) = CROSS_PRODUCT_3D(AB, BC)
                    
                    IF ( ALL(DABS(cross) < 1.D-12) ) THEN
                        info_ = 1 ! �����߼�����
                        RETURN
                    END IF
                    
                    ! ֱ����һ��������ĵ�
                    DO i = 1, SIZE(pAll_, 1), 1
                        ! �����봫������p123�غϵĵ������
                        IF ( ALL(DABS(pAll_(i, :) - p1_p2_p3_(1,:)) < 1.D-12) .OR. &
                             ALL(DABS(pAll_(i, :) - p1_p2_p3_(2,:)) < 1.D-12) .OR. &
                             ALL(DABS(pAll_(i, :) - p1_p2_p3_(3,:)) < 1.D-12) ) THEN
                            CYCLE
                        END IF
                             
                        ! �жϸõ��Ƿ���p123��ɵ��湲��
                        IF ( DABS( DIST_PF_SIGN(pAll_(i,:), p1_p2_p3_(1:3,:)) ) < 1.D-12 ) THEN
                            CYCLE
                        ELSE
                            p4_ = pAll_(i,:)
                        END IF 
                    END DO
            
                    RETURN
                END SUBROUTINE forciblySpecifyPoint4
        
                ! �ж��Ƿ��ĵ㹲��
                FUNCTION is4PointsOnCommonPlane(p_) RESULT(res_)
                    IMPLICIT NONE
                    REAL*8, INTENT(IN) :: p_(4,3)
                    LOGICAL*1 :: res_ ! 
                    REAL*8, ALLOCATABLE :: p_filted(:,:)
                    INTEGER*4 :: info
                    
                    CALL FILT_DUPL(p_, p_filted, info)
                    IF ( SIZE(p_filted, 1) /= 4 ) THEN
                        WRITE(UNIT=6, FMT="(A)") "ERROR - FUNCTION is4PointsOnCommonPlane(p_) RESULT(res_)" ! 
                        READ(UNIT=5, FMT=*)   
                        STOP
                    END IF
                    
                    IF ( DABS( DIST_PF_SIGN(p_(1,:), p_(2:4,:)) ) < 1.D-12 ) THEN
                        res_ = .TRUE.
                    ELSE
                        res_ = .FALSE.
                    END IF
                    RETURN
                END FUNCTION is4PointsOnCommonPlane
                
            
                ! Ѱ�ҵ㼯�и÷�������Զ�ĵ�
                FUNCTION support_point(p_, dir_) RESULT(res_)
                    IMPLICIT NONE
                    REAL(8), INTENT(IN) :: p_(:,:), dir_(3)  ! ������������㼯����������
                    REAL(8) :: res_(3)                       ! ������������õ�����Զ��
                    INTEGER*4 :: i, max_index! 
                    REAL*8 :: max_dot_product, dot_product_temp
                    
                    ! ����͹�����p1��dir_�����ϵ���Զ��
                    max_dot_product = - HUGE(1.0D0)
                    max_index = 1
                    DO i = 1, SIZE(p_, 1)
                        dot_product_temp = DOT_PRODUCT(dir_, p_(i,:))
                        IF (dot_product_temp > max_dot_product) THEN
                            max_dot_product = dot_product_temp
                            max_index = i
                        END IF
                    END DO
                    res_ = p_(max_index,:)
                    RETURN
                END FUNCTION support_point
        END SUBROUTINE Get_Start_Convex


        !_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
        !_/                                                                    _/
        !_/                   Math tools for GCLIB_QuickHull                         _/
        !_/                                                                    _/
        !_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
        !---------------------------------------------------------------
        !
        ! 3d�������  
        PURE FUNCTION CROSS_PRODUCT_3D(A_, B_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: A_(3), B_(3)
            REAL*8 :: res_(3) ! return value
            res_ = [A_(2) * B_(3) - A_(3) * B_(2), &
                    A_(3) * B_(1) - A_(1) * B_(3), &
                    A_(1) * B_(2) - A_(2) * B_(1)]
            RETURN
        END FUNCTION CROSS_PRODUCT_3D
        
        !---------------------------------------------------------------
        !
        ! �жϵ��Ƿ��ڵ������ڲ�
        FUNCTION isPointInSimplex(p_, simplex_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(3)
            REAL*8, INTENT(IN) :: simplex_(4,3)
            LOGICAL*1 :: res_ ! �������ʾԭ���Ƿ��ڵ������ڲ����߼�����
            LOGICAL*1 :: is ! 
            REAL*8 :: M(3)                     ! �����ε�����ƽ��
            REAL*8 :: MO(3)                     ! �����ε�����ƽ�� -> ԭ��
            REAL*8 :: AB(3), BC(3), dist(4), nml(4,3) 
            REAL*8 :: tri(3,3)
            REAL*8 :: vertexOnPlane(3,3)
            INTEGER*4 :: i ! 
            !---------------------------------------------------------------
            ! ��SITU 1���������ڵ����
            FORALL(i = 1:3) M(i) = SUM( simplex_(:,i) ) / 4.D0      
            
            ! ����������4���泯��ĵ�λ��ʸ
            ! face 1
            AB = simplex_(1,:) - simplex_(3,:)
            BC = simplex_(3,:) - simplex_(4,:)
            nml(1,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
            IF ( DOT_PRODUCT( nml(1,:), simplex_(1,:) - M  ) < .0D0 ) nml(1,:) = - nml(1,:)
            
            ! face 2
            AB = simplex_(1,:) - simplex_(2,:)
            BC = simplex_(2,:) - simplex_(4,:)
            nml(2,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
            IF ( DOT_PRODUCT( nml(2,:), simplex_(1,:) - M ) < .0D0 ) nml(2,:) = - nml(2,:)
            
            ! face 3
            AB = simplex_(1,:) - simplex_(2,:)
            BC = simplex_(2,:) - simplex_(3,:)
            nml(3,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
            IF ( DOT_PRODUCT( nml(3,:), simplex_(1,:) - M ) < .0D0 ) nml(3,:) = - nml(3,:)
            
            ! face 4
            AB = simplex_(2,:) - simplex_(3,:)
            BC = simplex_(3,:) - simplex_(4,:)
            nml(4,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
            IF ( DOT_PRODUCT( nml(4,:), simplex_(2,:) - M ) < .0D0 ) nml(4,:) = - nml(4,:)

            DO i = 1, 4, 1
                dist(i) = DOT_PRODUCT(simplex_(i,:) - p_, nml(i,:))        
            END DO
            
            ! ��SITU 2������������ϵ����
            DO i = 1, 4, 1
                IF ( isZero8( dist(i) ) ) THEN
                    SELECT CASE( i )
                      CASE (1)
                        vertexOnPlane(1,:) = simplex_(1,:)
                        vertexOnPlane(2,:) = simplex_(3,:)
                        vertexOnPlane(3,:) = simplex_(4,:)
                      CASE (2)
                        vertexOnPlane(1,:) = simplex_(1,:)
                        vertexOnPlane(2,:) = simplex_(2,:)
                        vertexOnPlane(3,:) = simplex_(4,:)
                      CASE (3)
                        vertexOnPlane(1,:) = simplex_(1,:)
                        vertexOnPlane(2,:) = simplex_(2,:)
                        vertexOnPlane(3,:) = simplex_(3,:)
                      CASE (4)
                        vertexOnPlane(1,:) = simplex_(2,:)
                        vertexOnPlane(2,:) = simplex_(3,:)
                        vertexOnPlane(3,:) = simplex_(4,:)
                    END SELECT
                      
                    is = IS_INSIDE_PF(vertexOnPlane, p_)
                    
                    IF ( is ) THEN
                        res_ = .TRUE.
                        RETURN
                    END IF
                      
                END IF
            END DO
            
            IF ( dist(1) > .0D0 .AND. dist(2) > .0D0 .AND. dist(3) > .0D0 .AND. dist(4) > .0D0 ) THEN
                res_  = .TRUE.

            ELSE
                res_  = .FALSE.
            END IF
            
            RETURN
        END FUNCTION isPointInSimplex
        
        !---------------------------------------------------------------
        !
        ! �ж�ʵ���Ƿ�Ϊ0 (˫����)
        FUNCTION isZero8(RN) ! intent(inout)
            IMPLICIT NONE
            LOGICAL*1 :: isZero8 ! return value
            REAL*8, INTENT(IN) :: RN
            !---------------------------------------------------------------
            IF ( DABS(RN) < 1.D-8 ) THEN
                isZero8 = .TRUE.
            ELSE    
                isZero8 = .FALSE.
            END IF
        END FUNCTION isZero8
        
        
        !---------------------------------------------------------------
        !
        ! ����ƽ�������ϵ��������㷨���������ڶ���α��ϣ�
        PURE FUNCTION IS_INSIDE_PF(vertexOnPlane_, arbitraryPoint_) RESULT(res_)
            IMPLICIT NONE
            LOGICAL*1 :: res_ ! 
            REAL*8, INTENT(IN) :: vertexOnPlane_( :, :)
            REAL*8, INTENT(IN) :: arbitraryPoint_(3)
            !---------------------------------------------------------------
            REAL*8 :: temp
            INTEGER*4 :: i, NNODE
            REAL*8 :: crossProdResu( SIZE(vertexOnPlane_, 1) )    ! �洢��˽���������ţ�
            REAL*8 :: V( SIZE(vertexOnPlane_, 1) , 3)
            LOGICAL*1 :: zeroMask( SIZE(vertexOnPlane_, 1) )
            !---------------------------------------------------------------
            V = vertexOnPlane_
            zeroMask = .FALSE. 
            NNODE = SIZE(vertexOnPlane_, 1)
            !---------------------------------------------------------------
            ! ! �����  �ʵ�-�ǵ� X ������
            !����ͶӰ�� XOYƽ�棬���һ��ƽ������
            DO i = 1, NNODE, 1
                IF ( i == NNODE) THEN
                    crossProdResu(i) = (V(1,1) - V(i,1)) * (arbitraryPoint_(2) - V(i,2)) &
                                    - (V(1,2) - V(i,2)) * (arbitraryPoint_(1) - V(i,1))
                    EXIT
                END IF

                crossProdResu(i) = (V(i+1,1) - V(i,1)) * (arbitraryPoint_(2) - V(i,2)) &
                                   - (V(i+1,2) - V(i,2)) * (arbitraryPoint_(1) - V(i,1))
            END DO

            ! ����
            FORALL ( i = 1 : NNODE, DABS(crossProdResu(i)) < 1.0D-12 ) crossProdResu(i) = .0D0
            
            ! Ϊ�˱��ⱻ�жϵ������ʾ�ĵ�Ͷ��XOYƽ�����ֶ�㹲�߶�����crossProdResu=0�����
            ! ������������crossProdResu = 0 ������Ͷ��XOZƽ��������һ���ж�
            DO i = 1, NNODE, 1
                IF ( crossProdResu(i) > 1.0D-15 ) THEN ! ���ַ���
                    zeroMask(i) = .TRUE.
                END IF
            END DO
            IF ( ANY(zeroMask) == .FALSE. ) THEN !ȫΪ0
                !����ͶӰ�� XOZƽ�棬���һ��ƽ������
                DO i = 1, NNODE, 1
                    IF ( i == NNODE) THEN
                        crossProdResu(i) = (V(1,1) - V(i,1)) * (arbitraryPoint_(3) - V(i,3)) &
                                        - (V(1,3) - V(i,3)) * (arbitraryPoint_(1) - V(i,1))
                        EXIT
                    END IF

                    crossProdResu(i) = (V(i+1,1) - V(i,1)) * (arbitraryPoint_(3) - V(i,3)) &
                                       - (V(i+1,3) - V(i,3)) * (arbitraryPoint_(1) - V(i,1))
                END DO
            END IF
            
            
            ! ��������ÿһ��Ԫ�ض����һ��Ԫ����ˣ� �����ָ��ģ�˵�����ֲ�ͬ��
            DO i = 1, NNODE, 1
                temp = crossProdResu(1) * crossProdResu(i)
                ! �����ָ��ģ�˵�����ֲ�ͬ�ţ���ô�õ㲻��ƽ����
                IF ( temp < .0D0 ) THEN
                    res_ = .FALSE.
                    RETURN
                END IF
            END DO

            res_ = .TRUE.
            RETURN
        END FUNCTION IS_INSIDE_PF
        
        
        !---------------------------------------------------------------
        !
        ! ��ȡ��λ����
        !
        !---------------------------------------------------------------
        PURE FUNCTION UTZVEC(vtr_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: vtr_(:)
            REAL*8, ALLOCATABLE :: res_(:) ! return value
            REAL*8 :: MD
            ALLOCATE( res_( SIZE(vtr_) ) )
            MD = NORM2( vtr_ )
            res_ = MERGE(.0D0, vtr_ / MD, MD < 1.D-12)
            RETURN
        END FUNCTION UTZVEC
        
        !---------------------------------------------------------------
        !
        ! �㵽����ƽ��Ĵ�ֱ����(�����ţ�����ж�)
        FUNCTION DIST_PF_SIGN(arbitraryPoint_, defi3PoinOnPlane_) RESULT(res_)
            IMPLICIT NONE
            REAL*8 :: res_ ! return value
            REAL*8, INTENT(IN) :: arbitraryPoint_(3), defi3PoinOnPlane_(3,3)
            !---------------------------------------------------------------
            REAL*8 :: x_xyz(3) , p_xyz(3) 
            REAL*8 :: n_vtr(3) 
            x_xyz = arbitraryPoint_
            p_xyz = defi3PoinOnPlane_(1,:)
            n_vtr = UNINML(defi3PoinOnPlane_)
            
            ! �������ĵ㲻��ȷ��һ���棬��ô����
            IF ( ALL( DABS(n_vtr) < 1.D-12 ) ) THEN
                WRITE(UNIT=6, FMT="(A)") "ERROR - PURE FUNCTION DIST_PF_SIGN(arbitraryPoint_, defi3PoinOnPlane_) RESULT(res_)" ! 
                READ(UNIT=5, FMT=*)   
                STOP
            END IF
            
            res_ = DOT_PRODUCT(x_xyz - p_xyz, n_vtr)            
            RETURN
        END FUNCTION DIST_PF_SIGN
        
        !---------------------------------------------------------------
        !
        ! ����ƽ��ĵ�λ��ʸ
        PURE FUNCTION UNINML(plane3Vertex_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: plane3Vertex_(3,3)
            REAL*8 :: res_(3)
            !---------------------------------------------------------------
            REAL*8 :: N1(3), N2(3), cross(3)
            !---------------------------------------------------------------
            N1 = plane3Vertex_(2,:) - plane3Vertex_(1,:)
            N2 = plane3Vertex_(3,:) - plane3Vertex_(2,:)
            cross = CROSS_PRODUCT_3D(N1, N2)
            res_ = MERGE(cross / NORM2(cross), .0D0, ANY( DABS(cross) > 1.D-12 ))
            RETURN
        END FUNCTION UNINML
        
        !---------------------------------------------------------------
        !
        ! �жϸ�����һϵ�пռ���Ƿ�ȫ���غ� 
        PURE FUNCTION OVERLAP(a_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, DIMENSION(:, :), INTENT(IN) :: a_
            LOGICAL :: res_
            INTEGER :: i, j, n_points
            REAL*8, PARAMETER :: TOLERANCE = 1.D-12

            n_points = SIZE(a_, 1)
            res_ = .TRUE.

            DO i = 1, n_points - 1
            DO j = i + 1, n_points
                IF (ANY(ABS(a_(i, :) - a_(j, :)) > TOLERANCE)) THEN
                res_ = .FALSE.
                RETURN
                END IF
            END DO
            END DO
            RETURN
        END FUNCTION OVERLAP
        
        !---------------------------------------------------------------
        !
        ! ����ֱ��ĳһ��Ĵ���ָ��õ�ĵ�λ����
        PURE FUNCTION VEC_PL(arbitraryPoint_, defi2PoinOnLine_) RESULT(res_)
            IMPLICIT NONE
            REAL*8 :: res_(3) ! return value
            REAL*8, INTENT(IN) :: arbitraryPoint_(3), defi2PoinOnLine_(2,3)
            !---------------------------------------------------------------
            REAL*8 :: A(3), B(3), AB(3), C(3), D(3), AC(3)
            REAL*8 :: vec(3)    ! ����
            !---------------------------------------------------------------
            A = defi2PoinOnLine_(1,:)
            B = defi2PoinOnLine_(2,:)
            C = arbitraryPoint_
            AB = B - A
            AC = C - A
            D = A + DOT_PRODUCT(AC, AB) / NORM2(AB) * UTZVEC(AB)
            !---------------------------------------------------------------
            res_ = UTZVEC(D - C)
            RETURN
        END FUNCTION VEC_PL

        !---------------------------------------------------------------
        !
        ! ��ȡ�ڿռ��ϵõ���ֱ�߼���̾��룬��Ӧ�ֱ�������ֱ���ϵĵ㣨2�����㣩
        FUNCTION FOOT_LL(lineDefiEP_1_, lineDefiEP_2_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: lineDefiEP_1_(2,3), lineDefiEP_2_(2,3)
            REAL*8 :: res_(2,3)
            !---------------------------------------------------------------
            REAL*8 :: P1(3) , Q1(3) , P2(3) , Q2(3) 
            REAL*8 :: d1_vec(3) , d2_vec(3) , r_vec(3) 
            REAL*8 :: a , b , c , d , e , f 
            REAL*8 :: s , t 
            REAL*8 :: L1_s(3) , L2_t(3) 
            !---------------------------------------------------------------
            P1 = lineDefiEP_1_(1,:)
            Q1 = lineDefiEP_1_(2,:)
            P2 = lineDefiEP_2_(1,:)
            Q2 = lineDefiEP_2_(2,:)

            d1_vec = Q1 - P1
            d2_vec = Q2 - P2
            r_vec = P1 - P2

            a = DOT_PRODUCT(d1_vec, d1_vec)
            b = DOT_PRODUCT(d1_vec, d2_vec)
            c = DOT_PRODUCT(d1_vec, r_vec)
            e = DOT_PRODUCT(d2_vec, d2_vec)
            f = DOT_PRODUCT(d2_vec, r_vec)

            d = a * e - b**2

            IF ( DABS(d) < 1.0D-12 ) THEN !��ֱ��ƽ�У���ôȡ��һ���ߵ��е���Ϊ����
                res_(1,:) = (P1 + Q1) / 2.D0 
                res_(2,:) = FOOT_PL(res_(1,:), lineDefiEP_2_)
            ELSE !��ֱ�߲�ƽ��
                s = (b * f - c * e) / d
                t = (a * f - b * c) / d
                L1_s = P1 + s * (Q1 - P1)
                L2_t = P2 + t * (Q2 - P2)

                res_(1,:) = L1_s
                res_(2,:) = L2_t
            END IF
            RETURN
        END FUNCTION FOOT_LL
        
        !---------------------------------------------------------------
        !
        ! �㵽����ֱ�ߵĴ���
        PURE FUNCTION FOOT_PL(arbitraryPoint_, defi2PoinOnLine_) RESULT(res_)
            IMPLICIT NONE
            REAL*8 :: res_(3) ! return value
            REAL*8, INTENT(IN) :: arbitraryPoint_(3), defi2PoinOnLine_(2,3)
            !---------------------------------------------------------------
            REAL*8 :: P(3)   ! ���жϵĵ�
            REAL*8 :: V(2,3)  ! ���жϵ�ֱ���ϵ���������
            !---------------------------------------------------------------
            P = arbitraryPoint_
            V = defi2PoinOnLine_
            !---------------------------------------------------------------
            res_ = V(1,:) + DOT_PRODUCT( P - V(1,:), UTZVEC(V(2,:) - V(1,:)) ) * UTZVEC(V(2,:) - V(1,:))
            RETURN
        END FUNCTION FOOT_PL
        
        FUNCTION SET_DIFF(A_, B_) RESULT(res_)
            INTEGER*4, INTENT(IN) :: A_(:), B_(:)
            INTEGER*4, ALLOCATABLE :: res_(:)
            INTEGER*4 :: i, j, count
            LOGICAL :: inB

            ! ����res_�Ĵ�С
            count = 0
            DO i = 1, SIZE(A_)
                inB = .FALSE.
                DO j = 1, SIZE(B_)
                    IF (A_(i) == B_(j)) THEN
                        inB = .TRUE.
                        EXIT
                    END IF
                END DO
                IF (.NOT. inB) THEN
                    count = count + 1
                END IF
            END DO

            ! ����res_���ڴ�ռ�
            ALLOCATE(res_(count))

            ! ���res_
            count = 0
            DO i = 1, SIZE(A_)
                inB = .FALSE.
                DO j = 1, SIZE(B_)
                    IF (A_(i) == B_(j)) THEN
                        inB = .TRUE.
                        EXIT
                    END IF
                END DO
                IF (.NOT. inB) THEN
                    count = count + 1
                    res_(count) = A_(i)
                END IF
            END DO

            RETURN
        END FUNCTION SET_DIFF
        
        
        
        !---------------------------------------------------------------
        !
        !  debug
        !
        !---------------------------------------------------------------
        SUBROUTINE getHullMeshesVertex_post_debug(Hull_, info_)
            IMPLICIT NONE
            TYPE(Face), ALLOCATABLE, TARGET, INTENT(IN) :: Hull_(:) !͹��
            INTEGER*4, INTENT(OUT) :: info_ ! 
            INTEGER*4 :: i, j, k, C, istat
            REAL*8, ALLOCATABLE :: meshes_(:,:,:)
            ! ������Ϣ
            info_ = 0
            
            IF(ALLOCATED(meshes_)) DEALLOCATE(meshes_)
            ALLOCATE( meshes_( SIZE(Hull_), 3, 3 ), STAT = istat )
            meshes_ = 0.D0
            DO i = 1, SIZE(Hull_)
                meshes_(i, :, :) = Hull_(i)%vertex(:,:)
            END DO
            
            
            ! �������������Ƿ�淶
            IF ( (SIZE(meshes_, 2) /= 3) .OR. (SIZE(meshes_, 3) /= 3) ) THEN
                info_ = 1 ! 1 - ������������
                RETURN
            END IF
            
            ! post ����
            C = 0
            
            ! �������ϵ
            WRITE(UNIT=6, FMT="(6I5,\)") 1, 2, 1, 3, 1, 4 ! 
            
            C = 4
            DO i = 1, SIZE(meshes_, 1), 1
                WRITE(UNIT=6, FMT="(6I5,\)") C + 1, C + 2, C + 2, C + 3, C + 3, C + 1 ! 
                C = C + 3
            END DO
            
            WRITE(UNIT=6, FMT="(A)")  ! 
            
            ! �����������һ��
            
            ! �������ϵ��s
            WRITE(UNIT=6, FMT="(12I5,\)") 0, 0, 0,  5, 0, 0,  0, 5, 0  , 0, 0 ,5! 
            
            DO i = 1, SIZE(meshes_, 1), 1
            DO j = 1, 3, 1
                WRITE(UNIT=6, FMT="(3F12.5,\)") meshes_(i,j,:) ! 
            END DO
            END DO
            
            RETURN
        END SUBROUTINE getHullMeshesVertex_post_debug

        
        !---------------------------------------------------------------
        PURE FUNCTION GET_RANDOM_UNIT_VECTOR(index_) RESULT(res_)
            IMPLICIT NONE
            INTEGER*4, INTENT(IN) :: index_ ! 
            REAL*8 :: res_(3)
            REAL*8, PARAMETER :: dataLib(3,100) = &
            RESHAPE(&
             [ 0.000001109357820885D0,  0.072093544214837393D0,  0.997397874913172555D0, &
               0.266483497218669374D0, -0.727347325988231153D0,  0.632417910157418883D0, &
               0.079214616132658941D0, -0.782543920607548071D0, -0.617535470164364719D0, &
              -0.993301267605208316D0,  0.106810772229378015D0,  0.044091390425458579D0, &
               0.082261341377368513D0,  0.991595302008176138D0, -0.099859044408155587D0, &
              -0.787452696781838490D0,  0.616178410256023601D0,  0.015569748404171571D0, &
              -0.247966562512464128D0,  0.750010049461640738D0, -0.613186357955148420D0, &
              -0.715817591888975313D0,  0.423804523888427931D0,  0.554972882827594716D0, &
               0.499764308041154848D0,  0.237809719054367125D0, -0.832875845448425078D0, &
               0.360748686617363812D0,  0.307777557994801998D0,  0.880416583157429655D0, &
               0.713138609686784886D0, -0.678418744074228530D0,  0.176582363396647901D0, &
               0.881992030996567422D0,  0.026379550968972942D0, -0.470525426039045791D0, &
              -0.267765386517834436D0,  0.464539693453386748D0, -0.844099858422679872D0, &
               0.513202226307113540D0,  0.794177664474205347D0,  0.325430963744568147D0, &
               0.266257765457365569D0,  0.689919118649417573D0,  0.673140707471819200D0, &
              -0.533214734590422568D0,  0.393416539739102400D0,  0.748936227642498564D0, &
              -0.623072641479377243D0, -0.654446770357797636D0,  0.428345547669355065D0, &
               0.584825748689458469D0,  0.437231667603634577D0,  0.683232985528625658D0, &
              -0.556342156780530561D0, -0.693940941632379182D0, -0.457087928209829908D0, &
               0.797251122953163582D0, -0.186816815361580540D0, -0.574012303394340728D0, &
               0.652717880922520921D0,  0.670487884243165855D0, -0.352711447230079078D0, &
              -0.119569576931363289D0, -0.933186657472575787D0, -0.338918542702544345D0, &
               0.662896092871913201D0, -0.734670864402726664D0,  0.144317327625279795D0, &
              -0.453865743569666802D0,  0.555714019359183631D0,  0.696554244478931106D0, &
               0.654083844194692787D0, -0.209153829113278511D0,  0.726931221320659904D0, &
               0.590510679076412859D0,  0.337909209878702432D0,  0.732880961531860775D0, &
               0.968625410428645917D0, -0.064469344047131227D0, -0.240017745073296679D0, &
              -0.836672384182689188D0, -0.337478629755403936D0,  0.431378599381644634D0, &
               0.415710848734430150D0,  0.722574771381445879D0, -0.552331594250728086D0, &
              -0.333326475889782536D0,  0.815058361243497620D0, -0.473891684077661635D0, &
              -0.652533192903382075D0, -0.591467557663984178D0,  0.473673474442383280D0, &
               0.394659527294562162D0, -0.550384256978558417D0, -0.735745218935055623D0, &
              -0.636304506189762753D0,  0.473703705794754570D0,  0.608868930492367122D0, &
              -0.719230459123433086D0, -0.158162890699728137D0,  0.676529413015133918D0, &
               0.629759138526492901D0, -0.491788561913722666D0,  0.601288148738358452D0, &
               0.584411917965700356D0, -0.367877772023600003D0,  0.723276333769192092D0, &
               0.870106618562407896D0, -0.204182999880998167D0,  0.448579730809907151D0, &
               0.529356795812083503D0, -0.718211329438827373D0,  0.451612520855297239D0, &
               0.733690094242977708D0, -0.622391387307088984D0, -0.272631264926984196D0, &
              -0.605777076602946218D0, -0.315061533953294726D0,  0.730595896022818714D0, &
              -0.761009425976650333D0, -0.636619547995314727D0,  0.124820690131605891D0, &
              -0.646761961270369112D0, -0.761942845893679443D0, -0.033794452875378959D0, &
               0.365154502536077674D0,  0.505749055061637143D0, -0.781588179658502025D0, &
               0.574247267419540908D0,  0.634851826576257938D0,  0.516917047652695638D0, &
               0.346341716472641781D0, -0.550932683186623917D0, -0.759289532410115098D0, &
              -0.794685184986554050D0, -0.055389826957407198D0,  0.604489391000797349D0, &
              -0.416259521322270454D0, -0.054995592820233065D0, -0.907581123469910711D0, &
               0.794777927582307919D0,  0.342095783921817331D0, -0.501296838660377997D0, &
              -0.338337965454608924D0, -0.286035970801144568D0, -0.896499216140138389D0, &
              -0.726532004741409887D0, -0.049688151104356579D0, -0.685333738937649595D0, &
              -0.603734615736470803D0, -0.585014438414317439D0,  0.541537275363678683D0, &
              -0.676560375498003186D0, -0.722348934167962309D0,  0.143101626868494480D0, &
               0.586582880385051575D0,  0.072766280975167824D0, -0.806613657702508258D0, &
              -0.755532705527683479D0, -0.071266043707085253D0, -0.651223066155029895D0, &
              -0.920701606636518566D0,  0.311540070620156373D0,  0.235056027225258340D0, &
               0.541712171882508864D0, -0.838526306892959261D0,  0.058494063654270075D0, &
              -0.408115455093796653D0, -0.092597310866135374D0, -0.908222171791651101D0, &
              -0.258240219479359101D0, -0.908622337155473581D0,  0.328203347736395479D0, &
              -0.061612129227968819D0, -0.446992987857170232D0,  0.892413141061087156D0, &
               0.788042672316281223D0, -0.496244917147545261D0,  0.364320914598434853D0, &
              -0.248619129130190686D0,  0.619445212796131295D0, -0.744631557869058658D0, &
               0.727207891810358387D0, -0.392604991169558049D0, -0.563054174123134521D0, &
              -0.730052156895783066D0,  0.157234865285751285D0,  0.665057174497340808D0, &
               0.600414670664006778D0,  0.750265884008508910D0,  0.276773059643389052D0, &
              -0.083928500830154310D0,  0.690568080639724524D0,  0.718381328230327632D0, &
               0.694831042024156353D0,  0.584804220606428005D0, -0.418585530806468986D0, &
              -0.111848450943919986D0, -0.781531383436509852D0, -0.613757786692161189D0, &
              -0.279182094755242194D0, -0.930461735000781665D0, -0.237272665234346397D0, &
              -0.689964963785805074D0, -0.305025070889099192D0, -0.656435872631251471D0, &
               0.633382581384791088D0,  0.583236672149216373D0,  0.508587740570539015D0, &
               0.466924244038473768D0, -0.606103736912413371D0,  0.643909939688702027D0, &
              -0.137658227056735444D0, -0.193627586092586290D0, -0.971369430457616478D0, &
               0.393853240338342958D0,  0.768953844741995574D0,  0.503576816117948800D0, &
              -0.132535470218959284D0,  0.729368436809752718D0, -0.671160213748950629D0, &
               0.159029880166712406D0,  0.267247506574191773D0,  0.950414787050390064D0, &
               0.585440601303706010D0, -0.650059126571057910D0,  0.484440331007677694D0, &
               0.086766095195569742D0, -0.926700911609081412D0,  0.365646092755564367D0, &
               0.404761320436991479D0, -0.409969869053845359D0, -0.817369549191842681D0, &
              -0.630382450683336315D0,  0.770188809015893039D0, -0.097093585458315548D0, &
              -0.042053492941287379D0, -0.611271645428856480D0, -0.790302776931813389D0, &
               0.929725661108754209D0,  0.077330619173836948D0, -0.360041900914436386D0, &
              -0.889604251783720934D0, -0.344981229410519730D0, -0.299319606044663511D0, &
               0.129702915764274479D0, -0.696106796017660678D0, -0.706124976318124986D0, &
              -0.796994723739967381D0, -0.420325416758673909D0, -0.433734889485847597D0, &
              -0.643021987392653815D0, -0.525087908251825164D0,  0.557499248732520325D0, &
               0.223259530927500754D0, -0.439307839166757808D0,  0.870151598456651798D0, &
               0.639217882809690274D0,  0.671377686488942249D0,  0.375036665382270096D0, &
               0.228323372420344811D0, -0.748223967023273318D0, -0.622920005119883879D0, &
              -0.632452534964462632D0,  0.397443937197173747D0, -0.664862472848508856D0, &
              -0.575267651846246730D0,  0.586755089131675400D0,  0.569899635126559057D0, &
               0.934572561750450670D0,  0.355419405776895792D0,  0.015848432742659273D0, &
              -0.122211293462219608D0,  0.261591882966958789D0,  0.957410093176425669D0, &
               0.418206651287156450D0, -0.714638510825073237D0,  0.560709368269252773D0, &
              -0.455037020713617735D0,  0.389115382040291002D0,  0.800956009553404180D0, &
               0.576937065595787169D0, -0.543479726634975457D0,  0.609732243758270287D0, &
              -0.094516770591717383D0,  0.753943490941892613D0,  0.650104447410771891D0, &
               0.489068888565033721D0, -0.424755340422356520D0,  0.761836283607213560D0, &
               0.986861350764715373D0,  0.139794765568494128D0,  0.081006776793618909D0, &
              -0.902962972513389861D0, -0.262938852206923646D0,  0.339883848203895222D0, &
              -0.712980642840275625D0,  0.087812143183863101D0,  0.695663446247195227D0], &
            [3,100])
                
            res_ = dataLib(:, index_)
            
            RETURN
        END FUNCTION GET_RANDOM_UNIT_VECTOR
        
        
    END MODULE GCLIB_QuickHull

            
            
            
            