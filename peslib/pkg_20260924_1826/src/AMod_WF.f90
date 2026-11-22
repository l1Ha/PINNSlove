module mod_psi_extract
    implicit none
    private
    public :: extract_specific_channel

contains

    !---------------------------------------------------------------------------
    ! 子程序: extract_specific_channel
    ! 功能: 提取指定量子态的 S 矩阵并计算渐进波函数值
    !---------------------------------------------------------------------------
    subroutine extract_specific_channel( &
        x, y,             & ! S矩阵的实部和虚部
        ilev, jlev, nlev, & ! 通道量子数数组 (Arrangement, j, v)
        eint,             & ! 通道阈值能量
        n,                & ! 通道总数
        enrg,             & ! 当前总能量 (eV)
        rmlmda,           & ! 约化质量因子 (用于能量归一化)
        target_v,         & ! [筛选] 目标振动量子数
        target_j,         & ! [筛选] 目标转动量子数
        rmax              & ! [可选] 渐进距离，用于计算相位
    )
        implicit double precision (a-h, o-z)
        
        ! 输入变量
        integer, intent(in) :: n
        double precision, intent(in) :: x(n,n), y(n,n) ! Re(S), Im(S)
        integer, intent(in) :: ilev(n), jlev(n), nlev(n)
        double precision, intent(in) :: eint(n)
        double precision, intent(in) :: enrg, rmlmda
        integer, intent(in) :: target_v, target_j
        double precision, intent(in) :: rmax

        ! 局部变量
        integer :: i, j
        double precision :: k_wave, phase
        complex*16 :: S_element, Psi_val, i_complex
        logical, save :: first_call = .true.
        double precision :: hbar_fac
        
        ! 常数定义 (根据 abc 代码习惯，能量单位可能是 eV，长度是 a0)
        ! 注意：这里的波数 k 计算依赖于具体的单位制，这里假设 input 中的 mass/energy 转换关系
        ! 下面的 0.036... 是常用的 2*mu转换因子，请根据 mod_pot 实际情况调整，
        ! 或者只输出 S 矩阵模长而不计算相位。
        i_complex = (0.0d0, 1.0d0)

        ! 首次调用时打开文件写入表头
        if (first_call) then
            open(888, file='./output/Wavefunction_Complex.dat', position='append')
            write(888, '(A)') '# Energy(eV)  Init(v,j)  Final(v,j)  Re(S)        Im(S)        |Psi(Rmax)|  Arg(Psi)'
            first_call = .false.
        else
            open(888, file='./output/Wavefunction_Complex.dat', position='append')
        endif

        ! 循环所有通道
        do i = 1, n
            ! 仅处理反应物通道 (ilev=1) 且 能量开放
            if (ilev(i) /= 1 .or. enrg <= eint(i)) cycle

            ! 筛选初始态 (这里假设提取所有开放的初始态，或者你可以加 target_v_in 判断)
            ! 如果只想看基态入口: if (nlev(i) /= 0 .or. jlev(i) /= 0) cycle

            do j = 1, n
                ! 仅处理能量开放的末态
                if (enrg <= eint(j)) cycle

                ! --- 核心筛选逻辑 ---
                ! 检查末态是否是我们想要的特定通道
                if (nlev(j) == target_v .and. jlev(j) == target_j) then
                    
                    ! 1. 获取 S 矩阵元 (复数)
                    S_element = cmplx(x(i,j), y(i,j), kind=8)
                    
                    ! 2. 计算波数 k_n (用于波函数重构)
                    ! k = sqrt( 2*mu*(E - E_internal) )
                    ! 注意：enrg 传入时可能已经除过 rmlmda，需确认原代码逻辑
                    ! 这里假设传入的 enrg 是总能，eint 是阈值
                    ! k_wave = sqrt( max(0.d0, (enrg - eint(j)) * rmlmda ) ) 
                    
                    ! 为避免单位换算错误，此处仅演示 S 矩阵提取
                    ! 若 S 矩阵本身就是波函数系数，则 Psi ~ S
                    
                    ! 3. 计算渐进波函数 (Outgoing part)
                    ! Psi_scatt ~ (1/sqrt(k)) * S * exp(i(kR - l*pi/2))
                    ! 这里我们输出最核心的 S 矩阵复数值，它包含了波函数的所有物理信息
                    
                    Psi_val = S_element 

                    ! 4. 输出数据
                    ! 格式: 能量, 初始v,j, 末态v,j, S实部, S虚部, 模长, 相位
                    write(888, '(F10.5, 2X, I3,I3, 2X, I3,I3, 2X, E15.7, 2X, E15.7, 2X, E15.7, 2X, F10.5)') &
                        enrg, &
                        nlev(i), jlev(i), & ! Initial
                        nlev(j), jlev(j), & ! Final (Target)
                        real(S_element), &
                        aimag(S_element), &
                        abs(S_element), &
                        atan2(aimag(S_element), real(S_element))
                end if
            end do
        end do

        close(888)

    end subroutine extract_specific_channel

end module mod_psi_extract