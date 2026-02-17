
module ConstitutiveModels_mod

    use PHYSICAL_MODULE, only: Xparam
    use Tools_mod, only: traceTensor, getParameterID

    implicit none

    contains

    subroutine VE_model_LPTT(C_Tensor, f)

        implicit none

        real(8), dimension(:,:), intent(in) :: C_Tensor
        real(8), dimension(3), intent(out) :: f

        real(8) :: trC, EpN

        EpN = Xparam(getParameterID('EpN'))

        trC = traceTensor(C_Tensor)
        f(1) = 1.0d0 + (EpN*(trC-3.0d0))
        f(2) = f(1)
        f(3) = 0.0d0

    end subroutine VE_model_LPTT

    !-----------------------------------------------------------------------

    subroutine VE_model_ePTT(C_Tensor, f)

        implicit none

        real(8), dimension(:,:), intent(in) :: C_Tensor
        real(8), dimension(3), intent(out) :: f

        real(8) :: trC, EpN

        EpN = Xparam(getParameterID('EpN'))

        trC = traceTensor(C_Tensor)
        f(1) = exp(EpN*(trC-3.0d0))
        f(2) = f(1)
        f(3) = 0.0d0

    end subroutine VE_model_ePTT

    !-----------------------------------------------------------------------

    subroutine VE_model_Giesekus(f)

        implicit none

        real(8), dimension(3), intent(out) :: f
        
        real(8) :: EpN

        EpN = Xparam(getParameterID('EpN'))

        f(1) = 1.0d0-EpN
        f(2) = 1.0d0-2.0d0*EpN
        f(3) = EpN

    end subroutine VE_model_Giesekus

    !-----------------------------------------------------------------------

    subroutine VE_model_FENE_CR(C_Tensor, f)

        implicit none

        real(8), dimension(:,:), intent(in) :: C_Tensor
        real(8), dimension(3), intent(out) :: f

        real(8) :: trC, EpN

        EpN = Xparam(getParameterID('EpN'))

        trC = traceTensor(C_Tensor)
        f(1) = (EpN-3.0d0)/(EpN-trC)
        f(2) = f(1)
        f(3) = 0.0d0

    end subroutine VE_model_FENE_CR

    !-----------------------------------------------------------------------

    subroutine VE_model_FENE_P(C_Tensor, f)

        implicit none

        real(8), dimension(:,:), intent(in) :: C_Tensor
        real(8), dimension(3), intent(out) :: f

        real(8) :: trC, EpN

        EpN = Xparam(getParameterID('EpN'))

        trC = traceTensor(C_Tensor)
        f(1) = 1.0d0
        f(2) = (EpN-3.0d0)/(EpN-trC)
        f(3) = 0.0d0

    end subroutine VE_model_FENE_P

    !-----------------------------------------------------------------------

    subroutine VE_tensor_regular(g)

        implicit none

        real(8), dimension(2), intent(out) :: g

        g(1) = 1.0d0
        g(2) = 1.0d0

    end subroutine VE_tensor_regular

    !-----------------------------------------------------------------------

    subroutine VE_tensor_FENE_CR(C_Tensor, g)

        implicit none

        real(8), dimension(:,:), intent(in) :: C_Tensor
        real(8), dimension(2), intent(out) :: g

        real(8) :: trC, EpN

        EpN = Xparam(getParameterID('EpN'))

        trC = traceTensor(C_Tensor)
        g(1) = (EpN-3.0d0)/(EpN-trC)
        g(2) = g(1)

    end subroutine VE_tensor_FENE_CR

    !-----------------------------------------------------------------------

    subroutine VE_tensor_FENE_P(C_Tensor, g)

        implicit none

        real(8), dimension(:,:), intent(in) :: C_Tensor
        real(8), dimension(2), intent(out) :: g

        real(8) :: trC, EpN

        EpN = Xparam(getParameterID('EpN'))

        trC = traceTensor(C_Tensor)
        g(1) = 1.0d0
        g(2) = (EpN-3.0d0)/(EpN-trC)

    end subroutine VE_tensor_FENE_P

    !-----------------------------------------------------------------------

    subroutine Elastic_tensor_Neo_Hookean(Ii, dwi)

        implicit none

        real(8), dimension(3), intent(in) :: Ii
        real(8), dimension(3), intent(out) :: dwi

        dwi(1) = 1.0d0
        dwi(2) = 0.0d0
        dwi(3) = 0.0d0

    end subroutine Elastic_tensor_Neo_Hookean

    !-----------------------------------------------------------------------

    subroutine Elastic_tensor_Mooney(Ii, dwi)

        implicit none

        real(8), dimension(3), intent(in) :: Ii
        real(8), dimension(3), intent(out) :: dwi

        real(8) :: BetaN

        BetaN = Xparam(getParameterID('BetaN'))

        dwi(1) = 1.0d0
        dwi(2) = BetaN
        dwi(3) = 0.0d0

    end subroutine Elastic_tensor_Mooney

    !-----------------------------------------------------------------------

    subroutine Elastic_tensor_Gent_Thomas(Ii, dwi)

        implicit none

        real(8), dimension(3), intent(in) :: Ii
        real(8), dimension(3), intent(out) :: dwi

        real(8) :: BetaN

        BetaN = Xparam(getParameterID('BetaN'))

        dwi(1) = 1.0d0
        dwi(2) = BetaN/Ii(2)
        dwi(3) = 0.0d0

    end subroutine Elastic_tensor_Gent_Thomas

    !-----------------------------------------------------------------------

    subroutine Elastic_tensor_Gent_I(Ii, dwi)

        implicit none

        real(8), dimension(3), intent(in) :: Ii
        real(8), dimension(3), intent(out) :: dwi

        real(8) :: JmN

        JmN = Xparam(getParameterID('JmN'))

        dwi(1) = JmN/(JmN-(Ii(1)-3.0d0))
        dwi(2) = 0.0d0
        dwi(3) = 0.0d0

    end subroutine Elastic_tensor_Gent_I

    !-----------------------------------------------------------------------

    subroutine Elastic_tensor_Gent_II(Ii, dwi)

        implicit none

        real(8), dimension(3), intent(in) :: Ii
        real(8), dimension(3), intent(out) :: dwi

        real(8) :: JmN, BetaN2

        JmN = Xparam(getParameterID('JmN'))
        BetaN2 = Xparam(getParameterID('BetaN2'))

        dwi(1) = JmN/(JmN-(Ii(1)-3.0d0))
        dwi(2) = BetaN2/Ii(2)
        dwi(3) = 0.0d0

    end subroutine Elastic_tensor_Gent_II

end module ConstitutiveModels_mod
